#include "LReceiver.h"

#include <array>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <unordered_map>
#include <iomanip>
#include <bitset>
#include <thread>

namespace scuPSI {

	void scuPSI::LReceiver::batchOutput(PRNG& prng, span<Channel> chls, block& commonSeed, const u64& senderSize, const u64& receiverSize, const u64& width, span<block> receiverSet, u64& hashLengthInBytes)
	{
		//------------------ Base OTs ----------------------
		//PSI接收方作为OT的发送方
		std::cout << "Receiver:base OT start" << std::endl;
		timer.setTimePoint("Sender start");
		IknpOtExtSender otExtSender;
		//std::cout << ch <<std::endl;
		std::cout << "Receiver:genBaseOts(prng, ch)前" << std::endl;
		otExtSender.genBaseOts(prng, chls[0]);
		std::cout << "Receiver:genBaseOts(prng, ch)后" << std::endl;
		std::vector<std::array<block, 2> > otMessages(width);//定义两个ot信息
		otExtSender.send(otMessages, prng, chls[0]);
		std::cout << "Receiver:base OT finished\n";
		timer.setTimePoint("Receiver base OT finished");

		//---------------- 哈希算法 并行 --------------------
		itemsToBins(chls, commonSeed, receiverSet, receiverSize, width, hashLengthInBytes, otMessages);
	}


	void scuPSI::LReceiver::itemsToBins(span<Channel> chls, block& commonSeed, span<block> receiverSet, const u64& receiverSize, const u64& width, u64& hashLengthInBytes, std::vector<std::array<block, 2> >& otMessages)
	{
		u64 numThreads(chls.size());
		const bool isMultiThreaded = numThreads > 1;

		//--------均衡哈希---------
		balance.init(receiverSize, maxBinSize, recvNumDummies);
		balance.insertItems(receiverSet);
		std::cout << "Receiver:balanced binning finished" << std::endl;
		timer.setTimePoint("Receiver balanced binning finished");

		//--------线程与Bin--------
		u64 height = 64;
		u64 logHeight = 8;

		u64* psiInBin = new u64[balance.mNumBins];
		u64 psiTotal;

		u64 sentData[20] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		u64 recvData[20] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		u64 totalSentData, totalRecvData, totalData;


		//--------线程与Bin--------
		auto routine = [&](u64 t)
		{
			auto& chl = chls[t];
			u64 binStartIdx = balance.mNumBins * t / numThreads;
			u64 tempBinEndIdx = (balance.mNumBins * (t + 1) / numThreads);
			u64 binEndIdx = std::min(tempBinEndIdx, balance.mNumBins);

			for (u64 i = binStartIdx; i < binEndIdx; i += stepSize)
			{
				auto curStepSize = std::min(stepSize, binEndIdx - i);

				u64 iterSend = 0;

				for (u64 k = 0; k < curStepSize; ++k)
				{
					u64 bIdx = i + k;
					u64 binCapacity = balance.mBins[bIdx].blks.size();
					span<block> itemsInBin = balance.mBins[bIdx].blks;

					//+----------------------------------- 计算单个Bin的交集 --------------------------------+
					std::cout << "Receiver:单次调用()函数" << std::endl;
					run(chl, commonSeed, binCapacity, receiverSize, height, logHeight, itemsInBin, bIdx, psiInBin[bIdx], width, hashLengthInBytes, otMessages);
				}
			}

			for (u64 i = 0; i < t; i++) {
				sentData[i] = chl.getTotalDataSent();
				recvData[i] = chl.getTotalDataRecv();
			}
		};

		//+-------------- 多线程 ----------------+
		std::vector<std::thread> thrds(chls.size());
		for (u64 i = 0; i < thrds.size(); ++i)
		{
			thrds[i] = std::thread([=] {
				routine(i);
				});
		}

		for (auto& thrd : thrds)
			thrd.join();

		//+--------- 验证结果是否正确 ---------+
		for (u64 i = 0; i < balance.mNumBins; i++) {
			psiTotal += psiInBin[i];
		}
		if (psiTotal == 100) {
			std::cout << "Receiver intersection computed - correct!\n";
		}
		timer.setTimePoint("Receiver intersection computed");
		std::cout << timer;

		for (u64 i = 0; i < 20; i++) {
			totalSentData += sentData[i];
			totalRecvData += recvData[i];
		}
		totalData = totalSentData + totalRecvData;

		std::cout << "Receiver sent communication: " << totalSentData / std::pow(2.0, 20) << " MB\n";
		std::cout << "Receiver received communication: " << totalRecvData / std::pow(2.0, 20) << " MB\n";
		std::cout << "Receiver total communication: " << totalData / std::pow(2.0, 20) << " MB\n";
	}


	void scuPSI::LReceiver::run(Channel& ch, block& commonSeed, const u64& senderSize, const u64& receiverSize, const u64& height, const u64& logHeight, span<block> receiverSet, u64 bIdx, u64& numsInBin, const u64& width, u64& hashLengthInBytes, std::vector<std::array<block, 2> >& otMessages)
	{
		const u64 bucket1 = 1 << 8; //(256)
		const u64 bucket2 = 1 << 8;
		auto heightInBytes = (height + 7) / 8;
		auto widthInBytes = (width + 7) / 8;
		auto locationInBytes = (logHeight + 7) / 8;
		auto receiverSizeInBytes = (receiverSize + 7) / 8;
		auto shift = (1 << logHeight) - 1;
		auto widthBucket1 = sizeof(block) / locationInBytes;

		//////////// Initialization ///////////////////

		PRNG commonPrng(commonSeed);
		block commonKey;
		AES commonAes;

		u8* matrixA[widthBucket1];//widthBucket1个列指针
		u8* matrixDelta[widthBucket1];
		for (auto i = 0; i < widthBucket1; ++i) {
			matrixA[i] = new u8[heightInBytes];//矩阵A_每个列指针指向一个长度为heightInBytes的数组
			matrixDelta[i] = new u8[heightInBytes];//矩阵D
		}

		u8* transLocations[widthBucket1];
		for (auto i = 0; i < widthBucket1; ++i) {
			transLocations[i] = new u8[receiverSize * locationInBytes + sizeof(u32)];//为什么要加sizeof(u32)？
		}

		block randomLocations[bucket1];

		u8* transHashInputs[width];
		for (auto i = 0; i < width; ++i) {
			transHashInputs[i] = new u8[receiverSizeInBytes];
			memset(transHashInputs[i], 0, receiverSizeInBytes);
		}

		// std::cout << "Receiver initialized\n";
		// timer.setTimePoint("Receiver initialized");

		/////////// Transform input /////////////////////

		commonPrng.get((u8*)&commonKey, sizeof(block));
		commonAes.setKey(commonKey);

		block* recvSet = new block[receiverSize];
		block* aesInput = new block[receiverSize];
		block* aesOutput = new block[receiverSize];
		const u64 h1LengthInBytes = 32;

		RandomOracle H1(h1LengthInBytes);
		u8 h1Output[h1LengthInBytes];

		for (auto i = 0; i < receiverSize; ++i) {//1. 随机预言机，先使用blake2函数哈希
			H1.Reset();
			H1.Update((u8*)(receiverSet.data() + i), sizeof(block));
			H1.Final(h1Output);

			aesInput[i] = *(block*)h1Output;
			recvSet[i] = *(block*)(h1Output + sizeof(block));
		}

		commonAes.ecbEncBlocks(aesInput, receiverSize, aesOutput);//2. 然后使用AES加密（PRG）

		for (auto i = 0; i < receiverSize; ++i) {//3. 最后recvSet与AES输出相异或
			recvSet[i] ^= aesOutput[i];
		}

		// std::cout << "Receiver set transformed\n";
		// timer.setTimePoint("Receiver set transformed");

		for (auto wLeft = 0; wLeft < width; wLeft += widthBucket1) {

			auto wRight = wLeft + widthBucket1 < width ? wLeft + widthBucket1 : width;//1. 横向，取较小值，处理边界情况
			auto w = wRight - wLeft;

			//////////// Compute random locations (transposed) ////////////////

			commonPrng.get((u8*)&commonKey, sizeof(block));
			commonAes.setKey(commonKey);

			for (auto low = 0; low < receiverSize; low += bucket1) {

				auto up = low + bucket1 < receiverSize ? low + bucket1 : receiverSize;//2. 纵向，取较小值，处理边界情况

				commonAes.ecbEncBlocks(recvSet + low, up - low, randomLocations);//3. Aes加密，生成randomLocations

				for (auto i = 0; i < w; ++i) {
					for (auto j = low; j < up; ++j) {
						memcpy(transLocations[i] + j * locationInBytes, (u8*)(randomLocations + (j - low)) + i * locationInBytes, locationInBytes);//4. 矩阵转置randomLocations -> transLocations ? 
					}
				}
			}

			//////////// Compute matrix Delta /////////////////////////////////

			for (auto i = 0; i < widthBucket1; ++i) {
				memset(matrixDelta[i], 255, heightInBytes);//使用全1填充矩阵Delta
			}

			for (auto i = 0; i < w; ++i) {
				for (auto j = 0; j < receiverSize; ++j) {
					auto location = (*(u32*)(transLocations[i] + j * locationInBytes)) & shift;

					matrixDelta[i][location >> 3] &= ~(1 << (location & 7));//非映射位置取反，然后按位与操作
				}
			}

			//////////////// Compute matrix A & sent matrix //////////////////仍发送较大的扩展矩阵?-与原始的概念有所不同

			u8* sentMatrix[w];

			for (auto i = 0; i < w; ++i) {
				PRNG prng(otMessages[i + wLeft][0]);//otMessages0所对应的matrixA0
				prng.get(matrixA[i], heightInBytes);

				sentMatrix[i] = new u8[heightInBytes];//otMessages1所对应的matrixA1
				prng.SetSeed(otMessages[i + wLeft][1]);
				prng.get(sentMatrix[i], heightInBytes);

				for (auto j = 0; j < heightInBytes; ++j) {
					sentMatrix[i][j] ^= matrixA[i][j] ^ matrixDelta[i][j];//sentMatrix参与异或运算，并作为最终的发送矩阵
				}

				ch.asyncSend(sentMatrix[i], heightInBytes);//按列发送
			}

			///////////////// Compute hash inputs (transposed) /////////////////////

			for (auto i = 0; i < w; ++i) {
				for (auto j = 0; j < receiverSize; ++j) {
					auto location = (*(u32*)(transLocations[i] + j * locationInBytes)) & shift;

					transHashInputs[i + wLeft][j >> 3] |= (u8)((bool)(matrixA[i][location >> 3] & (1 << (location & 7)))) << (j & 7);
				}
			}
		}

		// std::cout << "Receiver matrix sent and transposed hash input computed\n";
		// timer.setTimePoint("Receiver matrix sent and transposed hash input computed");

		/////////////////// Compute hash outputs ///////////////////////////

		RandomOracle H(hashLengthInBytes);
		u8 hashOutput[sizeof(block)];
		std::unordered_map<u64, std::vector<std::pair<block, u32>>> allHashes;
		u8* hashInputs[bucket2];
		for (auto i = 0; i < bucket2; ++i) {
			hashInputs[i] = new u8[widthInBytes];
		}

		for (auto low = 0; low < receiverSize; low += bucket2) {
			auto up = low + bucket2 < receiverSize ? low + bucket2 : receiverSize;

			for (auto j = low; j < up; ++j) {
				memset(hashInputs[j - low], 0, widthInBytes);
			}

			for (auto i = 0; i < width; ++i) {
				for (auto j = low; j < up; ++j) {
					hashInputs[j - low][i >> 3] |= (u8)((bool)(transHashInputs[i][j >> 3] & (1 << (j & 7)))) << (i & 7);
				}
			}

			for (auto j = low; j < up; ++j) {
				H.Reset();
				H.Update(hashInputs[j - low], widthInBytes);
				H.Final(hashOutput);

				allHashes[*(u64*)hashOutput].push_back(std::make_pair(*(block*)hashOutput, j));
			}
		}

		// std::cout << "Receiver hash outputs computed\n";
		// timer.setTimePoint("Receiver hash outputs computed");

		///////////////// Receive hash outputs from sender and compute PSI ///////////////////

		u8* recvBuff = new u8[bucket2 * hashLengthInBytes];

		auto psi = 0;

		for (auto low = 0; low < senderSize; low += bucket2) {
			auto up = low + bucket2 < senderSize ? low + bucket2 : senderSize;

			ch.recv(recvBuff, (up - low) * hashLengthInBytes);

			for (auto idx = 0; idx < up - low; ++idx) {
				u64 mapIdx = *(u64*)(recvBuff + idx * hashLengthInBytes);

				auto found = allHashes.find(mapIdx);
				if (found == allHashes.end()) continue;

				for (auto i = 0; i < found->second.size(); ++i) {
					if (memcmp(&(found->second[i].first), recvBuff + idx * hashLengthInBytes, hashLengthInBytes) == 0) {
						++psi;
						break;
					}
				}
			}
		}

		/*if (psi == 100) {
			std::cout << "Receiver intersection computed - correct!\n";
		}
		timer.setTimePoint("Receiver intersection computed");

		std::cout << timer;*/

		//////////////// Output communication /////////////////

		/*u64 sentData = ch.getTotalDataSent();
		u64 recvData = ch.getTotalDataRecv();
		u64 totalData = sentData + recvData;*/

		/*std::cout << "Receiver sent communication: " << sentData / std::pow(2.0, 20) << " MB\n";
		std::cout << "Receiver received communication: " << recvData / std::pow(2.0, 20) << " MB\n";
		std::cout << "Receiver total communication: " << totalData / std::pow(2.0, 20) << " MB\n";*/
	}
}

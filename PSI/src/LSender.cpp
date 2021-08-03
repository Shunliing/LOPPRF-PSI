#include "LSender.h"

#include <iostream>
#include <cstdlib>
#include <bitset>
#include <cmath>

namespace scuPSI {

	//批处理OPPRF
	void LSender::batchOutput(PRNG& prng, span<Channel> chls, block& commonSeed, const u64& senderSize, const u64& width, span<block> senderSet, u64& hashLengthInBytes)
	{
		//////////////////// Base OTs ///////////////////////
		////std::cout << "Sender:base OT start" << std::endl;
		//timer.setTimePoint("Sender start");
		//IknpOtExtReceiver otExtReceiver;
		////std::cout << "Sender:genBaseOts(prng, ch)前" << std::endl;
		//otExtReceiver.genBaseOts(prng, chls[0]);
		////std::cout << "Sender:genBaseOts(prng, ch)后" << std::endl;
		//BitVector choices(width);
		//std::vector<block> otMessages(width);
		//prng.get(choices.data(), choices.sizeBytes());
		//otExtReceiver.receive(choices, otMessages, prng, chls[0]);

		////std::cout << "Sender:base OT finished\n";
		//timer.setTimePoint("Sender base OT finished");

		//-----------------参考SpOT实现------------------------
		std::vector<std::array<block, 2>> baseOtSend(128);
		NaorPinkas baseOTs;
		baseOTs.send(baseOtSend, prng, chls[0], chls.size());
		IknpOtExtReceiver recvIKNP;
		recvIKNP.setBaseOts(baseOtSend);
		
		BitVector otChoices;
		otChoices.resize(width);
		otChoices.randomize(prng);
		std::vector<block> otMessages(width);
		recvIKNP.receive(otChoices, otMessages, prng, chls[0]);
		std::cout << "Sender:base OT finished\n";

		//////////////// 哈希算法 并行 ///////////////////////
		itemsToBins(chls, commonSeed, senderSet, senderSize, width, hashLengthInBytes, otMessages, otChoices);
	}


	//调用简单哈希，将元素映射到Bin中
	void LSender::itemsToBins(span<Channel> chls, block& commonSeed, span<block> senderSet, const u64& senderSize, const u64& width, u64& hashLengthInBytes, std::vector<block>& otMessages, BitVector& choices)
	{
		u64 numThreads(chls.size());
		const bool isMultiThreaded = numThreads > 1;

		//======简单哈希======
		simple.init(senderSize, maxBinSize, recvNumDummies);
		simple.insertItems(senderSet);
		std::cout << "Sender:simple binning finished" << std::endl;
		timer.setTimePoint("Sender simple binning finished");

		//====线程与Bins======
		auto routine = [&](u64 t)
		{
			auto& chl = chls[t];
			u64 binStartIdx = simple.mNumBins * t / numThreads;
			u64 tempBinEndIdx = (simple.mNumBins * (t + 1) / numThreads);//计算endIdx的中间变量
			u64 binEndIdx = std::min(tempBinEndIdx, simple.mNumBins);//处理边界条件

			u64 height = 128;
			u64 logHeight = 7;

			for (u64 i = binStartIdx; i < binEndIdx; i += stepSize)// 以stepSize个Bin作为循环单位
			{
				auto curStepSize = std::min(stepSize, binEndIdx - i);

				u64 iterSend = 0, iterRecv = 0;

				for (u64 k = 0; k < curStepSize; ++k)// 每次处理一个Bin,对应矩阵Q的一行
				{
					u64 bIdx = i + k;
					u64 binCapacity = simple.mBins[bIdx].blks.size();
					span<block> itemsInBin = simple.mBins[bIdx].blks;
					//u64 recvBinCapacity = binCapacity;

					//------------------------- 计算单个Bin的oprf输出 -----------------------
					std::cout << "Sender:run() " << "Thread:" << t << "; Bin:" << i + k << std::endl;
					run(chl, commonSeed, binCapacity, height, logHeight, itemsInBin, bIdx, width, hashLengthInBytes, otMessages, choices);

					//rowQ[k].resize(simple.mBins[bIdx].blks.size());// 矩阵Q的行数 = 元素个数(使用原有协议中的OT扩展方法，并令矩阵行数大于元素个数)
					//
					////=====================Compute OT row=====================
					//prfOtRows(simple.mBins[bIdx].blks, rowQ[k], mAesQ);//得到Bin中元素对应的OT扩展矩阵Q，其中mAesQ(OT扩展矩阵，与OT选择位有关)
				}
				//std::vector<u8> recvBuff;
				//chl.recv(recvBuff);// (接收R的矩阵，经过OT的选择得到Q,在OPRF-PSI中)
			}
		};

		//+------------ 多线程 --------------+
		std::vector<std::thread> thrds(chls.size());
		for (u64 i = 0; i < thrds.size(); ++i)
		{
			thrds[i] = std::thread([=] {
				routine(i);
				});
		}

		for (auto& thrd : thrds)
			thrd.join();

		//
		std::cout << "Sender:hash outputs computed and sent" << std::endl;
		timer.setTimePoint("Sender hash outputs computed and sent");
		std::cout << timer;
	}

	//在每对Bin中执行PSI
	void LSender::run(Channel& ch, block& commonSeed, const u64& senderSize, const u64& height, const u64& logHeight, span<block> senderSet, u64 bIdx, const u64& width, u64& hashLengthInBytes, std::vector<block>& otMessages, BitVector& choices)
	{//将std::<vector>& senderSet修改为span<block> senderSet

		const u64 bucket1 = 1 << 8;//(256)
		const u64 bucket2 = 1 << 8;
		auto heightInBytes = (height + 7) / 8;
		auto widthInBytes = (width + 7) / 8;
		auto locationInBytes = (logHeight + 7) / 8;
		auto senderSizeInBytes = (senderSize + 7) / 8;
		auto shift = (1 << logHeight) - 1;
		auto widthBucket1 = sizeof(block) / locationInBytes;

		////////////// Initialization //////////////////////

		PRNG commonPrng(commonSeed);
		block commonKey;
		AES commonAes;

		u8* transLocations[widthBucket1];
		for (auto i = 0; i < widthBucket1; ++i) {
			transLocations[i] = new u8[senderSize * locationInBytes + sizeof(u32)];
		}

		block randomLocations[bucket1];

		u8* matrixC[widthBucket1];
		for (auto i = 0; i < widthBucket1; ++i) {
			matrixC[i] = new u8[heightInBytes];
		}

		u8* transHashInputs[width];
		for (auto i = 0; i < width; ++i) {
			transHashInputs[i] = new u8[senderSizeInBytes];
			memset(transHashInputs[i], 0, senderSizeInBytes);
		}

		/////////// Transform input /////////////////////

		commonPrng.get((u8*)&commonKey, sizeof(block));
		commonAes.setKey(commonKey);

		block* sendSet = new block[senderSize];
		block* aesInput = new block[senderSize];
		block* aesOutput = new block[senderSize];
		const u64 h1LengthInBytes = 32;

		RandomOracle H1(h1LengthInBytes);
		u8 h1Output[h1LengthInBytes];

		for (auto i = 0; i < senderSize; ++i) {
			H1.Reset();
			H1.Update((u8*)(senderSet.data() + i), sizeof(block));
			H1.Final(h1Output);

			aesInput[i] = *(block*)h1Output;
			sendSet[i] = *(block*)(h1Output + sizeof(block));
		}

		commonAes.ecbEncBlocks(aesInput, senderSize, aesOutput);
		for (auto i = 0; i < senderSize; ++i) {
			sendSet[i] ^= aesOutput[i];
		}

		// std::cout << "Sender set transformed\n";
		// timer.setTimePoint("Sender set transformed");

		for (auto wLeft = 0; wLeft < width; wLeft += widthBucket1) {
			auto wRight = wLeft + widthBucket1 < width ? wLeft + widthBucket1 : width;
			auto w = wRight - wLeft;

			//////////// Compute random locations (transposed) ////////////////

			commonPrng.get((u8*)&commonKey, sizeof(block));
			commonAes.setKey(commonKey);

			for (auto low = 0; low < senderSize; low += bucket1) {

				auto up = low + bucket1 < senderSize ? low + bucket1 : senderSize;

				commonAes.ecbEncBlocks(sendSet + low, up - low, randomLocations);

				for (auto i = 0; i < w; ++i) {
					for (auto j = low; j < up; ++j) {
						memcpy(transLocations[i] + j * locationInBytes, (u8*)(randomLocations + (j - low)) + i * locationInBytes, locationInBytes);
					}
				}
			}

			//////////////// Extend OTs and compute matrix C ///////////////////

			u8* recvMatrix;
			recvMatrix = new u8[heightInBytes];

			for (auto i = 0; i < w; ++i) {
				PRNG prng(otMessages[i + wLeft]);
				prng.get(matrixC[i], heightInBytes);
				std::cout << "Sender:ch.recv(recvMatrix, heightInBytes);(212)" << std::endl;
				ch.recv(recvMatrix, heightInBytes);

				if (choices[i + wLeft]) {
					for (auto j = 0; j < heightInBytes; ++j) {
						matrixC[i][j] ^= recvMatrix[j];
					}
				}
			}

			///////////////// Compute hash inputs (transposed) /////////////////////
			for (auto i = 0; i < w; ++i) {
				for (auto j = 0; j < senderSize; ++j) {
					auto location = (*(u32*)(transLocations[i] + j * locationInBytes)) & shift;

					transHashInputs[i + wLeft][j >> 3] |= (u8)((bool)(matrixC[i][location >> 3] & (1 << (location & 7)))) << (j & 7);
				}
			}
		}

		// std::cout << "Sender transposed hash input computed\n";
		// timer.setTimePoint("Sender transposed hash input computed");

		/////////////////// Compute hash outputs ///////////////////////////
		RandomOracle H(hashLengthInBytes);
		u8 hashOutput[sizeof(block)];
		u8* hashInputs[bucket2];

		for (auto i = 0; i < bucket2; ++i) {
			hashInputs[i] = new u8[widthInBytes];
		}

		for (auto low = 0; low < senderSize; low += bucket2) {
			auto up = low + bucket2 < senderSize ? low + bucket2 : senderSize;

			for (auto j = low; j < up; ++j) {
				memset(hashInputs[j - low], 0, widthInBytes);
			}

			for (auto i = 0; i < width; ++i) {
				for (auto j = low; j < up; ++j) {
					hashInputs[j - low][i >> 3] |= (u8)((bool)(transHashInputs[i][j >> 3] & (1 << (j & 7)))) << (i & 7);
				}
			}

			u8* sentBuff = new u8[(up - low) * hashLengthInBytes];

			for (auto j = low; j < up; ++j) {
				H.Reset();
				H.Update(hashInputs[j - low], widthInBytes);
				H.Final(hashOutput);

				memcpy(sentBuff + (j - low) * hashLengthInBytes, hashOutput, hashLengthInBytes);
			}
			std::cout<< "Sender:ch.asyncSend(sentBuff, (up - low) * hashLengthInBytes);(266)" << std::endl;
			ch.asyncSend(sentBuff, (up - low) * hashLengthInBytes);
		}

		// std::cout << "Sender hash outputs computed and sent\n";
		// timer.setTimePoint("Sender hash outputs computed and sent");

		std::cout << timer;
	}
}

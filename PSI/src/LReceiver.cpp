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
		////PSI接收方作为OT的发送方
		//std::cout << "Receiver:base OT start" << std::endl;
		//timer.setTimePoint("Sender start");
		//IknpOtExtSender otExtSender;
		////std::cout << ch <<std::endl;
		//std::cout << "Receiver:genBaseOts(prng, ch)前" << std::endl;
		//otExtSender.genBaseOts(prng, chls[0]);
		//std::cout << "Receiver:genBaseOts(prng, ch)后" << std::endl;
		//std::vector<std::array<block, 2> > otMessages(width);//定义两个ot
		//otExtSender.send(otMessages, prng, chls[0]);
		//std::cout << "Receiver:base OT finished\n";
		//timer.setTimePoint("Receiver base OT finished");

		//--------------------- 参考SpOT实现 -------------------------
		std::vector<block> baseOtRecv(128);
		BitVector baseOtChoices(128);
		baseOtChoices.randomize(prng);
		NaorPinkas baseOTs;
		baseOTs.receive(baseOtChoices, baseOtRecv, prng, chls[0], chls.size());
		IknpOtExtSender sendIKNP;
		sendIKNP.setBaseOts(baseOtRecv, baseOtChoices);
		std::vector<std::array<block, 2>> otMessages(width);
		sendIKNP.send(otMessages, prng, chls[0]);
		std::cout << "Receiver:base OT finished\n";
		
		//----------------------- 均衡哈希 ----------------------------
		balance.init(receiverSize, maxBinSize, recvNumDummies);
		balance.insertItems(receiverSet);
		std::cout << "Receiver:balanced binning finished" << std::endl;
		timer.setTimePoint("Receiver balanced binning finished");

		//--------------------- megaBin定义 -------------------------
		u64 binStepSize = 1336;
		u64 mNumMegaBins;
		if(balance.mNumBins % binStepSize == 0)
			mNumMegaBins = balance.mNumBins / binStepSize;
		else
			mNumMegaBins = balance.mNumBins / binStepSize + 1;
		std::vector<std::vector<block>> megaBin(mNumMegaBins);
		u64 iterMega = 0;
		
		for(u64 i = 0; i < balance.mNumBins; i += binStepSize)
		{
			u64 curBinStepSize = std::min(binStepSize, balance.mNumBins - i);
			for(u64 j = 0; j < curBinStepSize; j++)
			{
				for(u64 k = 0; k < balance.mBins[i + j].blks.size(); k++)
				{
					megaBin[iterMega].push_back(balance.mBins[i + j].blks[k]);
				}
			}
			iterMega++;
		}
		//---------------------- 哈希算法 并行 ------------------------
		//itemsToBins(chls, commonSeed, receiverSet, receiverSize, width, hashLengthInBytes, otMessages);

		//------------------- 多线程1：OPRF矩阵 ----------------------
		u64 height = 131072;
		u64 logHeight = 17;
		u64 numThreads(chls.size());
		u64 sentData[8] = {0, 0, 0, 0, 0, 0, 0, 0};
		u64 recvData[8] = {0, 0, 0, 0, 0, 0, 0, 0};
		u64 binPsi[8] = {0, 0, 0, 0, 0, 0, 0, 0};
		u64 psiTotal = 0;
		
		const bool isMultiThreaded = numThreads > 1;
		
		const u64 bucket1 = 1 << 8; //(256)
		const u64 bucket2 = 1 << 8;
		auto heightInBytes = (height + 7) / 8;
		auto widthInBytes = (width + 7) / 8;
		auto locationInBytes = (logHeight + 7) / 8;
		auto binReceiverSizeInBytes = (receiverSize + 7) / 8;
		auto shift = (1 << logHeight) - 1;
		auto widthBucket1 = sizeof(block) / locationInBytes;
		
		std::unordered_map<u64, std::vector<std::pair<block, u32>>> allHashes;
		
		u64 totalSentData, totalRecvData, totalData;
		
		//+------------ 线程定义 --------------+
		auto routine = [&](u64 t)
		{
			auto& chl = chls[t];
			u64 megaBinStartIdx = mNumMegaBins * t / numThreads;
			u64 tempMegaBinEndIdx = (mNumMegaBins * (t + 1) / numThreads);
			u64 megaBinEndIdx = std::min(tempMegaBinEndIdx, mNumMegaBins);

			for (u64 i = megaBinStartIdx; i < megaBinEndIdx; i += stepSize)
			{
				auto curStepSize = std::min(stepSize, megaBinEndIdx - i);

				u64 iterSend = 0;

				for (u64 k = 0; k < curStepSize; ++k)
				{
					u64 bIdx = i + k;
					u64 megaBinReceiverSize = megaBin[bIdx].size();
					span<block> binReceiverSet = megaBin[bIdx];//测试
					
					std::cout << "Receiver:run() " << "Thread:" << t << "; MegaBin:" << i + k << std::endl;

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
						transLocations[i] = new u8[megaBinReceiverSize * locationInBytes + sizeof(u32)];//为什么要加sizeof(u32)？
					}

					block randomLocations[bucket1];

					u8* transHashInputs[width];
					for (auto i = 0; i < width; ++i) {
						transHashInputs[i] = new u8[binReceiverSizeInBytes];
						memset(transHashInputs[i], 0, binReceiverSizeInBytes);
					}

					// std::cout << "Receiver initialized\n";
					// timer.setTimePoint("Receiver initialized");

					/////////// Transform input /////////////////////

					commonPrng.get((u8*)&commonKey, sizeof(block));
					commonAes.setKey(commonKey);

					block* recvSet = new block[megaBinReceiverSize];
					block* aesInput = new block[megaBinReceiverSize];
					block* aesOutput = new block[megaBinReceiverSize];
					const u64 h1LengthInBytes = 32;

					RandomOracle H1(h1LengthInBytes);
					u8 h1Output[h1LengthInBytes];

					for (auto i = 0; i < megaBinReceiverSize; ++i) {//1. 随机预言机，先使用blake2函数哈希
						H1.Reset();
						H1.Update((u8*)(binReceiverSet.data() + i), sizeof(block));
						H1.Final(h1Output);

						aesInput[i] = *(block*)h1Output;
						recvSet[i] = *(block*)(h1Output + sizeof(block));
					}

					commonAes.ecbEncBlocks(aesInput, megaBinReceiverSize, aesOutput);//2. 然后使用AES加密（PRG）

					for (auto i = 0; i < megaBinReceiverSize; ++i) {//3. 最后recvSet与AES输出相异或
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

						for (auto low = 0; low < megaBinReceiverSize; low += bucket1) {

							auto up = low + bucket1 < megaBinReceiverSize ? low + bucket1 : megaBinReceiverSize;//2. 纵向，取较小值，处理边界情况

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
							for (auto j = 0; j < megaBinReceiverSize; ++j) {
								auto location = (*(u32*)(transLocations[i] + j * locationInBytes)) & shift;

								matrixDelta[i][location >> 3] &= ~(1 << (location & 7));//非映射位置取反，然后按位与操作
							}
						}

						//////////////// Compute matrix A & sent matrix //////////////////

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
							//std::cout<< "Receiver:ch.asyncSend();(210)" <<  std::endl;
							//std::cout << "sentMatrix[" << i << "]" <<sentMatrix[i] << "\n";
							chl.asyncSend(std::move(sentMatrix[i]), heightInBytes);//按列发送 chl.asyncSend(std::move(sentMatrix[i]));
						}

						///////////////// Compute hash inputs (transposed) /////////////////////

						for (auto i = 0; i < w; ++i) {
							for (auto j = 0; j < megaBinReceiverSize; ++j) {
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
					
					u8* hashInputs[bucket2];
					for (auto i = 0; i < bucket2; ++i) {
						hashInputs[i] = new u8[widthInBytes];
					}

					for (auto low = 0; low < megaBinReceiverSize; low += bucket2) {
						auto up = low + bucket2 < megaBinReceiverSize ? low + bucket2 : megaBinReceiverSize;

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
					
					
					// u8* recvBuff;//u8* recvBuff = new u8[bucket2 * hashLengthInBytes],移除大小限制，并在线程外的公共区域定义

					
					
					
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
			
		//++++++++++++++++++++++++++++ 多线程2:计算PSI ++++++++++++++++++++++++++	
		auto computingPsi = [&](u64 t)
		{
			///////////////// Receive hash outputs from sender and compute PSI ///////////////////
		#if 1
			auto& chl = chls[t];
			u64 mOutputNums = 2 * senderSize / numThreads;
			u8* recvBuff = new u8[mOutputNums * hashLengthInBytes];//recvBuff指向的空间过小
			chl.recv(recvBuff, mOutputNums * hashLengthInBytes);

			// u64 startIdx = senderSize * 2 * t / numThreads;
			// u64 tmpEndIdx = senderSize * 2 * (t + 1) / numThreads;
			// u64 endIdx = std::min(tmpEndIdx, senderSize * 2);
			
			//std::cout << "测试R：是否开始计算PSI" << std::endl;
			for (auto idx = 0; idx < mOutputNums; idx++) 
			{
				// 计算PSI
				u64 mapIdx = *(u64*)(recvBuff + idx * hashLengthInBytes);//取前64位

				// 查找1：对比前64位，有可能存在多个相同的值
				auto found = allHashes.find(mapIdx);
				if (found == allHashes.end()) continue;

				// 查找2：对比所有位，在多个相同的值中查找
				for (auto i = 0; i < found->second.size(); ++i) 
				{
					if (memcmp(&(found->second[i].first), recvBuff + idx * hashLengthInBytes, hashLengthInBytes) == 0) 
					{
						++binPsi[t];//todo:BinPSI修改为线程共享
						std::cout << "binPSI: " << binPsi[t] << std::endl;
						break;
					}
				}
			}
		#endif

		#if 0
			// 在megaBin中执行元素比较
			auto& chl = chls[t];
			u64 megaBinStartIdx = mNumMegaBins * t / numThreads;
			u64 tempMegaBinEndIdx = (mNumMegaBins * (t + 1) / numThreads);
			u64 megaBinEndIdx = std::min(tempMegaBinEndIdx, balance.mNumBins);
			
			for (u64 i = megaBinStartIdx; i < megaBinEndIdx; i++)
			{
				auto curStepSize = std::min(stepSize, megaBinEndIdx - i);
				for (u64 k = 0; k < curStepSize; ++k)
				{
					

					/*if (psi == 100) {
						std::cout << "Receiver intersection computed - correct!\n";
					}
					timer.setTimePoint("Receiver intersection computed");
					std::cout << timer;*/
				}
				
				//	//////////////// Output communication /////////////////
				//	u64 sentData = ch.getTotalDataSent();
				//	u64 recvData = ch.getTotalDataRecv();
				//	u64 totalData = sentData + recvData;
				//	std::cout << "Receiver sent communication: " << sentData / std::pow(2.0, 20) << " MB\n";
				//	std::cout << "Receiver received communication: " << recvData / std::pow(2.0, 20) << " MB\n";
				//	std::cout << "Receiver total communication: " << totalData / std::pow(2.0, 20) << " MB\n";
			}
		#endif
		};
		
		//+------------ 多线程等待 --------------+
		for (u64 i = 0; i < thrds.size(); ++i)//thrds.size()
		{
			thrds[i] = std::thread([=] {
				computingPsi(i);
				});
		}

		for (auto& thrd : thrds)
			thrd.join();
			
		//++++++++++++++++++++++++++++++ 验证结果是否正确 +++++++++++++++++++++++++++++++++
		for (u64 i = 0; i < numThreads; i++) {
			psiTotal += binPsi[i];
		}
		if (psiTotal == 100) {
			std::cout << "Receiver intersection computed - correct!\n";
		}
		std::cout << "Receiver intersection: " << psiTotal << std::endl;
		timer.setTimePoint("Receiver intersection computed");
		std::cout << timer;

		for (u64 i = 0; i < 8; i++) {
			totalSentData += sentData[i];
			totalRecvData += recvData[i];
		}
		totalData = totalSentData + totalRecvData;

		std::cout << "Receiver sent communication: " << totalSentData / std::pow(2.0, 20) << " MB\n";
		std::cout << "Receiver received communication: " << totalRecvData / std::pow(2.0, 20) << " MB\n";
		std::cout << "Receiver total communication: " << totalData / std::pow(2.0, 20) << " MB\n";
		
	}
}


#include "BalancedIndex.h"
#include "cryptoTools/Crypto/PRNG.h"
#include <random>
#include "cryptoTools/Common/Log.h"
#include "cryptoTools/Common/CuckooIndex.h"
#include <numeric>
#include <fstream>

using namespace osuCrypto;
using namespace std;

namespace scuPSI
{
	void BalancedIndex::print(span<block> items)
	{
		std::cout << "numIters=" << numIters << std::endl;
		std::cout << "mNumDummies=" << mNumDummies << std::endl;
		std::cout << "mNumBins=" << mNumBins << std::endl;
		for (u64 i = 0; i < mBins.size(); ++i)
		{
			std::cout << "BBin #" << i << " contains " << mBins[i].cnt << " elements" << std::endl;

			for (auto it = mBins[i].values.begin(); it != mBins[i].values.end(); ++it)//for each bin, list all alter light bins
			{
				for (u64 j = 0; j < it->second.size(); j++)
				{
					//std::cout << "\t" << it->second[j] << "\t" << items[it->second[j]] << std::endl;
					std::cout << "\t" << items[it->second[j].mIdx] << "\t" << it->second[j].mHashIdx << std::endl;
				}
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}

	void BalancedIndex::toFile(span<block> items)
	{
		ofstream outfile;
		outfile.open("ReceiverHashing.txt");
		outfile << "mNumDummies=" << mNumDummies << endl;
		outfile << "mNumBins=" << mNumBins << endl;
		for (u64 i = 0; i < mBins.size(); ++i)
		{
			outfile << "BBin #" << i << " contains " << mBins[i].blks.size() << " elements" << endl;

			for (u64 j = 0; j < mBins[i].blks.size(); j++)
				outfile << "\t" << mBins[i].blks[j] << "\t" << mBins[i].hashIdxs[j] << endl;

			outfile << endl;
		}
		outfile << std::endl;
		outfile.close();
	}

	//初始化：参数赋值
	void BalancedIndex::init(u64 inputSize, u64 maxBinSize, u64 numDummies, u64 statSecParam)
	{
		numIters = 0;
		mMaxBinSize = maxBinSize;
		mNumBins = 1 + inputSize / (maxBinSize - numDummies);
		mNumDummies = numDummies;
		//mMaxBinSize = mNumDummies + inputSize / mNumBins;
		mHashSeed = _mm_set_epi32(4253465, 3434565, 234435, 23987025); //hardcode hash
		mAesHasher.setKey(mHashSeed);
		mBins.resize(mNumBins);
	}

	//检查：Bin是否超出了最大容量
	void BalancedIndex::check()
	{
		for (u64 idxBin = 0; idxBin < mNumBins; ++idxBin)
		{
			if (mBins[idxBin].cnt > mMaxBinSize)
				std::cout << idxBin << "\t cnt:" << mBins[idxBin].cnt << "\t" << "mMaxBinSize " << mMaxBinSize << "\n";
		}
	}

	//插入：将元素插入到哈希表中
	void BalancedIndex::insertItems(span<block> items)
	{
		u64 inputSize = items.size();
		std::vector<u64> heavyBins;

		block cipher;
		u64 b1, b2; //2 bins index

		//1st pass：循环次数为元素个数，每次处理一个元素
		for (u64 idxItem = 0; idxItem < inputSize; ++idxItem)
		{
			cipher = mAesHasher.ecbEncBlock(items[idxItem]);

			b1 = _mm_extract_epi64(cipher, 0) % mNumBins; //1st 64 bits for finding bin location
			b2 = _mm_extract_epi64(cipher, 1) % mNumBins; //2nd 64 bits for finding alter bin location

			if (mBins[b1].cnt < mBins[b2].cnt)
			{
				auto iterB2 = mBins[b1].values.find(b2); //find alter b2 in BIN b1

				if (iterB2 != mBins[b1].values.end())  //if bins[b1] has b2 as unordered_map<b2,...>, insert index of this item，在查找过程中找到了
					iterB2->second.push_back({ 0,idxItem });//为什么都是second.push_back()? -- 第二个参数列表
				else
					mBins[b1].values.emplace(std::make_pair(b2, std::vector<item>{ {0, idxItem}})); //if not, insert new map

				mBins[b1].cnt++;
			}
			else
			{
				auto iterB1 = mBins[b2].values.find(b1); //find alter b1 in BIN b2

				if (iterB1 != mBins[b2].values.end())  //if bins[b1] has b2 as unordered_map<b2,...>, insert index of this item
					iterB1->second.push_back({ 1,idxItem });
				else
					mBins[b2].values.emplace(std::make_pair(b1, std::vector<item>{ {1, idxItem}})); //if not, insert new map

				mBins[b2].cnt++;
			}
		}

		//2nd pass：find light/heavy bins after 1st pass 循环次数为Bin的个数
		for (u64 idxBin = 0; idxBin < mNumBins; ++idxBin)
		{
			//for each bin, list all alter light bins
			for (auto it = mBins[idxBin].values.begin(); it != mBins[idxBin].values.end(); ++it)
			{
				if (mBins[it->first].cnt < mMaxBinSize)
					mBins[idxBin].lightBins.emplace_back(it->first);
			}

			if (mBins[idxBin].cnt > mMaxBinSize)
				heavyBins.emplace_back(idxBin);
		}
		//std::cout << "heavyBins.size() " << heavyBins.size() << "\n";

		//3rd pass：循环次数为heavyBins的大小，对heavyBins进行调整
		//=====================Self-Balacing-Step==========================
		while (heavyBins.size() > 0)
		{
			//std::cout << numIters << "\t " << heavyBins.size() << "\t";

			u64 b1 = heavyBins[rand() % heavyBins.size()]; //choose random bin that is heavy

			//std::cout << mBins[b1].cnt << "\t";

			if (mBins[b1].lightBins.size() > 0)
			{
				u64 i2 = rand() % mBins[b1].lightBins.size(); //choose random alter bin, that is light, to balance
				u64 b2 = mBins[b1].lightBins[i2];
				//std::cout << mBins[b2].cnt;

				if (mBins[b2].cnt < mBins[b1].cnt) //if true, do the balance (double check in some unexpected cases)
				{
					auto curSubBin = mBins[b1].values.find(b2);

					u64 rB = rand() % 2;

					if (rB == 1 || mBins[b2].cnt + 1 != mBins[b1].cnt) //if not tie
					{

						if (curSubBin->second.size() == 0)
							continue;

						u64 idxBalanced = rand() % curSubBin->second.size(); //choose random item of bin b1 that can move to b2

						//Remove item
						item idxBalancedItem = curSubBin->second[idxBalanced];
						curSubBin->second.erase(curSubBin->second.begin() + idxBalanced); //remove that item from b1
						mBins[b1].cnt--;

						if (mBins[b1].cnt < mMaxBinSize) //b1 may no longer a heavy bin 
						{
							auto it = std::find(heavyBins.begin(), heavyBins.end(), b1);
							heavyBins.erase(it);
						}

						auto newSubBin = mBins[b2].values.find(b1); //place idxBalancedItem into b2

						if (newSubBin != mBins[b2].values.end()) {  //if bins[b2] has b1 as unordered_map<b1,...>, emplace it into <b1,...>
							newSubBin->second.emplace_back(idxBalancedItem);

							if (mBins[b1].cnt < mMaxBinSize)
								mBins[b2].lightBins.emplace_back(b1);
						}
						else
						{
							mBins[b2].values.emplace(std::make_pair(b1, std::vector<item>{idxBalancedItem}));
							mBins[b2].lightBins.emplace_back(b1);
						}

						mBins[b2].cnt++;

						//b2 may become an heavy bin
						if (mBins[b2].cnt > mMaxBinSize)
						{
							heavyBins.emplace_back(b2);
							mBins[b1].lightBins.erase(mBins[b1].lightBins.begin() + i2);
						}
					}
					numIters++;
				}
			}
			//std::cout << "\n";
		}
		//check();

		for (u64 idxBin = 0; idxBin < mNumBins; ++idxBin)
		{
			for (auto it = mBins[idxBin].values.begin(); it != mBins[idxBin].values.end(); ++it)//for each bin, list all alter light bins
			{
				for (u64 idx = 0; idx < it->second.size(); idx++)
				{
					if (it->second[idx].mHashIdx == 0)
						mBins[idxBin].blks.push_back(items[it->second[idx].mIdx]);
					else
						mBins[idxBin].blks.push_back(items[it->second[idx].mIdx] ^ osuCrypto::OneBlock);

					mBins[idxBin].idxs.push_back(it->second[idx].mIdx);
					mBins[idxBin].hashIdxs.push_back(it->second[idx].mHashIdx);
				}
			}
		}
	}
}

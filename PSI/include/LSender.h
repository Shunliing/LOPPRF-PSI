#pragma once
#include "LDefines.h"
#include <vector>
#include "../tools/SimpleIndex.h"

namespace scuPSI {

	class LSender {
	public:
		LSender() {}

		SimpleIndex simple;

		Timer timer;

		// std::vector<AES> mAesQ;//参考SpOT,可以删除

		void batchOutput(PRNG& prng, span<Channel> chls, block& commonSeed, const u64& senderSize, const u64& width, span<block> senderSet, u64& hashLengthInBytes);

		void itemsToBins(span<Channel> chls, block& commonSeed, span<block> senderSet, const u64& senderSize, const u64& width, u64& hashLengthInBytes, std::vector<block>& otMessages, BitVector& choices);

		void run(Channel& ch, block& commonSeed, const u64& senderSize, const u64& height, const u64& logHeight, span<block> senderSet, u64 bIdx, const u64& width, u64& hashLengthInBytes, std::vector<block>& otMessages, BitVector& choices);
	};
}

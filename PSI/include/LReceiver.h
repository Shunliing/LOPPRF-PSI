#pragma once
#include "LDefines.h"
#include <vector>
#include "../tools/BalancedIndex.h"

namespace scuPSI {
	class LReceiver {
	public:
		LReceiver() {}

		BalancedIndex balance;

		u64 mMyInputSize;

		Timer timer;

		void batchOutput(PRNG& prng, span<Channel> chls, block& commonSeed, const u64& senderSize, const u64& receiverSize, const u64& width, span<block> receiverSet, u64& hashLengthInBytes);// prng用于baseOT

		void itemsToBins(span<Channel> chls, block& commonSeed, span<block> receiverSet, const u64& receiverSize, const u64& width, u64& hashLengthInBytes, std::vector<std::array<block, 2> >& otMessages);

		void run(Channel& ch, block& commonSeed, const u64& senderSize, const u64& receiverSize, const u64& height, const u64& logHeight, span<block> receiverSet, u64 bIdx, u64& numsInBin, const u64& width, u64& hashLengthInBytes, std::vector<std::array<block, 2> >& otMessages);
	};
}

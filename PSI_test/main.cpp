#include "../PSI/include/LDefines.h"
#include "../PSI/include/LUtils.h"
#include "../PSI/include/LSender.h"
#include "../PSI/include/LReceiver.h"
#include "../PSI/tools/BalancedIndex.h"
#include "../PSI/tools/SimpleIndex.h"

#include <cryptoTools/Network/IOService.h>
#include <cryptoTools/Network/Endpoint.h>
#include <cryptoTools/Network/Channel.h>
#include <cryptoTools/Common/CLP.h>
#include <cryptoTools/Common/Log.h>
#include <cryptoTools/Common/Defines.h>
#include <vector>
#include <thread>
//using namespace osuCrypto;

//using namespace osuCrypto;

using namespace std;
using namespace scuPSI;

block commonSeed = oc::toBlock(123456);
u64 senderSize;
u64 receiverSize;
u64 width;//定义时为变量，而使用时作为常量引用
u64 height;
u64 logHeight;
u64 hashLengthInBytes;
u64 bucket, bucket1, bucket2;
u64 numThreads;
string ip;


void runSender() {
	// set up networking
	IOService ios;
	Endpoint ep(ios, ip, EpMode::Server, "test-psi");
	std::vector<Channel> sendChls(numThreads);
	for (u64 i = 0; i < numThreads; ++i)
		sendChls[i] = ep.addChannel("chl" + std::to_string(i), "chl" + std::to_string(i));

	//随机生成发送方的元素
	//std::cout << "Sender:random set generating start" << std::endl;
	std::vector<block> senderSet(senderSize);
	PRNG prng(oc::toBlock(123));
	for (auto i = 0; i < senderSize; ++i) {
		senderSet[i] = prng.get<block>();
	}
	std::cout << "Sender:random set generating finished" << std::endl;

	LSender lSender;
	lSender.batchOutput(prng, sendChls, commonSeed, senderSize, width, senderSet, hashLengthInBytes);
	//lSender.run(prng, ch, commonSeed, senderSize, receiverSize, height, logHeight, width, senderSet, hashLengthInBytes, 32, bucket1, bucket2);

	for (u64 i = 0; i < numThreads; ++i)
		sendChls[i].close();
	ep.stop();
	ios.stop();
}

void runReceiver() {
	// set up networking
	IOService ios;
	Endpoint ep(ios, ip, EpMode::Client, "test-psi");

	std::vector<Channel> recvChls(numThreads);
	for (u64 i = 0; i < numThreads; ++i)
		recvChls[i] = ep.addChannel("chl" + std::to_string(i), "chl" + std::to_string(i));

	//生成100个相同元素
	//std::cout << "Receiver:random set generating start(100 same items)" << std::endl;
	vector<block> receiverSet(receiverSize);
	PRNG prng(oc::toBlock(123));
	for (auto i = 0; i < 100; ++i) {
		receiverSet[i] = prng.get<block>();
	}

	//其余元素不同
	PRNG prng2(oc::toBlock(456));
	for (auto i = 100; i < receiverSize; ++i) {
		receiverSet[i] = prng2.get<block>();
	}
	std::cout << "Receiver:random set generating finished" << std::endl;

	LReceiver lReceiver;
	lReceiver.batchOutput(prng, recvChls, commonSeed, senderSize, receiverSize, width, receiverSet, hashLengthInBytes);
	//lReceiver.run(prng, ch, commonSeed, senderSize, receiverSize, height, logHeight, width, receiverSet, hashLengthInBytes, 32, bucket1, bucket2);

	for (u64 i = 0; i < numThreads; ++i)
		recvChls[i].close();
	ep.stop();
	ios.stop();
}

void Balanced_Hashing_Test()
{
	osuCrypto::setThreadName("Sender");
	u64 setSize = 1 << 20, psiSecParam = 40, numThreads(1);
	PRNG prng(_mm_set_epi32(4253465, 3434565, 234435, 23987045));

	std::vector<block> set(setSize);
	for (u64 i = 0; i < set.size(); ++i)
		set[i] = prng.get<block>();

	BalancedIndex balance;
	osuCrypto::gTimer.reset();
	osuCrypto::gTimer.setTimePoint("start");
	balance.init(setSize, 40, 1);
	balance.insertItems(set);
	osuCrypto::gTimer.setTimePoint("end");
	std::cout << osuCrypto::gTimer << std::endl;
	balance.check();
	balance.toFile(set);
}

void Simple_Hashing_Test() {
	osuCrypto::setThreadName("Sender");
	u64 setSize = 1 << 20, psiSecParam = 40, numThreads(1);
	PRNG prng(_mm_set_epi32(4253465, 3434565, 234435, 23987045));

	std::vector<block> set(setSize);
	for (u64 i = 0; i < set.size(); ++i)
		set[i] = prng.get<block>();

	SimpleIndex simple;
	osuCrypto::gTimer.reset();
	osuCrypto::gTimer.setTimePoint("start");
	simple.init(setSize, 40, 1);
	simple.insertItems(set);
	osuCrypto::gTimer.setTimePoint("end");
	std::cout << osuCrypto::gTimer << std::endl;
	simple.toFile(set);
}


int main(int argc, char** argv) {

#if	0
	Simple_Hashing_Test();
#endif // ENABLE_SIMPLE_HASHING

#if	0	
	Balanced_Hashing_Test();
#endif // ENABLE_BALANCED_HASHING

#if	1
	oc::CLP cmd;
	cmd.parse(argc, argv);

	cmd.setDefault("ss", 20);
	senderSize = 1 << cmd.get<u64>("ss");//(2^20)

	cmd.setDefault("rs", 20);
	receiverSize = 1 << cmd.get<u64>("rs");//(2^20)

	cmd.setDefault("w", 632);
	width = cmd.get<u64>("w");

	cmd.setDefault("h", 20);
	logHeight = cmd.get<u64>("h");
	height = 1 << cmd.get<u64>("h");//(2^20)

	cmd.setDefault("hash", 10);
	hashLengthInBytes = cmd.get<u64>("hash");

	cmd.setDefault("ip", "localhost");
	ip = cmd.get<string>("ip");

	cmd.setDefault("thds", 4);
	numThreads = cmd.get<u64>("thds");

	bucket1 = bucket2 = 1 << 8;//(2^8=256)

	bool noneSet = !cmd.isSet("r");
	if (noneSet) {
		std::cout
			<< "=================================\n"
			<< "||  Private Set Intersection   ||\n"
			<< "=================================\n"
			<< "\n"
			<< "This program reports the performance of the private set intersection protocol.\n"
			<< "\n"
			<< "Experimenet flag:\n"
			<< " -r 0    to run a sender.\n"
			<< " -r 1    to run a receiver.\n"
			<< "\n"
			<< "Parameters:\n"
			<< " -ss     log(#elements) on sender side.\n"
			<< " -rs     log(#elements) on receiver side.\n"
			<< " -w      width of the matrix.\n"
			<< " -h      log(height) of the matrix.\n"
			<< " -hash   hash output length in bytes.\n"
			<< " -ip     ip address (and port).\n"
			;
	}
	else {
		if (cmd.get<u64>("r") == 0) {
			std::cout << "runSender()" << std::endl;
			runSender();
		}
		else if (cmd.get<u64>("r") == 1) {
			std::cout << "runReceiver()" << std::endl;
			runReceiver();
		}
	}
	return 0;
#endif // ENABLE_LPSI
}

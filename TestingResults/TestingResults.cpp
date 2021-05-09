
#include <fstream>
#include "AlgoritmicHeights.hpp"
#include "BasicGraph.hpp"
#include "HelperFunctions.hpp"

#include "CppUnitTest.h"


using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace TestingResults
{
	TEST_CLASS(TestingResults)
	{
		template<typename T>
		bool areSame(std::vector<T> v1, std::vector<T> v2) {
			if (v1.size() != v2.size()) return false;
			for (int i = 0; i < v1.size(); ++i) {
				if (v1[i] != v2[i]) return false;
			}
			return true;
		}

	public:
		
		TEST_METHOD(TestMethod1)
		{
			std::string input, result;
			input = "../Rosalind/AlgorithmicHeights/rosalind_bfs.txt";
			result = "../Rosalind/AlgorithmicHeights/rosalind_bfs_r.txt";
			std::ifstream inputStream(input, std::ios::binary);
			std::ifstream resultStream(result, std::ios::binary);
			auto graph = directedGraph<int, int>(inputStream);
			graph.BFS(1);
			auto test = graph.getDistances();
			auto correct = list<int>(resultStream);
			Assert::IsTrue(areSame(correct, test));



		}
	};
}

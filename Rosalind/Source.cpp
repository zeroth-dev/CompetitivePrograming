#include <fstream>
#include <algorithm>
#include <istream>
#include <chrono>
#include "BasicGraph.hpp"
#include "BioStronghold.hpp"
#include "AlgoritmicHeights.hpp"
#include "HelperFunctions.hpp"
//#include "Networking.hpp"


int main(int argc, char** argv) {
    std::string filename;
    if (argc == 1)
    #if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
        filename = "AlgorithmicHeights/rosalind_sdag.txt";
    #else
        filename = "../AlgorithmicHeights/rosalind_bfs.txt";
    #endif
//        filename = "AlgorithmicHeights/rosalind_dij.txt";
    else if (argc == 2)
        filename = argv[1];
    else{
            cerr << "Provide only input file as argument." << endl;
            return -1;
        }
    // open file for reading
    std::ifstream istrm(filename, std::ios::binary);
    if (!istrm.is_open()) {
        std::cout << "failed to open " << filename << '\n';
    }
    else {
        int nodeNumber;
        //auto adjList = directedAdjListUnity(istrm, nodeNumber);
        //auto output = ts(adjList, nodeNumber);
        
        auto adjList = directedAdjList(istrm, nodeNumber);
        auto output = sdag(adjList, nodeNumber);
        for (auto it : output) {
            if (it == INT_MAX) {
                std::cout << "x ";
            }
            else {
                std::cout << it << " ";
            }
        }
    }
    istrm.close();
	return 0;
}

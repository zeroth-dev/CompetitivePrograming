#include <fstream>
#include <algorithm>
#include <istream>
#include <chrono>
#include "BasicGraph.hpp"
#include "BioStronghold.hpp"
#include "AlgoritmicHeights.hpp"
//#include "Networking.hpp"


struct tokens : std::ctype<char>
{
    tokens(std::vector<char> chr) : std::ctype<char>(get_table(chr)) {}

    static std::ctype_base::mask const* get_table(std::vector<char> chr)
    {
        typedef std::ctype<char> cctype;
        static const cctype::mask* const_rc = cctype::classic_table();

        static cctype::mask rc[cctype::table_size];
        std::memcpy(rc, const_rc, cctype::table_size * sizeof(cctype::mask));
        for (auto c : chr) {
            rc[c] = std::ctype_base::space;
        }
        
        return &rc[0];
    }
};

pair < vector<string>, vector<string>> processFasta(istream &istr);
std::vector<int> processSet(istream& istrm);

int main(int argc, char** argv) {
    std::string filename;
    if (argc == 1)
        filename = "AlgorithmicHeights/rosalind_bfs.txt";
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
        // input for adjacency list
        
        
        int k, n, m, ver1, ver2, dist;
        k = 1;
        std::multimap<int, std::pair<int, int>> adjList;
        BasicGraph<int> graph;
        for (int i = 0; i < k; ++i) {
            istrm >> n >> m;
            for (int i = 0; i < m; ++i) {
                istrm >> ver1 >> ver2;
                graph.addUndirectedEdge(ver1, ver2);
            }
            graph.BFS(1);
            graph.printDistances();
        }
        
        /*
        int n, temp;
        istrm >> n;
        std::vector<int> input;
        for (int i = 0; i < n; ++i) {
            istrm >> temp;
            input.push_back(temp);
        }
        auto output = hs(input); 
        for (auto it : output) {
            std::cout << it << " ";
        
        }
        */
        /*
        int n, m, temp;
        istrm >> n >> m;
        std::vector<int> input;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                istrm >> temp;
                input.push_back(temp);
            }
            auto output = sum3(input);
            for (auto it : output) {
                std::cout << it << " ";
            }
            std::cout << "\n";
            input.clear();
        }
        */
       
    }
    
    istrm.close();
	return 0;
}


// TODO:
// Make another function just for strings

pair < vector<string>, vector<string>> processFasta(istream &istr) {
    vector<string> ids;
    vector<string> allDNA;
    string id;
    getline(istr, id);
    // Erase '>' character at the beginning
    id.erase(0, 1);
    ids.push_back(id);
    string line;
    vector<string> dnaStrings;
    while (getline(istr, line)) {
        if (line.find(">") != string::npos) {
            // Save the DNA strand
            string DNAstrand = "";
            for (auto it : dnaStrings) {
                DNAstrand.append(it);
            }
            allDNA.push_back(DNAstrand);
            dnaStrings.clear();
            // Erase '>' character at the beginning
            line.erase(0, 1);
            ids.push_back(line);
            getline(istr, line);
        }
        if (line[line.size() - 1] == '\r') line.pop_back();
        dnaStrings.push_back(line);
    }
    if (dnaStrings.size()) {
        string DNAstrand = "";
        for (auto it : dnaStrings) {
            DNAstrand.append(it);
        }
        allDNA.push_back(DNAstrand);
        dnaStrings.clear();
    }
    return make_pair(ids, allDNA);
}

// Process set of numbers in form {1, 2, 3, ...}

std::vector<int> processSet(istream &istrm) {
    std::vector<char> chars = { '{', ',' };
    istrm.imbue(std::locale(std::locale(), new tokens(chars)));
    std::vector<int> output;
    std::string temp;
    istrm >> temp;
    while (temp.back() != '}') {
        output.push_back(std::stoi(temp));
        istrm >> temp;
    }
    if (temp.size() > 1) {
        temp.pop_back();

        output.push_back(std::stoi(temp));
    }
    return output;
}
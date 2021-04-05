#include <fstream>
#include <algorithm>
#include <istream>
#include "Functions.hpp"
#include "Networking.hpp"

pair < vector<string>, vector<string>> processFasta(istream &istr);
int main(int argc, char** argv) {
    std::string filename;
    if (argc == 1)
        filename = "../datasets/rosalind_mrna.txt";
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
        std::string codon;
        istrm >> codon;
        std::cout << RNACombinations(codon) << std::endl;
        
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
        if (line.find("Rosalind") != string::npos) {
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
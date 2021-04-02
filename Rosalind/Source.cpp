#include <fstream>
#include <algorithm>
#include "Functions.hpp"

pair < vector<string>, vector<string>> processFasta(istream &istr);
int main() {


    std::string filename = "D:\\Competitive programming\\CompetitivePrograming\\Rosalind\\datasets\\rosalind_mprt.txt";
   
    // open file for reading
    std::ifstream istrm(filename, std::ios::binary);
    if (!istrm.is_open()) {
        std::cout << "failed to open " << filename << '\n';
    }
    else {
        vector<string> ids, dnaseq;
        auto processed = processFasta(istrm);
        ids = processed.first;
        dnaseq = processed.second;
        dnaseq[0][83] = ' ';
        cout << dnaseq[0];
    }
    istrm.close();
	return 0;
}





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
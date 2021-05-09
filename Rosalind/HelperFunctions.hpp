// Place for all helper functions for loading data
#include <fstream>
#include <algorithm>
#include <istream>
#include <vector>
#include <string>
#include "BasicGraph.hpp"

#ifndef HELPERFUNCTIONS_H
#define HELPERFUNCTIONS_H


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



// TODO:
// Make another function just for strings

std::pair < std::vector<std::string>, std::vector<std::string>> processFasta(std::istream& istr) {
    std::vector<std::string> ids;
    std::vector<std::string> allDNA;
    std::string id;
    std::getline(istr, id);
    // Erase '>' character at the beginning
    id.erase(0, 1);
    ids.push_back(id);
    std::string line;
    std::vector<std::string> dnaStrings;
    while (std::getline(istr, line)) {
        if (line.find(">") != std::string::npos) {
            // Save the DNA strand
            std::string DNAstrand = "";
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
        std::string DNAstrand = "";
        for (auto it : dnaStrings) {
            DNAstrand.append(it);
        }
        allDNA.push_back(DNAstrand);
        dnaStrings.clear();
    }
    return make_pair(ids, allDNA);
}

// Process set of numbers in form {1, 2, 3, ...}

std::vector<int> processSet(std::istream& istrm) {
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


std::multimap<int, std::pair<int, int>> directedAdjList(std::istream& istrm, int& nodeNumber) {

    int m, ver1, ver2, dist;
    
    std::multimap<int, std::pair<int, int>> list;
    
    istrm >> nodeNumber >> m;
    for (int i = 0; i < m; ++i) {
        istrm >> ver1 >> ver2 >> dist;
        list.insert({ ver1, { ver2, dist } });
    }
    return list;
}

std::multimap<int, int> directedAdjListUnity(std::istream& istrm, int& nodeNumber) {

    int m, ver1, ver2;

    std::multimap<int, int> list;

    istrm >> nodeNumber >> m;
    for (int i = 0; i < m; ++i) {
        istrm >> ver1 >> ver2;
        list.insert({ ver1, ver2 });
    }
    return list;
}

std::multimap<int, int> undirectedAdjList(std::istream& istrm, int &nodeNumber) {

    int m, ver1, ver2;

    std::multimap<int, int> list;

    istrm >> nodeNumber >> m;
    for (int i = 0; i < m; ++i) {
        istrm >> ver1 >> ver2;
        list.insert({ ver1, ver2 });
        list.insert({ ver2, ver1 });
    }
    return list;
}

template <typename T, typename I>
BasicGraph<T, I> directedGraph(std::istream& istrm, bool unity = true) {

    int n, m;
    T ver1, ver2;
    I dist = 1;
    BasicGraph<T, I> graph;
    istrm >> n >> m;
    for (int i = 0; i < m; ++i) {
        istrm >> ver1 >> ver2;
        if (!unity) istrm >> dist;
        graph.addDirectedEdge(ver1, ver2, dist);
    }
    return graph;
    
}

template <typename T, typename I>
BasicGraph<T, I> undirectedGraph(std::istream& istrm) {

    int n, m;
    T ver1, ver2;
    BasicGraph<T, I> graph;
    istrm >> n >> m;
    for (int i = 0; i < m; ++i) {
        istrm >> ver1 >> ver2;
        graph.addUnirectedEdge(ver1, ver2);
    }
    return graph;

}

template<typename T>
std::vector<T> list(std::istream& istrm) {
    T temp;
    std::vector<T> output;
    while(istrm >> temp){
        output.push_back(temp);
    }
    return output;
}
template<typename T>
std::vector<T> nodeDistFormat(std::istream& istrm) {
    T temp;
    std::vector<T> output;
    while (istrm >> temp) {
        istrm >> temp;
        output.push_back(temp);
    }
    return output;
}

#endif // !HELPERFUNCTIONS_H

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <queue>
#include <algorithm>
using namespace std;

uint64_t NcK(int n, int k) {
    uint64_t answer = 1;
    if (k > n / 2) k = n - k;

    for (long long i = 1; i <= k; i++) {
        answer *= n - (uint64_t)(i - 1);
        answer /= i;
    }
    return answer;
}

double binomialProbability(int n, int k, double p) {
    return (double)NcK(n, k) * pow(p, k) * pow(1.0 - p, n - k);
}

int countCG(string str) {
    int count = 0;
    for (auto it : str) {
        if (it == 'G' || it == 'C') count++;
    }
    return count;
}

vector<int> countChars(string str) {
    vector<int> number_of_chars;
    map<char, int> appearance;
    for (auto ch : str) {
        if (ch == '\r') continue;
        if (appearance.find(ch) == appearance.end()) {
            appearance.insert({ ch, 1 });
        }
        else {
            appearance[ch]++;
        }
    }
    for (auto it : appearance) {
        number_of_chars.push_back(it.second);
    }
    return number_of_chars;
}

string DNA2RNA(string dna) {
    string rna;
    for (char c : dna) {
        if (c == '\r') continue;
        if (c == 'T') {
            rna.push_back('U');
        }
        else rna.push_back(c);
    }
    return rna;
}

string reverseComplementDNA(string dna) {
    string reverseC;
    for (auto ch = dna.rbegin(); ch != dna.rend(); ch++) {
        if (*ch == '\r') continue;
        switch (*ch) {
        case 'A':
            reverseC.push_back('T');
            break;
        case 'T':
            reverseC.push_back('A');
            break;
        case 'C':
            reverseC.push_back('G');
            break;
        case 'G':
            reverseC.push_back('C');
            break;
        }
    }
    return reverseC;
}

// Each adult rabbit gives birth to k new rabbits
uint64_t getFib_n_k(int n, int k) {
    int curr = 1;

    vector<uint64_t> Fibonacci = { 0, 1 };
    while (Fibonacci.size() < n + 1) {
        uint64_t next = Fibonacci[curr] + k * Fibonacci[curr - 1];
        curr++;
        Fibonacci.push_back(next);
    }
    return Fibonacci[n];
}

// Each rabbit dies after k months, assume k>2
// NOT SOLVED YET
/*
uint64_t getFib_n_dk(int n, int k) {

    int curr = 1;
    int months_passed = 0;
    vector<uint64_t> bornRabbits = {1};
    vector<uint64_t> Fibonacci = { 0, 1 };
    while (Fibonacci.size() < n + 1) {
        uint64_t next = Fibonacci[curr] + Fibonacci[curr - 1];
        curr++;
        bornRabbits.push_back(Fibonacci[curr-1])
        if (months_passed < k) months_passed++;
        Fibonacci.push_back(next);
    }
    return Fibonacci[n];
}
*/

double GC_content(string dna) {
    string cpy = dna;
    cpy.erase(remove(cpy.begin(), cpy.end(), '\r'), cpy.end());
    double size = cpy.size();
    double gc = 0;
    for (char c : cpy) {
        if (c == 'G' || c == 'C') gc++;
    }
    return 100.0 * gc / size;
}

int hammingDistance(string first, string second) {
    int count = 0;
    for (int i = 0; i < first.size(); i++) {
        if (first[i] != second[i]) count++;
    }
    return count;
}
// k is dominant homogeneus, m is heterogeneous, n is recesive homogeneus
double dominantProb(int n, int m, int k) {
    double sum = n + m + k;
    int nn = n * (n - 1) / 2;
    int mm = m * (m - 1) / 2;
    int kk = k * (k - 1) / 2;
    int nm = n * m;
    int nk = n * k;
    int mk = m * k;

    double prob = (2.0 / (sum*(sum-1))) * (1 * kk + 3.0 / 4.0 * mm + 0 * nn + 1 * nk + 0.5 * nm + 1 * mk);
    return prob;

}

map<string, string> RNACodon;
void initRNACodon(); 
string RNA2prot(string RNA) {
    string prot = "";
    initRNACodon();
    for (int i = 0; i < RNA.size(); i+=3) {
        string rna = "aaa";
        rna[0] = RNA[i];
        rna[1] = RNA[i + 1];
        rna[2] = RNA[i + 2];
        if (RNACodon[rna] == "Stop") break;
        prot.append(RNACodon[rna]);
    }
    return prot;
}

vector<int> findSubstring(string str, string substr) {
    vector<int> indexes;
    int position = 0;
    while (true)
    {
        position = str.find(substr, position);
        if (position == string::npos) break;
        indexes.push_back(++position);
    }
    return indexes;
}

pair<vector<vector<int>>, string> DNAConcensus(vector<string> dnaseq) {
    int dnaSize = dnaseq[0].size();
    string concensus = "";
    vector<vector<int>> allProteins;
    for (int i = 0; i < dnaSize; i++) {
        int A = 0, C = 0, G = 0, T = 0, maximum;
        vector<int> numberOfProteins; // Order is ACGT
        for (int j = 0; j < dnaseq.size(); j++) {
            switch (dnaseq[j][i]) {
            case 'A':
                A++;
                break;
            case 'C':
                C++;
                break;
            case 'G':
                G++;
                break;
            case 'T':
                T++;
                break;
            default:
                break;
            }
        }
        numberOfProteins.push_back(A);
        numberOfProteins.push_back(C);
        numberOfProteins.push_back(G);
        numberOfProteins.push_back(T);
        allProteins.push_back(numberOfProteins);
        maximum = max({ A, C, G, T });
        if (maximum == A) {
            concensus.append("A");
        }
        else if (maximum == C) {
            concensus.append("C");
        }
        else if (maximum == G) {
            concensus.append("G");
        }
        else if (maximum == T) {
            concensus.append("T");
        }

    }

    return make_pair(allProteins, concensus);
}

vector<pair<int, int>> adjacencyO3(vector<string> dnaseq) {
    vector<pair<int, int>> adjacencyList;
    vector<string> start;
    vector<string> end;
    for (auto dna : dnaseq) {
        start.push_back(dna.substr(0, 3));
        end.push_back(dna.substr(dna.size() - 3, 3));
    }
    for (int i = 0; i < start.size(); i++) {
        for (int j = 0; j < end.size(); j++) {
            if (i == j) continue;
            if (start[i] == end[j]) adjacencyList.push_back(make_pair(j, i));
        }
    }
    return adjacencyList;
}

// Regex exercise 

void initRNACodon() {

    RNACodon.insert({ "UUU", "F" });    RNACodon.insert({ "CUU", "L" }); RNACodon.insert({ "AUU", "I" }); RNACodon.insert({ "GUU", "V" });
    RNACodon.insert({ "UUC", "F" });    RNACodon.insert({ "CUC", "L" }); RNACodon.insert({ "AUC", "I" }); RNACodon.insert({ "GUC", "V" });
    RNACodon.insert({ "UUA", "L" });    RNACodon.insert({ "CUA", "L" }); RNACodon.insert({ "AUA", "I" }); RNACodon.insert({ "GUA", "V" });
    RNACodon.insert({ "UUG", "L" });    RNACodon.insert({ "CUG", "L" }); RNACodon.insert({ "AUG", "M" }); RNACodon.insert({ "GUG", "V" });
    RNACodon.insert({ "UCU", "S" });    RNACodon.insert({ "CCU", "P" }); RNACodon.insert({ "ACU", "T" }); RNACodon.insert({ "GCU", "A" });
    RNACodon.insert({ "UCC", "S" });    RNACodon.insert({ "CCC", "P" }); RNACodon.insert({ "ACC", "T" }); RNACodon.insert({ "GCC", "A" });
    RNACodon.insert({ "UCA", "S" });    RNACodon.insert({ "CCA", "P" }); RNACodon.insert({ "ACA", "T" }); RNACodon.insert({ "GCA", "A" });
    RNACodon.insert({ "UCG", "S" });    RNACodon.insert({ "CCG", "P" }); RNACodon.insert({ "ACG", "T" }); RNACodon.insert({ "GCG", "A" });
    RNACodon.insert({ "UAU", "Y" });    RNACodon.insert({ "CAU", "H" }); RNACodon.insert({ "AAU", "N" }); RNACodon.insert({ "GAU", "D" });
    RNACodon.insert({ "UAC", "Y" });    RNACodon.insert({ "CAC", "H" }); RNACodon.insert({ "AAC", "N" }); RNACodon.insert({ "GAC", "D" });
    RNACodon.insert({ "UAA", "Stop" }); RNACodon.insert({ "CAA", "Q" }); RNACodon.insert({ "AAA", "K" }); RNACodon.insert({ "GAA", "E" });
    RNACodon.insert({ "UAG", "Stop" }); RNACodon.insert({ "CAG", "Q" }); RNACodon.insert({ "AAG", "K" }); RNACodon.insert({ "GAG", "E" });
    RNACodon.insert({ "UGU", "C" });    RNACodon.insert({ "CGU", "R" }); RNACodon.insert({ "AGU", "S" }); RNACodon.insert({ "GGU", "G" });
    RNACodon.insert({ "UGC", "C" });    RNACodon.insert({ "CGC", "R" }); RNACodon.insert({ "AGC", "S" }); RNACodon.insert({ "GGC", "G" });
    RNACodon.insert({ "UGA", "Stop" }); RNACodon.insert({ "CGA", "R" }); RNACodon.insert({ "AGA", "R" }); RNACodon.insert({ "GGA", "G" });
    RNACodon.insert({ "UGG", "W" });    RNACodon.insert({ "CGG", "R" }); RNACodon.insert({ "AGG", "R" }); RNACodon.insert({ "GGG", "G" });
    
      
}
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <queue>
#include <algorithm>
#include <math.h>
#include <regex>
using namespace std;

uint64_t NcK(int n, int k, bool modMilion = false) {
    uint64_t answer = 1;
    if (k > n / 2) k = n - k;

    for (long long i = 1; i <= k; i++) {
        answer *= n - (uint64_t)(i - 1);
        answer /= i;
        answer = modMilion ? answer % 1000000 : answer;
    }
    return answer;
}


uint64_t fact(int n, bool modMilion = false){
    if(n <= 1){
        return 1;
    }else{
        if (modMilion){
            return n * (fact(n-1) % 1000000);
        }else{
            return n * fact(n-1);
        }
    }
}

/// Given a map from keys to values, creates a new map from values to keys 
template<typename K, typename V>
static map<V, K> reverse_map(const map<K, V>& m) {
    map<V, K> r;
    for (const auto& kv : m)
        r[kv.second] = kv.first;
    return r;
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
std::map<std::string, std::string> DNACodon;
void initDNACodon(); 
std::map<std::string, double> AminoMassTable;
void initAminoMassTable();


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
std::string DNA2prot(std::string DNA) {
    std::string prot = "";
    initDNACodon();
    for (int i = 0; i < DNA.size(); i+=3) {
        std::string dna = "aaa";
        dna[0] = DNA[i];
        dna[1] = DNA[i + 1];
        dna[2] = DNA[i + 2];
        if (DNACodon[dna] == "Stop") break;
        prot.append(DNACodon[dna]);
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

// TODO:
// Construct regex based on bioinf text

std::vector<int> findMotifs(std::string dna, std::string motif){
    std::regex rgx(motif);
    std::vector<int> positions;
    for(std::sregex_iterator it (dna.begin(), dna.end(), rgx), it_end; it != it_end; it++){
        positions.push_back(it->position()+1);
    }
    return positions;
}

int RNACombinations(std::string protein){
    std::map <std::string, int> numberOfCombinations;
    initRNACodon();

    for (auto it : RNACodon){
        if (numberOfCombinations.find(it.second) == numberOfCombinations.end()){
            numberOfCombinations.insert({it.second, 1});
        }else{
            numberOfCombinations[it.second]++;
        }
    } 

    int combo = 1;
    int milion = 1000000;
    for(auto codon : protein){
        
        combo*=numberOfCombinations[std::string(1, codon)];
        combo%= milion;
    }
    combo*=numberOfCombinations["Stop"];
    combo%=milion;
    return combo;
}


std::vector<std::vector<int>> allPermutations(int n){
    uint64_t permutations = fact(n);
    std::cout << permutations <<std::endl;

    std::vector<int> currentPermutation;
    std::vector<std::vector<int>> allPermutations;
    for(int i = 1; i <=n; i++){
        currentPermutation.push_back(i);
    }
    do{
        allPermutations.push_back(currentPermutation);
    }while(std::next_permutation(currentPermutation.begin(), currentPermutation.end()));

    return allPermutations;
}

double proteinMass(std::string protein){
    initAminoMassTable();
    double totalMass = 0;
    for (auto aminoAcid : protein){
        totalMass+=AminoMassTable[std::string(1, aminoAcid)];
    }
    // Total mass is te sum of masses of amino acids plus a single molecule of water
    //totalMass+=AminoMassTable["H2O"];
    return totalMass;
}

std::vector<std::pair<int, int>> findReversePalin(std::string DNA){
    std::vector<std::pair<int, int>> result;
    for(int i = 0; i < DNA.size(); i++){
        for(int j = 4; j <=12 && i+j-1 < DNA.size(); j++){
            std::string temp = "";
            for(int k = 0; k < j; k++){
                temp.push_back(DNA[i+k]);
            }
            if(temp == reverseComplementDNA(temp)) result.push_back({i+1, j}); 
        }
    }
    return result;
}

std::vector<std::string> proteinCandidates(std::string dna){
    std::vector<std::string> candidates;
    std::string reverse = reverseComplementDNA(dna);
    initDNACodon();
    for(int i = 0; i+2 < dna.size(); i++){
        std::string triplet = "aaa";
        triplet[0] = dna[i];
        triplet[1] = dna[i+1];
        triplet[2] = dna[i+2];
        if (triplet == "ATG"){
            std::string protein = "";
            protein.append(DNACodon["ATG"]);
            for (int j = i+3; j < dna.size(); j+=3){
                std::string triplet = "aaa";
                triplet[0] = dna[j];
                triplet[1] = dna[j+1];
                triplet[2] = dna[j+2];
                if(DNACodon[triplet] == "Stop"){
                    if(std::find(candidates.begin(), candidates.end(), protein) == candidates.end())
                        candidates.push_back(protein);
                    break;
                }
                protein.append(DNACodon[triplet]);
            }
        }
    }
    for(int i = 0; i+2 < dna.size(); i++){
        std::string triplet = "aaa";
        triplet[0] = reverse[i];
        triplet[1] = reverse[i+1];
        triplet[2] = reverse[i+2];
        if (triplet == "ATG"){
            std::string protein = "";
            protein.append(DNACodon["ATG"]);
            for (int j = i+3; j < dna.size(); j+=3){
                std::string triplet = "aaa";
                triplet[0] = reverse[j];
                triplet[1] = reverse[j+1];
                triplet[2] = reverse[j+2];
                if(DNACodon[triplet] == "Stop"){
                    if(std::find(candidates.begin(), candidates.end(), protein) == candidates.end())
                        candidates.push_back(protein);
                    break;
                }
                protein.append(DNACodon[triplet]);
            }
        }
    }



    return candidates;
}

std::string splicedProtein(std::vector<std::string> totalDNA){
    std::string encodingDNA = totalDNA[0];
    for(int i = 1; i < totalDNA.size(); i++){
        int location = encodingDNA.find(totalDNA[i]);
        if(location != std::string::npos){
            int size = totalDNA[i].size();
            auto firstPart = 
            encodingDNA = encodingDNA.substr(0, location).append(encodingDNA.substr(location+size, totalDNA.size()-size-location));
        }
    }
    
    return DNA2prot(encodingDNA);
}

std::vector<int> longestIncSubseq(std::vector<int>){
    return std::vector<int>(5);
}

std::vector<int> longestDecSubseq(std::vector<int>){
    return std::vector<int>(5);
}

// Not done, something is wrong
std::string superString(std::vector<std::string> allDNA){
    std::string dnaSequence = allDNA[0];
    for(int i = 1; i < allDNA.size(); i++){
        std::string temp = allDNA[i];
        for(int j = temp.size(); j>temp.size()/2; j--){
            if(j > dnaSequence.size() ) continue;
            if(dnaSequence.find(temp.substr(0, temp.size()))!= std::string::npos) break;
            if(temp.substr(0, j) == dnaSequence.substr(dnaSequence.size()-j, j)){
                dnaSequence.append(temp.substr(j, temp.size()-j));
                break;
            }if(j == temp.size()/2+1){
                dnaSequence.append(temp);
            }
            
        }
        //cout << endl <<"dna sequence: " << endl << dnaSequence << endl << endl;
    }
    return dnaSequence;
}

// TODO 
// Implement BigNumbers

double numberOfPerfectMatchings(std::string rna){
    double A = 0, G = 0;
    for (auto it : rna){
        if (it == 'A') A++;
        else if (it == 'G') G++;
    }
    cout << "A: " << A << " G: " << G << endl;
    cout << (double)fact(A) << endl << (double)fact(G) << endl;
    return (double)fact(A)*(double)fact(G);
}

int numberOfPartialPermutations(int n, int k){
    cout << (fact(n, true) % 1000000) << endl;
    cout << (NcK(n, k) % 1000000) << endl;
    return (fact(n, true) % 1000000) * (NcK(n, k) % 1000000) % 1000000;
}



void initAminoMassTable(){

    AminoMassTable.insert({"A", 71.03711});     AminoMassTable.insert({"C", 103.00919});
    AminoMassTable.insert({"D", 115.02694});    AminoMassTable.insert({"E", 129.04259});
    AminoMassTable.insert({"F", 147.06841});    AminoMassTable.insert({"G", 57.02146});
    AminoMassTable.insert({"H", 137.05891});    AminoMassTable.insert({"I", 113.08406});
    AminoMassTable.insert({"K", 128.09496});    AminoMassTable.insert({"L", 113.08406});
    AminoMassTable.insert({"M", 131.04049});    AminoMassTable.insert({"N", 114.04293});
    AminoMassTable.insert({"P", 97.05276});     AminoMassTable.insert({"Q", 128.05858});
    AminoMassTable.insert({"R", 156.10111});    AminoMassTable.insert({"S", 87.03203});
    AminoMassTable.insert({"T", 101.04768});    AminoMassTable.insert({"V", 99.06841});
    AminoMassTable.insert({"W", 186.07931});    AminoMassTable.insert({"Y", 163.06333});
    AminoMassTable.insert({"H2O", 18.01056});                  

}

void initDNACodon(){
    DNACodon.insert({ "TTT", "F" });    DNACodon.insert({ "CTT", "L" }); DNACodon.insert({ "ATT", "I" }); DNACodon.insert({ "GTT", "V" });
    DNACodon.insert({ "TTC", "F" });    DNACodon.insert({ "CTC", "L" }); DNACodon.insert({ "ATC", "I" }); DNACodon.insert({ "GTC", "V" });
    DNACodon.insert({ "TTA", "L" });    DNACodon.insert({ "CTA", "L" }); DNACodon.insert({ "ATA", "I" }); DNACodon.insert({ "GTA", "V" });
    DNACodon.insert({ "TTG", "L" });    DNACodon.insert({ "CTG", "L" }); DNACodon.insert({ "ATG", "M" }); DNACodon.insert({ "GTG", "V" });
    DNACodon.insert({ "TCT", "S" });    DNACodon.insert({ "CCT", "P" }); DNACodon.insert({ "ACT", "T" }); DNACodon.insert({ "GCT", "A" });
    DNACodon.insert({ "TCC", "S" });    DNACodon.insert({ "CCC", "P" }); DNACodon.insert({ "ACC", "T" }); DNACodon.insert({ "GCC", "A" });
    DNACodon.insert({ "TCA", "S" });    DNACodon.insert({ "CCA", "P" }); DNACodon.insert({ "ACA", "T" }); DNACodon.insert({ "GCA", "A" });
    DNACodon.insert({ "TCG", "S" });    DNACodon.insert({ "CCG", "P" }); DNACodon.insert({ "ACG", "T" }); DNACodon.insert({ "GCG", "A" });
    DNACodon.insert({ "TAT", "Y" });    DNACodon.insert({ "CAT", "H" }); DNACodon.insert({ "AAT", "N" }); DNACodon.insert({ "GAT", "D" });
    DNACodon.insert({ "TAC", "Y" });    DNACodon.insert({ "CAC", "H" }); DNACodon.insert({ "AAC", "N" }); DNACodon.insert({ "GAC", "D" });
    DNACodon.insert({ "TAA", "Stop" }); DNACodon.insert({ "CAA", "Q" }); DNACodon.insert({ "AAA", "K" }); DNACodon.insert({ "GAA", "E" });
    DNACodon.insert({ "TAG", "Stop" }); DNACodon.insert({ "CAG", "Q" }); DNACodon.insert({ "AAG", "K" }); DNACodon.insert({ "GAG", "E" });
    DNACodon.insert({ "TGT", "C" });    DNACodon.insert({ "CGT", "R" }); DNACodon.insert({ "AGT", "S" }); DNACodon.insert({ "GGT", "G" });
    DNACodon.insert({ "TGC", "C" });    DNACodon.insert({ "CGC", "R" }); DNACodon.insert({ "AGC", "S" }); DNACodon.insert({ "GGC", "G" });
    DNACodon.insert({ "TGA", "Stop" }); DNACodon.insert({ "CGA", "R" }); DNACodon.insert({ "AGA", "R" }); DNACodon.insert({ "GGA", "G" });
    DNACodon.insert({ "TGG", "W" });    DNACodon.insert({ "CGG", "R" }); DNACodon.insert({ "AGG", "R" }); DNACodon.insert({ "GGG", "G" });
    
}

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
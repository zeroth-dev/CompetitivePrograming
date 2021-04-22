#include <vector>
#include <map>
#include <set>
#include <queue>
#include <stack>
// Here I am writing the function names according to the problem code

int fibo(int n) {
	std::vector<int> fibs;
	fibs.push_back(0);
	fibs.push_back(1);
	for (int i = 2; i <= n; ++i) {
		fibs.push_back(fibs[i - 1] + fibs[i - 2]);
	}
	return fibs[n];
}

std::vector<int> bins(std::vector<int> input, std::vector<int> search) {
	std::vector<int> positions;
	int initialPosition = input.size() / 2;
	int currPosition = 0;
	for (auto num : search) {
		currPosition = input.size() / 2;
		bool notFound = false;
		int leftLimit = 0;
		int rightLimit = input.size();
		while (input[currPosition] != num) {
			if (input[currPosition] < num && input[currPosition + 1] > num) {
				notFound = true;
				break;
			}
			if (input[currPosition] < num) {
				leftLimit = currPosition;
				currPosition += (rightLimit - currPosition) / 2;
			}
			else {
				rightLimit = currPosition;
				currPosition -= (currPosition - leftLimit) / 2;
				if (currPosition - leftLimit == 1) currPosition = leftLimit;
			}
		}
		if (notFound) positions.push_back(-1);
		else positions.push_back(currPosition + 1);
	}
	return positions;
}

// O(n^2) algorithm 
// Can be improved by using a map for adjacency list
std::vector<int> deg(std::vector<std::vector<int>> adjList) {
	int n = adjList.size();
	std::vector<int> degree;
	for (int i = 1; i < n; ++i) {
		int countEdges = 0;
		for (int j = 1; j < n; ++j) {
			if (adjList[i][j] == 1) ++countEdges;
		}
		if (countEdges) degree.push_back(countEdges);
	}
	return degree;
}
// Count the number of swaps
int ins(vector<int> unsortedArray) {
	int countSwaps=0;
	for (int i = 1; i < unsortedArray.size(); ++i) {
		int k = i;
		while (k > 0 && unsortedArray[k] < unsortedArray[k - 1]) {
			std::swap(unsortedArray[k], unsortedArray[k - 1]);
			++countSwaps;
			--k;
		}
	}
	return countSwaps;
}

// First we find the degree of each vertex and then we sum degrees of neighboring nodes
std::vector<int> ddeg(std::multimap<int, int> adjList) {
	int n = adjList.size();
	int nodes = adjList.rbegin()->first;
	std::vector<int> degrees(nodes + 1);
	std::vector<int> doubleDegrees(nodes + 1);
	for (int i = 1; i <= nodes; ++i) {
		auto range = adjList.equal_range(i);
		for (auto it = range.first; it != range.second; ++it) {
			if (it->second == 0) continue;
			degrees[it->first]++;
		}
	}
	for (int i = 1; i <= nodes; ++i) {
		auto range = adjList.equal_range(i);
		for (auto it = range.first; it != range.second; ++it) {
			doubleDegrees[i]+=degrees[it->second];
		}
	}
	return doubleDegrees;
}

std::vector<int> maj(std::vector<std::vector<int>> lists) {
	vector<int> solution;
	for (int i = 0; i < lists.size(); ++i) {
		std::sort(lists[i].begin(), lists[i].end());
		int currNum = lists[i][0];
		int count = 1;
		for (int j = 1; j < lists[i].size(); ++j) {
			if (lists[i][j] == currNum) ++count;
			else {
				currNum = lists[i][j];
				count = 1;
			}
			if (count > lists[i].size() / 2) {
				break;
			}
		}
		if (count > lists[i].size() / 2) {
			solution.push_back(currNum);
		}
		else solution.push_back(-1);
	}
	return solution;
}	 

std::vector<int> mer(std::vector<int> first, std::vector<int> second) {
	std::vector<int> sol;

	for (int i = 0, j = 0; i < first.size() || j < second.size();) {
		if (j == second.size()) {
			sol.push_back(first[i]);
			++i;
		}
		else if (i == first.size()) {
			sol.push_back(second[j]);
			++j;
		}
		else if (first[i] < second[j]) {
			sol.push_back(first[i]);
			++i;
		}
		else {
			sol.push_back(second[j]);
			++j;
		}
	}
	return sol;
}

// code is "2sum"
// finding items in a set in c++ is logarithmic in terms of size of set
std::vector<std::pair<int, int>> sum2(std::vector<std::vector<int>> input) {
	std::vector<std::pair<int, int>> sol;
	for (int i = 0; i < input.size(); ++i) {
		std::set<int> temp;
		bool found = false;
		for (int j = 0; j < input[i].size(); ++j) {
			if (input[i][j] < 0) {
				auto it = temp.find(-input[i][j]);
				if (it != temp.end()) {
					int first = std::distance(input[i].begin(), std::find(input[i].begin(), input[i].end(), -input[i][j]));
					sol.push_back({first+1, j+1});
					found = true;
					break;
				}
				else {
					temp.insert(input[i][j]);
				}
			}
			else {
				auto it = temp.find(-input[i][j]);
				if (it != temp.end()) {
					int first = std::distance(input[i].begin(), std::find(input[i].begin(), input[i].end(), -input[i][j]));
					sol.push_back({ first + 1, j + 1 });
					found = true;
					break;
				}
				else {
					temp.insert(input[i][j]);
				}
			}
			
		}
		if (!found) {
			sol.push_back({ 0, 0 });
		}
		
	}
	return sol;
}

std::vector<int> bfs(std::multimap<int, int> adjList) {
	int nodes = adjList.rbegin()->first;
	std::vector<int> distFromFirst(nodes + 1, -1);
	std::vector<int> visited(nodes + 1, 0);
	std::queue<int> q;
	distFromFirst[1] = 0;
	q.push(1);
	visited[1] = 1;
	while (!q.empty()) {
		auto currNode = q.front();
		auto range = adjList.equal_range(q.front());
		q.pop();
		for (auto i = range.first; i != range.second; ++i) {
			if (!visited[i->second]) {
				q.push(i->second);
				distFromFirst[i->second] = distFromFirst[currNode]+1;
				visited[i->second] = 1;
			}
		}

	}
	return distFromFirst;
}


int cc(std::multimap<int, int> adjList) {
	return 0;
}

// Here we have to build a binary heap
std::vector<int> hea(vector<int> input) {
	vector<int> heap;
	heap.push_back(0);
	for (int i = 0; i < input.size(); ++i) {
		heap.push_back(input[i]);
		int j = heap.size() - 1;
		while (j > 1 && heap[j] > heap[j / 2]) {
			std::swap(heap[j], heap[j / 2]);
			j /= 2;
		}
	}
	return heap;
}
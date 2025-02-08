#include <vector>
#include <map>
#include <set>
#include <queue>
#include <deque>
#include <stack>
#include <iterator>
#include <chrono>


using namespace std::chrono;

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
int ins(std::vector<int> unsortedArray) {
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

template<typename T>
std::vector<T> insertionSort(std::vector<T> unsortedArray) {
	
	for (int i = 1; i < unsortedArray.size(); ++i) {
		int k = i;
		while (k > 0 && unsortedArray[k] < unsortedArray[k - 1]) {
			std::swap(unsortedArray[k], unsortedArray[k - 1]);
			
			--k;
		}
	}
	return unsortedArray;
}

// More generalized insertion sort
template<typename Iterator>
void insertionSort(Iterator itBegin, Iterator itEnd) {
	std::vector<typename Iterator::value_type> output;
	auto start = itBegin++;
	for (; itBegin != itEnd; ++itBegin) {
		auto currIt = itBegin;
		while (itBegin != start && *itBegin < *(itBegin - 1)) {
			std::iter_swap(itBegin, itBegin - 1);
			itBegin--;
		}
		itBegin = currIt;
	}
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
	std::vector<int> solution;
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

// Merge function for merge sort
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

template <typename Iterator>
std::vector<typename Iterator::value_type> merge(const Iterator it1begin, const Iterator it1end, const Iterator it2begin, const Iterator it2end) {
	std::vector<typename Iterator::value_type> merged(std::distance(it1begin, it1end) + std::distance(it2begin, it2end));
	auto it1 = it1begin;
	auto it2 = it2begin;
	while (it1 != it1end && it2 != it2end) {
		merged.push_back(std::move((*it1 <= *it2) ? *it1++ : *it2++));
	}
	merged.insert(merged.end(), it1, it1end);
	merged.insert(merged.end(), it2, it2end);
	return merged;
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


int cc(std::multimap<int, int> adjList, int nodeNumber) {
	std::vector<int> visited(nodeNumber+1, 0);
	std::stack<int> nodes;
	int islands = 0;
	for (int i = 1; i <= nodeNumber; ++i) {
		if (!visited[i]) {
			visited[i] = 1;
			nodes.push(i);
			while (!nodes.empty()) {
				auto currNode = nodes.top();
				auto range = adjList.equal_range(currNode);
				nodes.pop();
				for (auto it = range.first; it != range.second; ++it) {
					if (!visited[it->second])
						nodes.push(it->second);
					visited[it->second] = 1;
				}
			}

			++islands;
		}
	}
	return islands;
}

void maxHeapify(std::vector<int>& input, int top, int size) {
	int left = 2 * top + 1;
	int right = 2 * top + 2;
	int largest;
	if (left < size && input[left] > input[top])
		largest = left;
	else largest = top;
	if (right < size && input[right] > input[largest])
		largest = right;
	if (largest != top) {
		std::swap(input[top], input[largest]);
		maxHeapify(input, largest, size);
	}
}

// Here we have to build a binary heap
std::vector<int> hea(std::vector<int> input) {
	
	for (int i = input.size() / 2; i > 0; --i) {
		maxHeapify(input, i-1, input.size());
	}
	return input;
} 

// Merge sort
template <typename T>
std::vector<T> ms(std::vector<T> input) {
	
	if (input.size() <=32) {
		return insertionSort(input);
	}
	std::vector<T> left(input.begin(), input.begin() + input.size() / 2);
	std::vector<T> right(input.begin() + input.size() / 2, input.end());

	return mer(ms(left), ms(right));
}



template <typename T>
int par(std::vector<T> &input, int low, int high) {
	std::deque<T> output;
	T firstElement = input[low];
	int i = low;
	output.push_back(firstElement);
	for (int j = low + 1; j <= high; ++j) {
		if (input[j] <= firstElement) {
			output.push_front(input[j]);

				++i;
		}
		else {
			output.push_back(input[j]);
		}
	}
	input = std::vector<T>(output.begin(), output.end());
	return i;
}

std::vector<int> sum3(std::vector<int> input) {
	auto copy = input;
	std::sort(copy.begin(), copy.end());
	std::vector<int> output;
	for (int i = 0; i < copy.size() - 2; ++i) {
		int first = i;
		int second = i + 1;
		int third = copy.size() - 1;

		while (second < third) {
			if (copy[first] + copy[second] + copy[third] == 0) {
				auto it1 = std::find(input.begin(), input.end(), copy[first]);
				auto it2 = std::find(input.begin(), input.end(), copy[second]);
				auto it3 = std::find(input.begin(), input.end(), copy[third]);
				output.push_back(std::distance(input.begin(), it1)+1);
				output.push_back(std::distance(input.begin(), it2)+1);
				output.push_back(std::distance(input.begin(), it3)+1);
				std::sort(output.begin(), output.end());
				return output;
			}
			else if (copy[first] + copy[second] + copy[third] < 0) {
				second++;
			}
			else third--;
		}
	}
	output.push_back(-1);
	return output;
}

int bip(std::multimap<int, int> adjList, int nodeNumber) {
	// I will use the fact that you can colour a graph with only two colors
	// if and only if it is bipartite
	// with simple bfs we can "color all the nodes and check for valid coloring

	std::queue<int> q;
	std::vector<int> color(nodeNumber + 1, 0); // 0 is unvisited, 1 is red, -1 is blue
	q.push(1);
	color[1] = 1;
	while (!q.empty()) {
		int currNode = q.front();
		q.pop();
		auto range = adjList.equal_range(currNode);
		for (auto it = range.first; it != range.second; ++it) {
			
			if (!color[it->second]) {
				q.push(it->second);
				color[it->second] = -color[currNode];
			}
			else {
				if (color[it->second] == color[currNode]) return -1;
			}
		}
	}
	return 1;
}

int dag(std::multimap<int, int> adjList, int nodeNumber) {

	std::stack<int> nodes;
	std::vector<int> visited(nodeNumber + 1, 0);
	nodes.push(1);
	visited[1] = 1;
	while (!nodes.empty()) {
		int currNode = nodes.top();
		nodes.pop();
		auto range = adjList.equal_range(currNode);
		for (auto it = range.first; it != range.second; ++it) {
			if (visited[it->second])
				return -1;
			else {
				nodes.push(it->second);
				visited[it->second] = 1;
			}
		}
	}
	return 1;
}

std::vector<int> dij(std::multimap<int, std::pair<int, int>> adjList, int nodeNumber, int startNode) {
	std::vector<int> visited(nodeNumber, 0);
	std::vector<int> dist(nodeNumber, -1);
	

	auto cmp = [](std::pair<int, int> left, std::pair<int, int> right) {return left.second > right.second; };
	std::priority_queue <std::pair<int, int>, std::vector<std::pair<int, int>>, decltype(cmp)> dij_que(cmp);
	dij_que.push({ startNode, 0 });
	dist[startNode-1] = 0;
	while (!dij_que.empty()) {
		auto currNode = dij_que.top();
		dij_que.pop();
		if (visited[currNode.first-1]) continue;
		visited[currNode.first-1] = 1;

		auto range = adjList.equal_range(currNode.first);
		for (auto it = range.first; it != range.second; ++it) {
			if (!visited[it->second.first-1]) {
			
				int distance = dist[currNode.first-1] + it->second.second;
				if (distance < dist[it->second.first-1] || dist[it->second.first-1] == -1)
					dist[it->second.first-1] = distance;
				dij_que.push({ it->second.first, dist[it->second.first-1] });
			}
		}
	}

	return dist;
}

std::vector<int> hs(std::vector<int> input) {
	input = hea(input);
	for (int i = input.size() - 1; i > 0; --i) {
		std::swap(input[0], input[i]);
		maxHeapify(input, 0, i-1);
	}
	return input;
}

std::vector<int> bf(std::multimap<int, std::pair<int, int>> adjList, int nodeNumber) {
	std::vector<int> visited(nodeNumber + 1, 0);
	std::vector<int> dist(nodeNumber + 1, INT_MAX);
	dist[1] = 0;
	for (int i = 0; i < nodeNumber; ++i) {
		for (const auto& it : adjList) {
			int edge = it.second.second;
			int secondNode = it.second.first;
			if (dist[it.first] == INT_MAX) continue;
			if (dist[it.first] + edge < dist[secondNode]) {
				dist[secondNode] = dist[it.first] + edge;
			}
		}
	}
	return dist;
}

int nwc(std::multimap<int, std::pair<int, int>> adjList, int nodeNumber) {
	std::vector<int> visited(nodeNumber + 1, 0);
	std::vector<int> dist(nodeNumber + 1, INT_MAX);
	
	dist[0] = 0;
	for (int i = 1; i <= nodeNumber; ++i) {
		adjList.insert({ 0, {i, 0} });
	}
	for (int i = 0; i < nodeNumber; ++i) {
		for (const auto& it : adjList) {
			int edge = it.second.second;
			int secondNode = it.second.first;
			if (dist[it.first] == INT_MAX) continue;
			if (dist[it.first] + edge < dist[secondNode]) {
				dist[secondNode] = dist[it.first] + edge;
			}
		}
	}

	for (const auto& it : adjList) {
		if (!it.first) continue;
		int edge = it.second.second;
		int secondNode = it.second.first;
		if (dist[it.first] == INT_MAX) {
			continue;
		}
		if (dist[it.first] != INT_MAX == dist[secondNode] > dist[it.first] + edge) {
			return 1;
		}
	}
	return -1;
}

void dfsVisit(std::multimap<int, int> &adjList, std::vector<int>& visited, int& time, std::vector<int>& visitTime, std::vector<int>& endTime, int node, std::vector<int> &topologicallySorted) {
	++time;
	visitTime[node] = time;
	auto range = adjList.equal_range(node + 1);
	visited[node] = 1;
	for (auto it = range.first; it != range.second; ++it) {
		if (!visited[it->second - 1]) {
			dfsVisit(adjList, visited, time, visitTime, endTime, it->second-1, topologicallySorted);
		}
	}
	++time;
	endTime[node] = time;
	topologicallySorted.push_back(node+1);
}

std::vector<int> dfs(std::multimap<int, int> adjList, int nodeNumber, std::vector<int> &visitTime, std::vector<int> &endTime) {
	std::vector<int> visited(nodeNumber, 0);
	std::vector<int> topologicallySorted;
	int time = 0;
	for (int i = 0; i < nodeNumber; ++i) {
		if (!visited[i]) {
			dfsVisit(adjList, visited, time, visitTime, endTime, i, topologicallySorted);
		}
	}
	std::reverse(topologicallySorted.begin(), topologicallySorted.end());
	return topologicallySorted;
}

std::vector<int> ts(std::multimap<int, int> adjList, int nodeNumber) {
	std::vector<int> visitTime(nodeNumber);
	std::vector<int> endTime(nodeNumber);
	auto sorted = dfs(adjList, nodeNumber, visitTime, endTime);
	return sorted;
}

std::vector<int> par3(std::vector<int> input) {
	int first = input[0];
	std::deque<int> newList;
	std::vector<int> output;
	int same = 0;
	for (auto it : input) {
		if (it < first) newList.push_front(it);
		else if (it > first) newList.push_back(it);
		else ++same;
	}
	bool less = true;
	for (auto it : newList) {
		if (less && it > first) {
			less = false;
			for (int i = 0; i < same; ++i) {
				output.push_back(first);
			}
		}
		output.push_back(it);
	}
	return output;
}

int cte(std::multimap<int, std::pair<int, int>> adjList, int nodeNumber, int edgeStartNode, int edgeEndNode) {
	auto dist = dij(adjList, nodeNumber, edgeEndNode);
	if (dist[edgeStartNode - 1] == -1) return -1;
	else {
		auto range = adjList.equal_range(edgeStartNode);
		for (auto it = range.first; it != range.second; ++it) {
			if (it->second.first == edgeEndNode) {
				return dist[edgeStartNode - 1] + it->second.second;
			}
		}
	}
}

std::multimap<int, int> dist2Unity(std::multimap<int, std::pair<int, int>> &adjList) {
	std::multimap <int, int> undirected;
	for (auto it : adjList) {
		undirected.insert({ it.first, it.second.first });
	}
	return undirected;
}

std::vector<int> sdag(std::multimap<int, std::pair<int, int>> adjList, int nodeNumber) {
	std::vector<int> distances(nodeNumber, INT_MAX);
	std::vector<bool> visited(nodeNumber);
	std::vector<int> sorted = ts(dist2Unity(adjList), nodeNumber);
	distances[0] = 0;
	bool start = false;
	for (auto it : sorted) {
		if (it != 1 && !start)continue;
		if (it == 1) start = true;
		auto range = adjList.equal_range(it);
		
		for (auto it2 = range.first; it2 != range.second; ++it2) {
			int dist = distances[it - 1] + it2->second.second;
			if (dist < distances[it2->second.first-1]) {
				distances[it2->second.first-1] = dist;
			}
		}
	}
	return distances;
}

int inv(std::vector<int> inputArray)
{
	// BAD SOLUTION
	int inverseNumber = 0;
	for (int i = 0; i < inputArray.size(); ++i)
	{
		int currNum = inputArray[i];
		for (int j = i + 1; j < inputArray.size(); ++j)
		{
			if (currNum > inputArray[j])
			{
				inverseNumber++;
			}
		}
	}
	return inverseNumber;
}
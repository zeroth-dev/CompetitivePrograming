
#include <vector>
#include <map>
#include <queue>
#include <functional>
#include "BasicNode.hpp"
#ifndef BASICGRAPH_H
#define BASICGRAPH_H

template <typename T, typename I = int>
class BasicGraph 
{
#ifndef INF
#define INF std::numeric_limits<I>::infinity();
#endif
public:
	BasicGraph(int n);
	BasicGraph() {};

private:
	std::vector<BasicNode<T, I>> _nodes;
	std::multimap<BasicNode<T, I>, std::pair< BasicNode<T, I>, I>> _edges;
	std::map<T, int> _location; // needs testing against vectod<pair<T, int>>
	std::vector<std::pair<T, int>> _loc;
public:
	void addNode(T& nodeID);
	void removeNode(T& nodeID);
	BasicNode<T, I> getNode(T& nodeID);

	void addDirectedEdge(T& nodeID1, T& nodeID2, I distance);
	void addUndirectedEdge(T& nodeID1, T& nodeID2);
	std::vector<BasicNode<T, I>> nodes() { return _nodes; };
	void BFS(T nodeID);
	void DFS(T nodeID);
	int nodeNumber = 0;
	void printDistances();
};

template <typename T, typename I>
BasicGraph<T, I>::BasicGraph(int n) {
	_nodes.reserve(n);
}

template <typename T, typename I>
void BasicGraph<T, I>::addNode(T& nodeID) {
	_location.insert({ nodeID, _nodes.size() });
	_nodes.push_back(nodeID);
	++nodeNumber;
}


template <typename T, typename I>
void BasicGraph<T, I>::removeNode(T& nodeID) {
	_nodes.erase(_nodes.begin() + _location[nodeID]);
	_location.erase(_location.find(nodeID));
	--nodeNumber;
}

template <typename T, typename I>
BasicNode<T, I> BasicGraph<T, I>::getNode(T& nodeID) {
	return _nodes[_location[nodeID]];
}

template <typename T, typename I>
void BasicGraph<T, I>::addDirectedEdge(T& nodeID1, T& nodeID2, I distance) {
	if (_location.find(nodeID1) == _location.end()) addNode(nodeID1);
	if (_location.find(nodeID2) == _location.end()) addNode(nodeID2);

	_nodes[_location[nodeID1]].connectedNodes.push_back(_nodes[_location[nodeID2]]);
	_edges.insert({ nodeID1, {nodeID2, distance} });
}

template <typename T, typename I>
void BasicGraph<T, I>::addUndirectedEdge(T& nodeID1, T& nodeID2) {
	if (_location.find(nodeID1) == _location.end()) addNode(nodeID1);
	if (_location.find(nodeID2) == _location.end()) addNode(nodeID2);

	_nodes[_location[nodeID1]].connectedNodes.push_back(_nodes[_location[nodeID2]]);
	_nodes[_location[nodeID2]].connectedNodes.push_back(_nodes[_location[nodeID1]]);
	_edges.insert({ nodeID1, {nodeID2, 1.0} });
	_edges.insert({ nodeID2, {nodeID1, 1.0} });

}




template <typename T, typename I>
void BasicGraph<T, I>::BFS(T nodeID) {
	std::vector<bool> visited(nodeNumber, false);
	int loc = _location[nodeID];
	std::queue<int> q;
	q.push(loc);
	visited[loc] = true;
	_nodes[loc].distance = 0;
	while (!q.empty()) {
		auto currNode = q.front();
		q.pop();
		auto range = _edges.equal_range(currNode);
		for (auto it = range.first; it != range.second; ++it) {
			loc = _location[it->second.first.id()];
			if (!visited[loc]) {
				q.push(loc);
				_nodes[loc].distance = _nodes[currNode].distance + 1;
				visited[loc] = true;
			}
		}
	}
}

template <typename T, typename I>
void BasicGraph<T, I>::printDistances() {
	auto nodes = _nodes;
	std::sort(nodes.begin(), nodes.end(), 
		[]( BasicNode<T, I> lhs,  BasicNode<T, I> rhs) -> bool {
			return (lhs < rhs); 
		});
	for (auto it : nodes) {
		std::cout << "Distance to node " << it.id() << " is " << it.distance << std::endl;

	}
}

#endif // BASICGRAPH_H

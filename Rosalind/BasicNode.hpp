
#include <vector>
#include <limits>
#include <iostream>
#ifndef BASICNODE_H
#define BASICNODE_H
template <typename T, typename I = int>
class BasicNode
{

#define INF std::numeric_limits<I>::infinity();
public:

	BasicNode(T id) : _id(id) {};
private:
	T _id;
	T _parent = NULL;
	
public:
	T id() { return _id; }
	T parent() { return _parent; }
	I distance = INF;

	std::vector<BasicNode<T>> connectedNodes;

	bool operator==(const BasicNode& rhs) {
		return _id == rhs._id;
	}
};

#endif // BASICNODE_H

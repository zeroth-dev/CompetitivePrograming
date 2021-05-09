
#include <vector>
#include <climits>
#include <iostream>
#ifndef BASICNODE_H
#define BASICNODE_H
template <typename T, typename I = int>
class BasicNode
{

#define INF std::numeric_limits<I>::max()
public:

	BasicNode(T id) : _id(id) {};
private:
	T _id;
	T _parent = NULL;
	
	
public:
	T id() const { return _id; }
	T parent() { return _parent; }
	I distance = INF;
	int startTime = 0;
	int finishTime = 0;
	std::vector<BasicNode<T>> connectedNodes;

	bool operator==(const BasicNode& rhs) {
		return _id == rhs._id;
	}
	
};
template <typename T, typename I = int>
inline bool operator<(const BasicNode<T, I>& lhs, const BasicNode<T, I>& rhs) {
	return lhs.id() < rhs.id();
}
#endif // BASICNODE_H

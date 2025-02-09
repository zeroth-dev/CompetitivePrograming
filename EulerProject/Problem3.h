#pragma once
#include "GenericProblem.h"
#include <vector>
class Problem3 :
    public GenericProblem
{

private:
	std::vector<long long int> sieveData;
	void generateSieve(long long int maxNumber);
	long long int mNumber = 0;

public:
	Problem3() {};
	virtual void solve(std::ostream&);
	virtual void readInput(std::istream&) {}
	virtual void setNumber(long long int numberToSet) { mNumber = numberToSet; }
};


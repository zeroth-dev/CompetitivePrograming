#pragma once
#include "GenericProblem.h"
#include <vector>
class Problem5 :
    public GenericProblem
{
private:
	std::vector<long long int> sieveData;
	void generateModifiedSieve(long long int maxNumber);
	long long int mNumber = 0;

public:
	Problem5() {};
	virtual void solve(std::ostream&);
	virtual void readInput(std::istream&) {}
	virtual void setNumber(long long int numberToSet) { mNumber = numberToSet; }
};


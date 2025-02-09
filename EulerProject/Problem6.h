#pragma once
#include "GenericProblem.h"
class Problem6 :
    public GenericProblem
{
	long long int mNumber = 0;

public:
	Problem6() {};
	virtual void solve(std::ostream&);
	virtual void readInput(std::istream&) {}
	virtual void setNumber(long long int numberToSet) { mNumber = numberToSet; }
};


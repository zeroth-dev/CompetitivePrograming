#pragma once
#include "GenericProblem.h"
class Problem2 :
    public GenericProblem
{
public:
	Problem2() {};
	virtual void solve(std::ostream&);
	virtual void readInput(std::istream&) {}
	virtual void setNumber(long long int) {}
};


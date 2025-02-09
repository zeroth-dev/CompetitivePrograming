#pragma once
#include "GenericProblem.h"

class Problem1 : public GenericProblem
{
public:
	Problem1() {};
	virtual void solve(std::ostream&);
	virtual void readInput(std::istream&) {}
	virtual void setNumber(long long int) {}
};


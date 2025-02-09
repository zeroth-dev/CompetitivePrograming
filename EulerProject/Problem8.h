#pragma once
#include "GenericProblem.h"
class Problem8 :
    public GenericProblem
{
private:
	std::string inputNumber;
	int mNumber;
public:
	Problem8() {};
	virtual void solve(std::ostream&);
	virtual void readInput(std::istream&);
	virtual void setNumber(long long int numberToSet) { mNumber = numberToSet; }
};


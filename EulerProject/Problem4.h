#pragma once
#include "GenericProblem.h"
class Problem4 :
    public GenericProblem
{
private:
	bool checkPalindrome(int number);
public:
	Problem4() {};
	virtual void solve(std::ostream&);
	virtual void readInput(std::istream&) {}
	virtual void setNumber(long long int numberToSet) { }
};


#pragma once
#include <iostream>

// Abstract class for all the problems
class GenericProblem
{

public:

	// Set a single number required for calculation
	virtual void setNumber(long long int) = 0;

	// Read the input from the specified input stream
	virtual void readInput(std::istream&) = 0;

	// Solve the problem and output the result to the specified output stream
	virtual void solve(std::ostream&) = 0;
};


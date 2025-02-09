#include <iostream>
#include <fstream>
#include "Problem8.h"

int main()
{
	std::ifstream input;
	input.open("Data/problem8.txt");
	GenericProblem *problem = new Problem8();
	problem->readInput(input);
	input.close();
	problem->setNumber(13);
	problem->solve(std::cout);
	return 0;
}
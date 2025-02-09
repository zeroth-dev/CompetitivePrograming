#include "Problem8.h"

void Problem8::solve(std::ostream& output)
{
	long long int largestProduct = 0;
	for (int i = 0; i < inputNumber.size()-mNumber; i++)
	{
		long long int currProduct = 1;
		// We start from the last number in the line because if we reach 0 that will advance us the fastest
		for (int currNum = mNumber-1; currNum >= 0; currNum--)
		{
			currProduct *= static_cast<int>(inputNumber[i + currNum] - '0');
			// If we got a 0, advance until the next digit
			if (currProduct == 0)
			{
				i += currNum;
				break;
			}
			if (currProduct > largestProduct)
			{
				largestProduct = currProduct;
			}
		}
	}
	output << largestProduct;
}

void Problem8::readInput(std::istream& input)
{
	std::string currentLine;
	while (input >> currentLine)
	{
		inputNumber += currentLine;
	}
}

#include "Problem3.h"

void Problem3::generateSieve(long long int maxNumber)
{
	std::vector<bool> markedList(maxNumber + 1, true);
	for (long long int i = 2; i <= maxNumber; i++)
	{
		if (markedList[i] == true)
		{
			
			for (long long int temp = i+i; temp <= maxNumber; temp += i)
			{
				markedList[temp] = false;
			}
			sieveData.push_back(i);
		}
	}
}

void Problem3::solve(std::ostream& output)
{
	// Create a vector of all the prime numbers up to sqrt(mNumber)
	generateSieve(std::sqrt(mNumber));

	// Start from the end of the vector and check if it divides mNumber
	// Since we are iterating from the end, we will get the largest prime divisor
	for (auto it = sieveData.rbegin(); it != sieveData.rend(); ++it)
	{
		if (mNumber % *it == 0)
		{
			output << *it;
			return;
		}
	}
}

#include "Problem5.h"

void Problem5::generateModifiedSieve(long long int maxNumber)
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
			long long int largestPrimePower = i;
			while (largestPrimePower <= maxNumber)
			{
				largestPrimePower *= i;
			}
			sieveData.push_back(largestPrimePower/i);
		}
	}
}

void Problem5::solve(std::ostream& output)
{
	generateModifiedSieve(mNumber);
	long long int solution = 1;
	for (auto num : sieveData)
	{
		solution *= num;
	}
	output << solution;

}

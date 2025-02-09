#include "Problem7.h"

void Problem7::generateSieve(long long int maxNumber)
{
	std::vector<bool> markedList(maxNumber + 1, true);
	for (long long int i = 2; i <= maxNumber; i++)
	{
		if (markedList[i] == true)
		{

			for (long long int temp = i + i; temp <= maxNumber; temp += i)
			{
				markedList[temp] = false;
			}
			sieveData.push_back(i);
			std::cout << "Prime found: " << i << std::endl;
			if (sieveData.size() > mNumber)
			{
				return;
			}
		}
	}
}

void Problem7::solve(std::ostream& output)
{
	long long int largestNumber = mNumber * mNumber;
	generateSieve(largestNumber);

	while (sieveData.size() < mNumber)
	{
		sieveData.clear();
		largestNumber *= 10;
		generateSieve(largestNumber);
	}
	output << sieveData[mNumber-1];
}

#include "Problem2.h"

void Problem2::solve(std::ostream& output)
{
	long long int finalSum = 2;
	long long int fiboN2 = 1;
	long long int fiboN1 = 2;
	long long int currFibo = 0;
	while (currFibo < 4e6)
	{
		currFibo = fiboN2 + fiboN1;
		fiboN2 = fiboN1;
		fiboN1 = currFibo;
		if (currFibo % 2 == 0)
		{
			finalSum += currFibo;
		}
	}
	output << finalSum;
}

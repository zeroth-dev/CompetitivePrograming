#include "Problem1.h"


void Problem1::solve(std::ostream& output)
{
	// The sum of all numbers less than a 100
	// that can be divided by 3 or 5 can be calculated directly
	// First we find the sum off all multiples of 3 and 5 less than 1000
	// Then we substract the sum of all multiples of 3*5=15 less than 100
	int sum3 = (999 / 3) * (3 + 999) / 2;
	int sum5 = (995 / 5) * (5 + 995) / 2;
	int sum15 = (990 / 15) * (15 + 990) / 2;

	int finalSum = sum3 + sum5 - sum15;
	output << finalSum << std::endl;
}

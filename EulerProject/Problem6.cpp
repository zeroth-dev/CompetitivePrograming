#include "Problem6.h"

void Problem6::solve(std::ostream& output)
{
	long long int sumOfNums = (mNumber * (mNumber + 1) / 2);
	long long sumOfNumsSquared = sumOfNums* sumOfNums;

	long long int sumOfSquaredNums = mNumber * (2 * mNumber + 1) * (mNumber + 1) / 6;
	output << sumOfSquaredNums << " " << sumOfNumsSquared << std::endl;
	output << sumOfNumsSquared - sumOfSquaredNums;

}

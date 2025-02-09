#include "Problem4.h"

bool Problem4::checkPalindrome(int number)
{
	// We don't need to construct a string in the correct order 
	// since we are checking for palindromes which are still the same when reversed
	std::string strNumber;
	while (number > 0)
	{
		char currDigit = '0' + (number % 10);
		strNumber.push_back(currDigit);
		number /= 10;
	}
	bool isPalindrome = true;
	for (int i = 0; i < strNumber.size() / 2; i++)
	{
		if (strNumber[i] != strNumber[strNumber.size() - 1 - i])
		{
			isPalindrome = false;
			break;
		}
	}
	return isPalindrome;
}

void Problem4::solve(std::ostream& output)
{
	// To find the largest palindrome you can get by multiplying 3-digit numbers
	// we first assume we will get a 6-digit number
	// Then we assume that the two 3-digit numbers are of a form 9xy
	// We check for each combination of two 2-digit numbers xy to find the largest palindrome
	// We start from the largest towards the lowest numbers
	for (int firstNum = 99; firstNum >= 10; firstNum--)
	{
		for (int secNum = firstNum; secNum >= 10; secNum--)
		{
			long potentialPalindrome = (900 + firstNum) * (900 + secNum);
			if (checkPalindrome(potentialPalindrome))
			{
				output << potentialPalindrome << " " << 900+firstNum << " " << 900+secNum << std::endl;
			}
		}
	}
}

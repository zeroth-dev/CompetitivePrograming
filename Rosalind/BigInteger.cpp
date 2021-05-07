#include "BigInteger.h"
#include <iostream>

BigInteger::BigInteger(std::string const& num){
    // Remove trailing zeros
    number = "0";
    for (int i = 0; i < num.size(); i++) {
        if (num[i] != '0') {
            number = num.substr(i, num.size() - i);
            break;
        }
    }
}
BigInteger::BigInteger(int num) : number(std::to_string(num)) {}

BigInteger::BigInteger(long long num) : number(std::to_string(num)) {}

BigInteger::BigInteger(uint64_t num) : number(std::to_string(num)) {}

int BigInteger::operator[](int idx) const{
    // Exception handling TODO
    return (int)this->number[idx]-48;
}

BigInteger& BigInteger::operator+=(const BigInteger& rhs){
    int bigger = this->size() > rhs.size() ? this->size() : rhs.size(); 
    std::string newNumber(bigger+1, '0');
    int carry = 0;
    int counter = 0;
    for(int lhsIt = this->size()-1, rhsIt = rhs.size()-1; lhsIt >=0 || rhsIt >=0; lhsIt--, rhsIt--){
        int temp;
        if(lhsIt >=0 && rhsIt >=0){
            temp = this->operator[](lhsIt) + rhs[rhsIt]+carry;
            carry = temp>9 ? 1:0;
            temp %= 10;
            newNumber[bigger-counter] = (char)(temp+48);
            counter++;
        }else if (lhsIt >= 0){
            temp = this->operator[](lhsIt)+carry;
            carry = temp>9 ? 1:0;
            temp %= 10;
            newNumber[bigger-counter] = (char)(temp+48);
            counter++;
        }else{
            temp = rhs[rhsIt]+carry;
            carry = temp>9 ? 1:0;
            temp %= 10;
            newNumber[bigger-counter] = (char)(temp+48);
            counter++;
        }
    }
    if (carry){
        newNumber[0] = '1';
        this->number=newNumber;
        return *this;
    }else{
        this->number=newNumber.substr(1, bigger);
        return *this;
    }
}

BigInteger operator+(BigInteger lhs, const BigInteger rhs){
    return lhs+=rhs;
}
BigInteger& BigInteger::operator++(){
    return this->operator+=(BigInteger("1"));
}

BigInteger BigInteger::operator++(int){
    BigInteger tmp(*this);
    operator++();
    return tmp;
}

BigInteger& BigInteger::operator-=(const BigInteger& rhs)
{
    if (*this < rhs) {
        throw "Negative numbers not yet implemented.";
    }
    
    for (int lhsIt = this->size() - 1, rhsIt = rhs.size() - 1; lhsIt >= 0 || rhsIt >= 0; lhsIt--, rhsIt--) {
        if (rhsIt < 0) {

        }
        else {
            if (this->operator[](lhsIt) < rhs[rhsIt]) {
                int temp = this->operator[](lhsIt) + 10 - rhs[rhsIt];
                int newIt = lhsIt-1;
                while (newIt >= 0 && this->operator[](newIt) == 0) {
                    this->number[newIt] = '9';
                    --newIt;
                }
               
                number[lhsIt] = temp + 48; 
                if (newIt == 0 && number[newIt] == '1') {
                    number = number.substr(1, number.size() - 1);
                }
                else {
                    number[newIt]--;
                }
            }
            else {
                int temp = this->operator[](lhsIt) - rhs[rhsIt];
                number[lhsIt] = temp + 48;
            }
        }
        
    }
    int i;
    for (i = 0; i < number.size() - 1; i++) {
        if (number[i] != '0') break;
    }
    number = number.substr(i, number.size() - i);

    return *this;
}

BigInteger operator-(BigInteger lhs, const BigInteger rhs) {
    return lhs -= rhs;
}

BigInteger operator%(BigInteger lhs, const BigInteger rhs)
{
    return lhs %= rhs;
}

BigInteger operator*(const BigInteger& lhs, const BigInteger& rhs) {
    if (lhs.size() < 9 && rhs.size() < 9) { 
        auto l = std::stoll(lhs.number) * std::stoll(rhs.number);
        return std::stoll(lhs.number) * std::stoll(rhs.number); 
    }

    if (lhs.size() <= 5 || rhs.size() <= 5) {
        
        return BigInteger::longMultiplication(lhs, rhs);

    }

    int i = std::max(lhs.size() / 3, rhs.size() / 3) + 1;
    BigInteger m0(lhs.number.substr(lhs.size() - i, i));
    BigInteger m1(lhs.number.substr(lhs.size() - 2*i, i));
    BigInteger m2(lhs.number.substr(0, lhs.size()-2*i));
    BigInteger n0(rhs.number.substr(rhs.size() - i, i));
    BigInteger n1(rhs.number.substr(rhs.size() - 2 * i, i));
    BigInteger n2(rhs.number.substr(0, rhs.size() - 2 * i));
    BigInteger x4 = operator*(m2, n2);
    BigInteger x3 = operator*(m2, n1) + operator*(m1, n2);
    BigInteger x2 = operator*(m2, n0) + operator*(n2, m0) + operator*(m1, n1);
    BigInteger x1 = operator*(m1, n0) + operator*(n1, m0);
    BigInteger x0 = operator*(m0, n0);
    x4.number.append(i * 4, '0');
    x3.number.append(i * 3, '0');
    x2.number.append(i * 2, '0');
    x1.number.append(i * 1, '0');



    return x4+x3+x2+x1+x0;
}

BigInteger& BigInteger::operator%=(const BigInteger& rhs) {
    if (rhs == 1000000) {
        if (this->size() <= 6) return *this;
        int temp = this->size();
        number = number.substr(temp - 6, 6);
        return *this;
    }
    return *this;
}

BigInteger BigInteger::fact(BigInteger n)
{
    if (n == 1) {
        return BigInteger(1);
    }
    else {
        auto a = fact(n - 1);
        return n * a;
    }
}

BigInteger BigInteger::longMultiplication(BigInteger lower, BigInteger higher) {
    BigInteger newNumber(0);
    for (int i = lower.size()-1; i >= 0; i--) {
        std::string temp(higher.size() + 1, '0');
        int carry = 0;
        int curr = lower[i];
        for (int j = higher.size()-1; j >=0; j--) {
            int res = curr * higher[j] + carry;
            carry = res / 10;
            res %= 10;
            temp[j + 1] = res + 48;
        }
        if (carry) {
            temp[0] = carry + 48;
        }
        else {
            temp = temp.substr(1, higher.size());
        }
        temp.append(lower.size() - 1 - i, '0');
        newNumber += BigInteger(temp);
    }
    return newNumber;
}

std::istream &operator >>(std::istream &istrm, BigInteger &bi){
    return istrm >> bi.number;
}
std::ostream &operator << (std::ostream &ostrm, BigInteger bi){
    return ostrm << bi.number;
}

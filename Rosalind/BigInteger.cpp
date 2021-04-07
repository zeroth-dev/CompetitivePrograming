#include "BigInteger.h"
#include <iostream>

BigInteger::BigInteger(std::string const& num) : number(num){}

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


BigInteger operator*(const BigInteger &lhs, const BigInteger &rhs){
    
}

std::istream &operator >>(std::istream &istrm, BigInteger &bi){
    return istrm >> bi.number;
}
std::ostream &operator << (std::ostream &ostrm, BigInteger bi){
    return ostrm << bi.number;
}
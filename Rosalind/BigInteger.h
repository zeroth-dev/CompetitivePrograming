
#ifndef BIGINTEGER_H
#define BIGINTEGER_H

#include <string>
#include <istream>
#include <ostream>


class BigInteger{
public:
    BigInteger(std::string const& num);
    BigInteger(int num);
private:
    std::string number;
public:
    static BigInteger longMultiplication(BigInteger lower, BigInteger higher);
public:

    friend BigInteger operator* (const BigInteger &lhs, const BigInteger &);
    BigInteger&  operator +=(const BigInteger& rhs);
    BigInteger& operator -=(const BigInteger& rhs);
    BigInteger& operator ++();
    BigInteger operator ++(int);
    BigInteger& operator %=(const BigInteger& rhs);
    int operator[](int idx) const;
    friend BigInteger operator+(BigInteger lhs, const BigInteger rhs);
    friend BigInteger operator-(BigInteger lhs, const BigInteger rhs);
    friend std::istream &operator >> (std::istream &strm, BigInteger &bi);
    friend std::ostream &operator << (std::ostream &strm, BigInteger bi);
   
    int size() const { return this->number.size(); }
    static BigInteger fact(BigInteger n);
};
BigInteger operator+(BigInteger lhs, const BigInteger rhs);
BigInteger operator*(const BigInteger &lhs, const BigInteger &rhs);
BigInteger operator-(BigInteger lhs, const BigInteger rhs);
BigInteger operator%(BigInteger lhs, const BigInteger rhs);
std::istream &operator >> (std::istream &istrm, BigInteger &bi);

std::ostream &operator << (std::ostream &ostrm, BigInteger bi);

inline bool operator ==(const BigInteger& lhs, const BigInteger& rhs) {
    if (lhs.size() != rhs.size()) return false;
    for (int i = 0; i < lhs.size(); i++) {
        if (lhs[i] != rhs[i]) return false;
    }
    return true;
}
inline bool operator!=(const BigInteger& lhs, const BigInteger& rhs) { return !operator==(lhs, rhs); }
inline bool operator<(const BigInteger& lhs, const BigInteger& rhs) {
    if (lhs.size() < rhs.size()) return true;
    if (lhs.size() > rhs.size()) return false;
    if (lhs == rhs) return false;
    for (int i = 0; i < lhs.size(); i++) {
        if (lhs[i] > rhs[i]) return false;
    }
    
    return true;
}
inline bool operator>(const BigInteger& lhs, const BigInteger& rhs) { return operator<(rhs, lhs);}
inline bool operator<=(const BigInteger& lhs, const BigInteger& rhs) { return !operator>(lhs, rhs); }
inline bool operator>=(const BigInteger& lhs, const BigInteger& rhs) { return !operator<(lhs, rhs); }


#endif // !BIGINTEGER_H
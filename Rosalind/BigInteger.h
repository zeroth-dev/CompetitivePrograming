
#include <string>
#include <istream>
#include <ostream>


class BigInteger{
public:
    BigInteger(std::string const& num);
private:
    std::string number;
public:

    friend BigInteger operator * (const BigInteger &lhs, const BigInteger &);
    BigInteger&  operator +=(const BigInteger& rhs);
    BigInteger& operator ++();
    BigInteger operator ++(int);

    int operator[](int idx) const;
    friend BigInteger operator+(BigInteger lhs, const BigInteger rhs);
    friend std::istream &operator >> (std::istream &strm, BigInteger &bi);
    friend std::ostream &operator << (std::ostream &strm, BigInteger bi);
    int size() const { return this->number.size(); }
   
};
BigInteger operator+(BigInteger lhs, const BigInteger rhs);
BigInteger operator*(const BigInteger &lhs, const BigInteger &rhs);
std::istream &operator >> (std::istream &istrm, BigInteger &bi);

std::ostream &operator << (std::ostream &ostrm, BigInteger bi);

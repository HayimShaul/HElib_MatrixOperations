#ifndef ___WORD___
#define ___WORD___

#include <assert.h>
#include <bitset>
#include <vector>

// Author: Hayim Shaul 2016.

// UnsignedWord.
// implements a positive number by its bit representation. Bit can be either PlainBit or EncryptedBit
// UnsignedWords are implemented to have a variable number of bits, with a maximum defined to be MAX_BIT_NUM
//
// Note: UnsignedWord can be only be positive, negative numbers can be implemented with 2-complement but require
// a fixed number of bit representation, or some other TIHKUMATION like keeping the sign in a separate bit
// and for every + operation compute both +/- and choose by multiplying with sign bit

template<class Bit>
class UnsignedWord {
private:
	static const int MAX_BIT_NUM = 3;

	std::vector<Bit> _bits;
public:
	UnsignedWord(int n = 0);
	UnsignedWord(const Bit &b);

	void operator+=(const UnsignedWord &w);
	void operator-=(const UnsignedWord &w);
	void operator*=(const UnsignedWord &w) { *this = *this * w; }
	void operator*=(const Bit &b) { *this = *this * b; }

	void operator<<=(int i);
	void operator>>=(int i);

	Bit &operator[](int i) { return _bits[i]; }
	const Bit &operator[](int i) const { return _bits[i]; }

	//neg should expand the bit representation to MAX_BIT_NUM and use 2-complement
	//void neg();

	UnsignedWord<Bit> operator+(const UnsignedWord<Bit> &b) const { UnsignedWord<Bit> c(*this); c+=b; return c; }
	UnsignedWord<Bit> operator-(const UnsignedWord<Bit> &b) const { UnsignedWord<Bit> c(*this); c-=b; return c; }
	UnsignedWord<Bit> operator*(const UnsignedWord<Bit> &b) const;
	UnsignedWord<Bit> operator*(const Bit &b) const;
	UnsignedWord<Bit> operator<<(int i) const { UnsignedWord<Bit> c(*this); c<<=i; return c; }
	UnsignedWord<Bit> operator>>(int i) const { UnsignedWord<Bit> c(*this); c>>=i; return c; }

	Bit operator>(const UnsignedWord<Bit> &b) const;
	Bit operator<(const UnsignedWord<Bit> &b) const { return b.operator>(*this); }

	void setBitLength(int i) { _bits.resize(i); }
	int bitLength() const { return _bits.size(); }
	unsigned long get() const;
	std::ostream &operator<<(std::ostream &out);

	template<int A>
	std::bitset<A> to_bitset() const;

	static int in_range(int a) { return a & (((unsigned int)-1) >> (32 - MAX_BIT_NUM)); }
	static int max_bit_num() { return MAX_BIT_NUM; }
};

template<class Bit>
const int UnsignedWord<Bit>::MAX_BIT_NUM;


// implementation

template<class Bit>
inline UnsignedWord<Bit>::UnsignedWord(int n) {
	while ((n > 0) && (_bits.size() < MAX_BIT_NUM)) {
		_bits.push_back(Bit(n & 1));
		n >>= 1;
	}
}

template<class Bit>
inline void UnsignedWord<Bit>::operator+=(const UnsignedWord<Bit> &w) {
	int len = std::max(bitLength(), w.bitLength());
	int min = std::min(bitLength(), w.bitLength());

	int i;
	Bit carry(0);

	for (i = 0; i < min; ++i) {
		Bit new_carry = w[i]*carry + _bits[i]*(w[i] + carry);
		_bits[i] += w[i] + carry;
		carry = new_carry;
	}

	for (i = min; i < bitLength(); ++i) {
		Bit new_carry =  _bits[i]*carry;
		_bits[i] += carry;
		carry = new_carry;
	}

	_bits.resize(std::min(len + 1, MAX_BIT_NUM));

	for (i = min; i < w.bitLength(); ++i) {
		Bit new_carry = w[i]*carry;
		_bits[i] = w[i] + carry;
		carry = new_carry;
	}

	if (len < MAX_BIT_NUM)
		_bits[len] = carry;

}

template<class Bit>
inline void UnsignedWord<Bit>::operator-=(const UnsignedWord<Bit> &w) {
	assert(w.bitLength() <= bitLength());

	int min = std::min(bitLength(), w.bitLength());

	int i;
	Bit borrow(0);

	for (i = 0; i < min; ++i) {
		Bit new_borrow = (!_bits[i])*(w[i] + borrow) + w[i]*borrow;
		_bits[i] += w[i] + borrow;
		borrow = new_borrow;
	}

	for (i = min; i < bitLength(); ++i) {
		Bit new_borrow = (!_bits[i])*borrow;
		_bits[i] += borrow;
		borrow = new_borrow;
	}
}

template<class Bit>
void UnsignedWord<Bit>::operator>>=(int i) {
	if (i >= bitLength()) {
		_bits.resize(1);
		_bits[0] = Bit(0);
		return;
	}

	int j;
	for (j = 0; j < bitLength() - i; ++j)
		_bits[j] = _bits[j + i];
	_bits.resize(bitLength() - i);
}

template<class Bit>
void UnsignedWord<Bit>::operator<<=(int i) {
	int j;

	_bits.resize(std::min(MAX_BIT_NUM, bitLength() + i));

	for (j = bitLength() - 1; j >= i; --j)
		_bits[j] = _bits[j - i];
	for (j = i - 1; j >= 0; --j)
		_bits[j] = Bit(0);
}

template<class Bit>
UnsignedWord<Bit> UnsignedWord<Bit>::operator*(const Bit &b) const {
	UnsignedWord<Bit> a(*this);
	for (int i = 0; i < bitLength(); ++i)
		a._bits[i] *= b;
	return a;
}

template<class Bit>
UnsignedWord<Bit> operator*(const Bit &b, const UnsignedWord<Bit> &w) { return w*b; }


template<class Bit>
UnsignedWord<Bit> UnsignedWord<Bit>::operator*(const UnsignedWord<Bit> &w) const {
	if (bitLength() > w.bitLength())
		return w * (*this);

	std::vector< UnsignedWord<Bit> > toAdd;
	toAdd.resize(w.bitLength());

	toAdd[0] = (*this) * w[0];

	for (int i = 1; i < w.bitLength(); ++i)
		toAdd[i] = (*this << i) * w[i];

	return addArray(toAdd);
}

template<class Bit>
unsigned long UnsignedWord<Bit>::get() const {
	unsigned long r = 0;
	for (int i = 0; i < bitLength(); ++i)
		r += _bits[i].get() << i;
	return r;
}

template<class Bit>
template<int A>
std::bitset<A> UnsignedWord<Bit>::to_bitset() const {
	std::bitset<A> b;
	unsigned long r = 0;
	for (int i = 0; i < bitLength(); ++i)
		b[i] = _bits[i].get();
	return b;
}

// compute the multiplication of every bit arr[i] for    start <= i < end
template<class Bit>
Bit mulArray(const std::vector<Bit> &arr, int start = 0, int end = -1) {
	if (end == -1)
		end = arr.size();

	if (start == end - 1)
		return arr[start];

	if (start == end - 2)
		return arr[start] * arr[start + 1];

	int mid = (start + end - 1) / 2;

	return mulArray(arr, start, mid) * mulArray(arr, mid, end);
}

// compute the and of every bit arr[i] for    start <= i < end
template<class Bit>
Bit addArray(const std::vector<Bit> &arr, int start = 0, int end = -1) {
	if (end == -1)
		end = arr.size();

	if (start == end - 1)
		return arr[start];

	if (start == end - 2)
		return arr[start] + arr[start + 1];

	int mid = (start + end - 1) / 2;

	return addArray(arr, start, mid) + addArray(arr, mid, end);
}

template<class Bit>
Bit UnsignedWord<Bit>::operator>(const UnsignedWord<Bit> &w) const {
	int i;

	if (bitLength() == 0)
		return Bit(0);
	if (w.bitLength() == 0)
		return Bit(1);

	Bit ourMsbIsAllZero(1);
	Bit hisMsbIsAllZero(1);

	//	compute this loop in a shallower circuit
	//	for (i = bitLength() - 1; i >= w.bitLength(); --i)
	//		ourMsbIsAllZero *= !((*this)[i]);
	if (bitLength() > w.bitLength()) {
		std::vector<Bit> arr;
		arr.resize(bitLength() - w.bitLength());
		for (int j = w.bitLength(); j < bitLength(); ++j)
			arr[j - w.bitLength()] = !((*this)[j]);
		ourMsbIsAllZero = mulArray(arr);
	}

	//	compute this loop in a shallower circuit
	//	for (i = w.bitLength() - 1; i >= bitLength(); --i)
	//		hisMsbIsAllZero *= !(w[i]);
	if (w.bitLength() > bitLength()) {
		std::vector<Bit> arr;
		arr.resize(w.bitLength() - bitLength());
		for (int j = bitLength(); j < w.bitLength(); ++j)
			arr[j - bitLength()] = !(w[j]);
		hisMsbIsAllZero = mulArray(arr);
	}

	int commonLength = std::min(bitLength(), w.bitLength());

	std::vector<Bit> sameBit;
	sameBit.resize(commonLength);
	for (i = 0; i < commonLength; ++i)
		sameBit[i] = !((*this)[i] + w[i]);

	std::vector<Bit> equalMsbVec;
	equalMsbVec.resize(commonLength);
	equalMsbVec[commonLength - 1] = Bit(1);
	for (i = 0; i < commonLength - 1; ++i)
		equalMsbVec[i] = mulArray(sameBit, i+1, -1);

// TODO: bodyIsGreaterVec should be bitDeterminesGreater or something
	std::vector<Bit> bodyIsGreaterVec;
	bodyIsGreaterVec.resize(commonLength);
	for (i = 0; i < commonLength - 1; ++i)
		bodyIsGreaterVec[i] = equalMsbVec[i]*((*this)[i])*(!w[i]);
	// We know equalMsbVec[commonLength-1] == Bit(1)
	bodyIsGreaterVec[commonLength - 1] = ((*this)[commonLength - 1])*(!w[commonLength - 1]);
	Bit bodyIsGreater = addArray(bodyIsGreaterVec);

//	Bit equalMSB(1);
//	Bit bodyIsGreater(0);
//	for (i = commonLength - 1; i >= 0; --i) {
////		bodyIsGreater = bodyIsGreater | (equalMSB*((*this)[i])*(!w[i]));
////		which equals
////		bodyIsGreater = bodyIsGreater + (equalMSB*((*this)[i])*(!w[i])) + bodyIsGreater*(equalMSB*((*this)[i])*(!w[i]));
////		but bodyIsGreater==1 =>  equalMSB==0,    and equalMSB==1 => bodyIsGreater==0,   so we can write:
////		bodyIsGreater = bodyIsGreater + (equalMSB*((*this)[i])*(!w[i]));
////		equalMSB *= !(((*this)[i]) + w[i]);
//		bodyIsGreater = bodyIsGreater + (equalMsbVec[i]*((*this)[i])*(!w[i]));
//	}

//	Bit isGreater = (!ourMsbIsAllZero) | (hisMsbIsAllZero*bodyIsGreater);
//	Bit isGreater = (!ourMsbIsAllZero) + (hisMsbIsAllZero*bodyIsGreater) + (!ourMsbIsAllZero)*hisMsbIsAllZero*bodyIsGreater;
//	Since ourMsbAllZero==0  => hisMsbIsAllZero==1  we can write
//	Bit isGreater = (!ourMsbIsAllZero) + (hisMsbIsAllZero*bodyIsGreater) + (!ourMsbIsAllZero)*bodyIsGreater;
//	but that is really:
//	Bit isGreater = (Bit(1) + ourMsbIsAllZero) + (hisMsbIsAllZero*bodyIsGreater) + (Bit(1) + ourMsbIsAllZero)*bodyIsGreater;
//	Bit isGreater = Bit(1) + ourMsbIsAllZero + (hisMsbIsAllZero*bodyIsGreater) + (ourMsbIsAllZero*bodyIsGreater) + bodyIsGreater;
//	Bit isGreater = Bit(1) + ourMsbIsAllZero + (Bit(1)+hisMsbIsAllZero+ourMsbIsAllZero)*bodyIsGreater;
	Bit isGreater = Bit(1) + ourMsbIsAllZero + (Bit(1)+hisMsbIsAllZero+ourMsbIsAllZero)*bodyIsGreater;
//	Bit isGreater = Bit(1) + ourMsbIsAllZero + (!(hisMsbIsAllZero+ourMsbIsAllZero))*bodyIsGreater;

	return isGreater;
}


template<class Bit>
UnsignedWord<Bit> min(const UnsignedWord<Bit> &a, const UnsignedWord<Bit> &b) {
	Bit bIsGreater = a < b;
	return bIsGreater*a + (Bit(1)+bIsGreater)*b;
}

template<class Bit>
UnsignedWord<Bit> min(const std::vector< UnsignedWord<Bit> > &arr, int start = 0, int end = -1) {
	if (end == -1)
		end = arr.size();

	if (start == end - 1)
		return arr[start];

	if (start == end - 2)
		return min(arr[start], arr[start + 1]);

	int mid = (start + end - 1) / 2;

	return min(min(arr, start, mid), min(arr, mid, end));
}

template<class Bit>
std::ostream &operator<<(std::ostream &out, const UnsignedWord<Bit> &w) {
	for (int i = w.bitLength() - 1; i >= 0; --i)
		out << w[i];
	return out;
}

#endif


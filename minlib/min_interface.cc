#include <time.h>
#include <stdlib.h>

#include <iostream>
#include <sstream>
#include <bitset>

#include "plain_bit.h"
#include "encrypted_bit.h"
#include "unsigned_word.h"
#include "polynomial.h"
#include "converter.h"
#include "print.h"

ZP zp;

#define DEBUG

#ifdef DEBUG
#	define EncryptedNumber int
#endif


EncryptedNumber doMin(const std::vector<EncryptedNumber> &numbers) {
	SumOfMulOfPolynomialTerms<int> _p;
	SumOfMulOfPolynomialTerms<EncryptedNumber> p;

	std::cerr << "constructing polynomial" << std::endl;
	_p.constructMinPolynomial(numbers.size(), zp);
	std::cerr << "polynomial constructed" << std::endl;
	p.encrypt(_p);

	EncryptedNumber min = p.compute(numbers);

	return min;
}

EncryptedNumber doMin2(const std::vector<EncryptedNumber> &numbers, int begin = 0, int end = -1) {
	static SumOfMulOfPolynomialTerms<int> _p;
	static SumOfMulOfPolynomialTerms<EncryptedNumber> p;

	if (end == -1) {
		_p.constructMinPolynomial(2, zp);
		p.encrypt(_p);
		end = numbers.size();
	}

	if (begin + 1 == end)
		return numbers[begin];
	if (begin + 2 == end) {
		std::vector<EncryptedNumber> _n;
		_n.resize(2);
		_n[0] = numbers[begin];
		_n[1] = numbers[begin + 1];
		EncryptedNumber r = p.compute(_n, NULL, 0, -1, &zp);
		return r;
	}

	int mid = (end - begin) / 2;

	std::vector<EncryptedNumber> _n;
	_n.resize(2);
	_n[0] = doMin2(numbers, begin, mid);
	_n[1] = doMin2(numbers, mid, end);

	return p.compute(_n, NULL, 0, -1, &zp);
}


int main(int, char **) {
	long R=1;
	long p = 31;
	long r = 1;
	long d = 1;
	long c = 2;
	long k=80;
	long L = 10;
	long s=0;
	long chosen_m = 0;
	Vec<long> gens;
	Vec<long> ords;

	// We need to initialize two structures.
	// Keys is an object that hold the keys for encrypting
	// Zp holds varius table related to Z_p - the fields over which we operate
	// 
	// These structs are init-ed here, but it is logical that you'd want to init
	// them in your code
	//
	// Init the FHE keys with all the parameters
#	ifndef DEBUG
	Keys::initKeys(s, R, p, r, d, c, k, 64, L, chosen_m, gens, ords);
#	endif
	// construct tables that depend on the prime that defines Z_p
	zp.set_p(p);

	std::vector<EncryptedNumber> numbers;
	numbers.resize(4);
	numbers[0] = Converter<EncryptedNumber>::fromInt(3);
	numbers[1] = Converter<EncryptedNumber>::fromInt(5);
	numbers[2] = Converter<EncryptedNumber>::fromInt(1);
	numbers[3] = Converter<EncryptedNumber>::fromInt(7);

	EncryptedNumber min = doMin2(numbers);

	std::cout << "Min = " << zp.mod(Converter<EncryptedNumber>::toInt(min)) << std::endl;
}

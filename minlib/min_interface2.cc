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

//#define DEBUG

#ifdef DEBUG
	typedef int Number;
	typedef PlainBit Bit;
#else
	typedef EncryptedNumber Number;
	typedef EncryptedBit Bit;
#endif

clock_t start;


UnsignedWord<Bit> min(const std::vector<UnsignedWord<Bit> > &numbers, int begin = 0, int end = -1) {
	if (end == -1)
		end = numbers.size();

	if (begin + 1 == end)
		return numbers[begin];
	if (begin + 2 == end) {
		UnsignedWord<Bit> r = min(numbers[begin], numbers[begin+1]);
		return r;
	}

	int mid = (end - begin) / 2;

	UnsignedWord<Bit> r1 = min(numbers, begin, mid);
	UnsignedWord<Bit> r2 = min(numbers, mid, end);

	return min(r1, r2);
}

int mylog2(int x) {
	for (int i = 1; i < 32; ++i)
		if (1 << i > x)
			return i;
	return -1;
}

template<class Bit, class Number>
void bitify(UnsignedWord<Bit> &out, const Number &in) {
	static std::vector< MulOfPolynomialTerms<Number> > bitPolynomials;

	if (bitPolynomials.size() == 0) {
		bitPolynomials.resize(mylog2(zp.p()));
		for (int i = 0; i < mylog2(zp.p()); ++i) {
			MulOfPolynomialTerms<int> p;
			std::vector<int> bits;

			for (int j = 1; j < zp.p(); ++j) {
				if (j & (1 << i))
					bits.push_back(j);
			}
			p.semiCharacteristicTerm(0, bits, zp);
			bitPolynomials[i].encrypt(p);
		}
	}


	out.setBitLength(bitPolynomials.size());
	for (int i = 0; i < bitPolynomials.size(); ++i) {
std::cerr << "Need to take power ^{zp.p()-1} of this\n";
exit(1);
		out[i] = bitPolynomials[i].compute(in);
	}
}

template<class Bit, class Number>
void unbitify(Number &out, const UnsignedWord<Bit> in) {
	out = in[0];
	for (int i = 1; i < in.bitSize(); ++i)
		out += in[i] * (1 << i);
}

UnsignedWord<Bit> doMin3(const std::vector<Number> &numbers) {
	std::vector< UnsignedWord<Bit> > _numbers;
	_numbers.resize(numbers.size());
	for (int i = 0; i < numbers.size(); ++i) {
		bitify(_numbers[i], numbers[i]);
		std::cout << "bitify 0-" << i << " took " << (clock()-start) << " clocks." << std::endl;
	}

	UnsignedWord<Bit> _min = min(_numbers);
	return _min;
//	Number min;
//	unbitify(min, _min);
//	return min;
}

int main(int, char **) {
	long R=1;
	long p = 5;
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

	std::vector<Number> numbers;
	numbers.resize(2);
	numbers[0] = Converter<Number>::fromInt(1);
	numbers[1] = Converter<Number>::fromInt(2);

//	Number min = doMin3(numbers);
//	std::cout << "Min = " << zp.mod(Converter<Number>::toInt(min)) << std::endl;

	start = clock();
	UnsignedWord<Bit> min = doMin3(numbers);
	std::cout << "Took " << (clock()-start) << " clocks." << std::endl;
//	std::cout << " Min = " << min << std::endl;

}

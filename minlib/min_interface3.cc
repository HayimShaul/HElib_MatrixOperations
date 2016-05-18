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
#include "settings.h"

ZP zp;
clock_t start;

//#define DEBUG

#ifdef DEBUG
#	define EncryptedNumber ZP
#endif

//inline EncryptedNumber power(EncryptedNumber x, int e) {
//	std::cerr << "coputing " << Converter<EncryptedNumber>::toInt(x) << "^" << e << std::endl;
//	if (e == 2) {
//		std::cerr << "before " << Converter<EncryptedNumber>::toInt(x) << "^2" << std::endl;
//		EncryptedNumber res = x*x;
////		std::cerr << "done " << e << std::endl;
//		std::cerr << "after x=" << Converter<EncryptedNumber>::toInt(x) << std::endl;
//		std::cerr << "after res=" << Converter<EncryptedNumber>::toInt(res) << std::endl;
//		std::cerr << "after " << Converter<EncryptedNumber>::toInt(x) << "^2" << " = " << Converter<EncryptedNumber>::toInt(res) << std::endl;
//		return res;
//	}
//	if (e == 1)
//		return x;
//
//	if (e & 1) {
//		EncryptedNumber res = x*power(x, e-1);
//		std::cerr << "done " << e << std::endl;
//		std::cerr << Converter<EncryptedNumber>::toInt(x) << "^" << e << " = " << Converter<EncryptedNumber>::toInt(res) << std::endl;
//		return res;
//	} else {
//		EncryptedNumber m = power(x, e/2);
//		EncryptedNumber res = m*m;
//		std::cerr << "done " << e << std::endl;
//		std::cerr << Converter<EncryptedNumber>::toInt(x) << "^" << e << " = " << Converter<EncryptedNumber>::toInt(res) << std::endl;
//		return res;
//	}
//}

inline EncryptedNumber power(EncryptedNumber x, int e) {
	std::cerr << "computing " << Converter<EncryptedNumber>::toInt(x) << "^" << e << std::endl;
	int msb = 0;
	for (int i = 0; i < 32; ++i)
		if (e & (1 << i))
			msb = i;
	--msb;
	int new_e = 1;
	EncryptedNumber res = x;
	while (msb >= 0) {
		res *= res;
		new_e *= 2;
		if (e & (1 << msb)) {
			res *= x;
			++new_e;
		}
		std::cerr << Converter<EncryptedNumber>::toInt(x) << "^" << new_e << " = " << Converter<EncryptedNumber>::toInt(res) << std::endl;
		--msb;
	}
	return res;
}

EncryptedNumber myMin(const EncryptedNumber &x, const EncryptedNumber &y) {
	EncryptedNumber z = x - y;
	// check if 0 < z < Keys::max_value
	std::cerr << "checking z = " << Converter<EncryptedNumber>::toInt(z) << std::endl;
	std::cerr << "starting copute clock=" << (clock() - start) << std::endl;
	EncryptedNumber semi = Settings<EncryptedNumber>::is_positive(z);
	std::cerr << "starting power clock=" << (clock() - start) << std::endl;
	EncryptedNumber r = power(semi, zp.p() - 1);
	std::cerr << "finished clock=" << (clock() - start) << std::endl;

	return x*r - y*(r-1);
}

EncryptedNumber doMin3(const std::vector<EncryptedNumber> &numbers, int begin = 0, int end = -1) {
	if (end == -1)
		end = numbers.size();

	if (begin + 1 == end)
		return numbers[begin];
	if (begin + 2 == end) {
		EncryptedNumber r = myMin(numbers[begin], numbers[begin+1]);
		return r;
	}

	int mid = (end - begin) / 2;

	EncryptedNumber r1 = doMin3(numbers, begin, mid);
	EncryptedNumber r2 = doMin3(numbers, mid, end);

	EncryptedNumber r = myMin(r1, r2);

	return r;
}


void testMany() {
#ifdef DEBUG
	for (int i = 0 ; i < 100; ++i) {
		std::vector<EncryptedNumber> numbers;
		numbers.resize(2);
		int a = random() % (zp.p() / 2);
		int b = random() % (zp.p() / 2);
		numbers[0] = Converter<EncryptedNumber>::fromInt(a);
		numbers[1] = Converter<EncryptedNumber>::fromInt(b);

		EncryptedNumber _min = doMin3(numbers);
		int realMin = min(a, b);

		if (_min != realMin) {
			std::cerr << "Error min(" << a << ", " << b << ") != " << Converter<EncryptedNumber>::toInt(_min) << std::endl;
			exit(1);
		}
		std::cerr << "OK i="<<i<< std::endl;
	}
#endif
}

int main(int, char **) {
	long R=1;
	long p = 7;
	long r = 1;
	long d = 1;
	long c = 2;
	long k=80;
	long L = 30;
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
	if (p > 200)
		Settings<EncryptedNumber>::max_value(100, zp);
	else
		Settings<EncryptedNumber>::max_value(p/2, zp);

	std::vector<EncryptedNumber> numbers;
	numbers.resize(2);
	numbers[0] = Converter<EncryptedNumber>::fromInt(1);
	numbers[1] = Converter<EncryptedNumber>::fromInt(2);
//	numbers[2] = Converter<EncryptedNumber>::fromInt(1);
//	numbers[3] = Converter<EncryptedNumber>::fromInt(7);

	start = clock();
	EncryptedNumber min = doMin3(numbers);
	std::cerr << "min took " << (clock() - start) << std::endl;

	std::cout << "Min = " << zp.mod(Converter<EncryptedNumber>::toInt(min)) << std::endl;

	testMany();
}

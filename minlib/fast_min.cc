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
#	define EncryptedNumber ZP
#endif

inline EncryptedNumber power(EncryptedNumber x, int e) {
//	std::cerr << "coputing x^" << e << std::endl;
	if (e == 2) {
		EncryptedNumber res = x*x;
		std::cerr << Converter<EncryptedNumber>::toInt(x) << "^2" << " = " << Converter<EncryptedNumber>::toInt(res) << std::endl;
		return res;
	}
	if (e == 1)
		return x;

	if (e & 1) {
		EncryptedNumber res = x*power(x, e-1);
		std::cerr << Converter<EncryptedNumber>::toInt(x) << "^" << e << " = " << Converter<EncryptedNumber>::toInt(res) << std::endl;
		return res;
	} else {
		EncryptedNumber m = power(x, e/2);
		EncryptedNumber res = m*m;
		std::cerr << Converter<EncryptedNumber>::toInt(x) << "^" << e << " = " << Converter<EncryptedNumber>::toInt(res) << std::endl;
		return res;
	}
}


EncryptedNumber mul(const std::vector<EncryptedNumber> &arr, int start = 0, int end = -1) {
	if (end == -1)
		end = arr.size();

	if (start == end - 1)
		return arr[start];

	if (start == end - 2)
		return arr[start] * arr[start + 1];

	int mid = (start + end) / 2;

	return mul(arr, start, mid) * mul(arr, mid, end);
}

inline int log2(long long x) {
	int l = 0;
	while ((x & 1) == 0) {
		++l;
		x >>= 1;
	}
	return l;
}

EncryptedNumber fastMin(EncryptedNumber x, EncryptedNumber y) {
//	std::cerr << "Computing D=x-y" << std::endl;
	// check that D is negative (i.e. D > p/2)
	// negative D => x<y
	EncryptedNumber D = x-y;

//	std::cerr << "Computing E=D*2" << std::endl;
	// check that D2 is odd (odd => x<y, even => y<x)
	EncryptedNumber E = D*2;

//	std::cerr << "Computing L,Linv" << std::endl;
	int L = log2(zp.p() + 1);
//	int Linv = need to be such:  Linv*L = 1 (mod \phi(n))   (\phi(n) = p-1)  which is possible only of L and p-1 are coprime
//	int Linv = findInverse(L, zp.p()-1);
	int Linv = 1;
	std::cerr << "L=" << L << "  Linv=" << Linv << std::endl;

//	std::cerr << "Computing Q\n";
	EncryptedNumber Q1 = power(E, L);
	EncryptedNumber Q = power(Q1, Linv);

	std::vector<EncryptedNumber> Pi;
	Pi.resize(L-1);

	for (int i = 1; i <= L-1; ++i) {
//		std::cerr << "Computing P["<<i<<"]\n";
		Pi[i-1] = E*(1 << i) - Q;
	}

//	std::cerr << "Computing Ptag\n";
	EncryptedNumber Ptag = mul(Pi);

//	std::cerr << "Computing P\n";
	EncryptedNumber P = power(Ptag, zp.p() - 1);

//	std::cerr << "computing min1\n";
	// P=0 => E is even => x-y>0  => x>y

	EncryptedNumber min1 = P*x;
//	std::cerr << "computing min2\n";
	EncryptedNumber min2 = (P-1)*y;

//	std::cerr << "computing min1 + min2\n";
	EncryptedNumber min = min1 - min2;

//	std::cerr << "done\n";

	return min;
}

EncryptedNumber doMin3(const std::vector<EncryptedNumber> &numbers, int begin = 0, int end = -1) {
	if ((zp.p() != (1 << 2)-1) && (zp.p() != (1 << 3)-1) && (zp.p() != (1 << 5)-1) && (zp.p() != (1 << 7)-1) && (zp.p() != (1 << 13)-1) && (zp.p() != (1 << 17)-1) && (zp.p() != (1 << 19)-1) && (zp.p() != (1LL << 31)-1)) {
		std::cerr << "p must be Mersenne prime!" << std::endl;
		exit(1);
	}

	if (end == -1)
		end = numbers.size();

	if (begin + 1 == end)
		return numbers[begin];
	if (begin + 2 == end) {
		EncryptedNumber min = fastMin(numbers[begin], numbers[begin + 1]);
		std::cerr << "min("<< Converter<EncryptedNumber>::toInt(numbers[begin]) << ", " << Converter<EncryptedNumber>::toInt(numbers[begin+1]) << ") = " << Converter<EncryptedNumber>::toInt(min) << std::endl;
		return min;
	}

	int mid = (end - begin) / 2;

	EncryptedNumber min1 = doMin3(numbers, mid, end);
	EncryptedNumber min2 = doMin3(numbers, begin, mid);
	EncryptedNumber min = fastMin(min1, min2);
	std::cerr << "min("<< Converter<EncryptedNumber>::toInt(min1) << ", " << Converter<EncryptedNumber>::toInt(min2) << ") = " << Converter<EncryptedNumber>::toInt(min) << std::endl;
	return min;
}


int many_mins() {
	for (int i = 0; i < 10000; ++i) {
		std::vector<EncryptedNumber> numbers;
		numbers.resize(2);
		numbers[0] = Converter<EncryptedNumber>::fromInt(random());
		numbers[1] = Converter<EncryptedNumber>::fromInt(random());

		EncryptedNumber min = doMin3(numbers);

		if ((min > numbers[0]) || (min > numbers[1]) || ((min != numbers[0]) && (min != numbers[1]))) {
			std::cerr << "Error when computing min(" << Converter<EncryptedNumber>::toInt(numbers[0]) << ", " << Converter<EncryptedNumber>::toInt(numbers[1]) << ") failed" << std::endl;
			exit(1);
		}
	}
}

int main(int, char **) {
	long R=1;
//	long p = (1<<17) - 1;
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
	numbers.resize(2);
	numbers[0] = Converter<EncryptedNumber>::fromInt(2);
	numbers[1] = Converter<EncryptedNumber>::fromInt(1);
//	numbers[2] = Converter<EncryptedNumber>::fromInt(1);
//	numbers[3] = Converter<EncryptedNumber>::fromInt(7);

	std::cerr << "starting min" << std::endl;
	clock_t start = clock();
	EncryptedNumber min = doMin3(numbers);
	std::cout << "Min3 took " << (clock()-start) << " clocks" << std::endl;

	std::cout << "Min = " << zp.mod(Converter<EncryptedNumber>::toInt(min)) << std::endl;

	many_mins();
}

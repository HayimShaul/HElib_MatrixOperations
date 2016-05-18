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

#define TEST_FHE 1

ZP zp;

void test1(const char *name, bool a(), int iter, int seed) {
	if (seed == -1)
		seed = time(NULL);
	srand(seed);
	int i = 0;
	bool ok = true;
	while ((i < iter) && (ok)) {
		if (!a())
			ok = false;
		++i;
	}
	std::cout << "Test " << name << "     ";
	if (ok) {
		std::cout << "OK" << std::endl;
	} else {
		std::cout << "FAILED    (seed=" << seed << "  iter=" << i << ")" << std::endl;
		exit(0);
	}
}

#define testClass(name, func, BIT, iter, seed)                                                           \
	do {                                                                                                 \
		srand(seed);                                                                                     \
		int i = 0;                                                                                       \
		bool ok = true;                                                                                  \
		clock_t start = clock();																		 \
		while ((i < iter) && (ok)) {                                                                     \
			if (!func<BIT>())                                                                            \
				ok = false;                                                                              \
			++i;                                                                                         \
		}                                                                                                \
		std::cout << "Test " << name << "     ";                                                         \
		if (ok) {                                                                                        \
			std::cout << "OK (" << (clock() - start) << " clocks)" << std::endl;                         \
		} else {                                                                                         \
			std::cout << "FAILED    (seed=" << seed << "  iter=" << i << ")" << std::endl;               \
			exit(0);                                                                                     \
		}                                                                                                \
	} while(false)


#if TEST_FHE == 1

#	define testBitClass(name, func, iter, seed)                                                             \
		testClass(name " <PlainBit>", func, PlainBit, iter, seed);										 \
		testClass(name " <EncrypedBit>", func, EncryptedBit, iter, seed);

#	define testNumberClass(name, func, iter, seed)                                                             \
		testClass(name " <int>", func, int, iter, seed);										 \
		testClass(name " <EncryptedNumber>", func, EncryptedNumber, iter, seed);

#else // of TEST_FHE == 1

#	define testBitClass(name, func, iter, seed)                                                             \
		testClass(name " <PlainBit>", func, PlainBit, iter, seed);

#	define testNumberClass(name, func, iter, seed)                                                             \
		testClass(name " <int>", func, int, iter, seed);

#endif // of TEST_FHE == 1







#define testBit(name, func, iter, _seed)                                                                 \
	do {                                                                                                 \
		int seed = _seed;                                                                                \
		if (seed == -1)                                                                                  \
			seed = time(NULL);                                                                           \
		testBitClass(name, func, iter, seed);																 \
	} while(false)

#define testNumber(name, func, iter, _seed)                                                                 \
	do {                                                                                                 \
		int seed = _seed;                                                                                \
		if (seed == -1)                                                                                  \
			seed = time(NULL);                                                                           \
		testNumberClass(name, func, iter, seed);																 \
	} while(false)




template<class Bit>
bool testUnsignedWord() {
	int b1 = UnsignedWord<Bit>::in_range(rand());
	UnsignedWord<Bit> a(b1);

	int b2 = a.get();
	return b2 == b1;
}

template<class Bit>
bool testUnsignedWordPrint() {
	int b1 = UnsignedWord<Bit>::in_range(rand());
	UnsignedWord<PlainBit> a(b1);

	std::bitset<64> out1 = a.to_bitset<64>();
	std::bitset<64> out2(b1);

	return out1 == out2;
}

template<class Bit>
bool testUnsignedWordAddition() {
	int a = UnsignedWord<Bit>::in_range(rand());
	int b = UnsignedWord<Bit>::in_range(rand());

	UnsignedWord<Bit> A(a);
	UnsignedWord<Bit> B(b);

	UnsignedWord<Bit> C = A + B;

	int c = C.get();

	return (c == UnsignedWord<Bit>::in_range(a+b));
}

template<class Bit>
bool testUnsignedWordSub() {
	int a = UnsignedWord<Bit>::in_range(rand());
	int b = UnsignedWord<Bit>::in_range(rand());

	if (a < b) {
		int temp = a;
		a = b;
		b = temp;
	}

	UnsignedWord<Bit> A(a);
	UnsignedWord<Bit> B(b);

	UnsignedWord<Bit> C = A - B;

	int c = C.get();

	return (c == a-b);
}

template<class Bit>
bool testUnsignedWordShr() {
	int a = UnsignedWord<Bit>::in_range(rand());
	int b = UnsignedWord<Bit>::max_bit_num();

	UnsignedWord<Bit> A(a);

	UnsignedWord<Bit> C = A >> b;

	int c = C.get();

	return (c == a >> b);
}

template<class Bit>
bool testUnsignedWordShl() {
	int a = UnsignedWord<Bit>::in_range(rand());
	int b = UnsignedWord<Bit>::max_bit_num();

	UnsignedWord<Bit> A(a);

	UnsignedWord<Bit> C = A << b;

	int c = C.get();

	return (c == UnsignedWord<Bit>::in_range(a << b));
}

template<class Bit>
bool testUnsignedWordMul() {
	int a = UnsignedWord<Bit>::in_range(rand());
	int b = UnsignedWord<Bit>::in_range(rand());

	UnsignedWord<Bit> A(a);
	UnsignedWord<Bit> B(b);

	UnsignedWord<Bit> C = A * B;

	int c = C.get();

	return (c == UnsignedWord<Bit>::in_range(a*b));
}

template<class Bit>
bool testUnsignedWordGreater() {
	int a = UnsignedWord<Bit>::in_range(rand());
	int b = UnsignedWord<Bit>::in_range(rand());

	UnsignedWord<Bit> A(a);
	UnsignedWord<Bit> B(b);

	Bit C = A < B;

	int c = C.get();

	return (c == (a < b));
}

template<class Bit>
bool testUnsignedWordMin() {
	int a = UnsignedWord<Bit>::in_range(rand());
	int b = UnsignedWord<Bit>::in_range(rand());

	UnsignedWord<Bit> A(a);
	UnsignedWord<Bit> B(b);

	UnsignedWord<Bit> C = min(A, B);

	int c = C.get();

	return (c == std::min(a, b));
}

template<class Number>
bool testPolynomial() {
	std::vector<Number> coef;
	coef.resize(5);
	coef[0] = Converter<Number>::fromInt(1);
	coef[1] = Converter<Number>::fromInt(2);
	coef[2] = Converter<Number>::fromInt(3);
	coef[3] = Converter<Number>::fromInt(4);
	coef[4] = Converter<Number>::fromInt(5);

	Polynomial<Number> p(coef);

	Number x;
	x = Converter<Number>::fromInt(2);
	Number result = p.compute(x);

	return zp.mod(Converter<Number>::toInt(result)) == zp.mod(129);
}

template<class Number>
bool testPolynomialTermMul() {
	MulOfPolynomialTerms<Number> p;
	int i;

	// build a polynomial that gets 0 over odd numbers

	for (i = 1; i < 9; i += 2) {
		p *= PolynomialTerm<Number>(typename PolynomialTerm<Number>::X(), -i);
	}

	int ok = true;
	i = 0;
	while ((i < 9) && ok) {
//		fprintf(stdout, "***************** i = %d, compute = %d\n", i, p.compute(i));
		if ((i & 1) == 0)
			ok &= (Converter<Number>::toInt(p.compute(i)) != 0);
		else
			ok &= (Converter<Number>::toInt(p.compute(i)) == 0);
		++i;
	}

	return ok;
}

template<class Number>
bool testSemiCharacteristicPolynomial() {
	MulOfPolynomialTerms<Number> p;

	int  i;

	if (zp.p() != 11) {
		cerr << "Test was written for p=11. Re-run with p=11" << std::endl;
		return false;
	}
	
	std::vector<int> odds;
	for (i = 0; i < zp.p(); ++i) {
		if ((i & 1) == 1)
			odds.push_back(i);
	}

	p.semiCharacteristicTerm(typename PolynomialTerm<Number>::X(), odds, zp);

	int ok = true;
	i = 0;
	while ((i < zp.p()) && ok) {
		fprintf(stderr, "i = %d, compute = %d\n", i, Converter<Number>::toInt(p.compute(i)));
		if ((i & 1) == 0)
			ok &= (Converter<Number>::toInt(p.compute(i)) != 0);
		else
			ok &= (Converter<Number>::toInt(p.compute(i)) == 0);
		++i;
	}

	return ok;
}

template<class Number>
bool testCharacteristicPolynomial() {
	MulOfPolynomialTerms<int> _p;
	MulOfPolynomialTerms<Number> p;

	int  i;

	if (zp.p() != 11) {
		cerr << "Test was written for p=11. Re-run with p=11" << std::endl;
		return false;
	}
	
	
	_p.characteristicTerm(typename PolynomialTerm<int>::X(), 1, zp);

	p.encrypt(_p);

	std::cerr << "poly:" << std::endl << _p << std::endl;
	std::cerr << "encrypted poly:" << std::endl << p << std::endl;

	int ok = true;
	i = 0;
//	while ((i < zp.p()) && ok) {
	while ((i < zp.p()) ) {
		fprintf(stderr, "i = %d, compute = %d\n", i, Converter<Number>::toInt(p.compute(i)));
		if (i == 1)
			ok &= (zp.mod(Converter<Number>::toInt(p.compute(i))) == 1);
		else
			ok &= (zp.mod(Converter<Number>::toInt(p.compute(i))) == 0);
		++i;
	}

	return ok;
}

template<class Number>
bool testCharacteristicMultiPolynomial() {
	MulOfPolynomialTerms<int> _p;
	MulOfPolynomialTerms<Number> p;
	int  i;

	if (zp.p() != 3) {
		cerr << "Test was written for p=3. Re-run with p=3" << std::endl;
		return false;
	}
	
	std::vector<int> ch;
	ch.resize(2);
	for (i = 0; i < 2; ++i)
		ch[i] = i;
	_p.characteristicTerm(ch, zp);

	p.encrypt(_p);

	std::vector<int> x;
	x.resize(ch.size());

	int ok = true;
	for (x[0] = 0; x[0] < zp.p(); ++x[0])
	for (x[1] = 0; x[1] < zp.p(); ++x[1]) {
		std::vector<Number> _x;
		_x.resize(x.size());
		for (i = 0; i < x.size(); ++i)
			_x[i] = Converter<Number>::fromInt(x[i]);
//		if ((x[0] == ch[0]) && (x[1] == ch[1]) && (x2 == ch[2]))
		if ((x[0] == ch[0]) && (x[1] == ch[1]))
			ok &= (zp.mod(Converter<Number>::toInt(p.compute(_x))) == 1);
		else
			ok &= (zp.mod(Converter<Number>::toInt(p.compute(_x))) == 0);
	}
		

	return ok;
}

template<class Number>
bool testMin() {
	SumOfMulOfPolynomialTerms<int> _p;
	SumOfMulOfPolynomialTerms<Number> p;

	if (zp.p() != 3) {
		cerr << "Test was written for p=3. Re-run with p=3" << std::endl;
		return false;
	}
	_p.constructMinPolynomial(3, zp);

	p.encrypt(_p);

	std::vector<int> x;
	x.resize(3);

	for (x[0] = 0; x[0] < zp.p(); ++x[0])
	for (x[1] = 0; x[1] < zp.p(); ++x[1])
	for (x[2] = 0; x[2] < zp.p(); ++x[2]) {
		std::vector<int>::iterator min_it = std::min_element(x.begin(), x.end());
		int min = *min_it;

		std::vector<Number> _x;
		_x.resize(x.size());
		for (int i = 0 ; i < x.size(); ++i)
			_x[i] = Number(x[i]);

		std::cout << std::endl << std::endl << "Computing min polynom with x_0=" << x[0] << " x_1=" << x[1] << " x_2=" << x[2] << std::endl;

		int _min = zp.mod(Converter<Number>::toInt(p.compute(_x)));
		std::cout << "_min=" << _min << "   min=" << min << std::endl;

		if (min != _min)
			return false;
	}
	return true;
}

int main(int argc, char **argv) {
	ArgMapping amap;

//	bool dry=false;
//	amap.arg("dry", dry, "dry=1 for a dry-run");

	long R=1;
	amap.arg("R", R, "number of rounds");

	long p = 2;
	amap.arg("p", p, "plaintext base");

	long r = 1;
	r = 1;
	amap.arg("r", r,  "lifting");

	long d = 1;
	amap.arg("d", d, "degree of the field extension");
	amap.note("d == 0 => factors[0] defines extension");

	long c = 2;
	amap.arg("c", c, "number of columns in the key-switching matrices");
	
	long k=80;
	amap.arg("k", k, "security parameter");

	long L = 10;
	amap.arg("L", L, "# of levels in the modulus chain",  "heuristic");

	long s=0;
	amap.arg("s", s, "minimum number of slots");

//	long repeat=1;
//	amap.arg("repeat", repeat,  "number of times to repeat the test");

	long chosen_m = 0;
	amap.arg("m", chosen_m, "use specified value as modulus", NULL);

//	Vec<long> mvec;
//	amap.arg("mvec", mvec, "use product of the integers as  modulus", NULL);
//	amap.note("e.g., mvec='[5 3 187]' (this overwrite the m argument)");

	Vec<long> gens;
	amap.arg("gens", gens, "use specified vector of generators", NULL);
	amap.note("e.g., gens='[562 1871 751]'");

	Vec<long> ords;
	amap.arg("ords", ords, "use specified vector of orders", NULL);
	amap.note("e.g., ords='[4 2 -4]', negative means 'bad'");

//	long seed=0;
//	amap.arg("seed", seed, "PRG seed");

	amap.parse(argc, argv);

#if TEST_FHE == 1
	Keys::initKeys(s, R, p, r, d, c, k, 64, L, chosen_m, gens, ords);
#endif

	zp.set_p(p);

//	testBit("to/from int", testUnsignedWord, 1, -1);
//	testBit("to/from bitset", testUnsignedWordPrint, 1, -1);
//	testBit("operator +", testUnsignedWordAddition, 1, -1);
//	testBit("operator -", testUnsignedWordSub, 1, -1);
//	testBit("operator <<", testUnsignedWordShl, 1, -1);
//	testBit("operator >>", testUnsignedWordShr, 1, -1);
//	testBit("operator *", testUnsignedWordMul, 1, -1);
//	testBit("operator <", testUnsignedWordGreater, 1, -1);
//	testBit("operator min", testUnsignedWordMin, 1, -1);

//	testNumber("polynomial", testPolynomial, 1, -1);
//	testNumber("MulOfPolynomialTerms", testPolynomialTermMul, 1, -1);
//	testNumber("semi characteristic polynomial", testSemiCharacteristicPolynomial, 1, -1);
//	testNumber("characteristic polynomial", testCharacteristicPolynomial, 1, -1);
//	testNumber("characteristic multi-polynomial", testCharacteristicMultiPolynomial, 1, -1);
	testNumber("min", testMin, 1, -1);
}



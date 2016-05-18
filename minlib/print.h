#ifndef ___PRINT___
#define ___PRINT___

#include <iostream>

template<class Number>
std::ostream &operator<<(std::ostream &out, const PolynomialTerm<Number> &p) {
	out << p.stringify();
	return out;
}

template<class Number>
std::ostream &operator<<(std::ostream &out, const MulOfPolynomialTerms<Number> &p) {
	out << p.stringify();
	return out;
}

template<class Number>
std::ostream &operator<<(std::ostream &out, const SumOfMulOfPolynomialTerms<Number> &p) {
	out << p.stringify();
	return out;
}

#endif

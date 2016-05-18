#ifndef ___SETTINGS___
#define ___SETTINGS___

#include "polynomial.h"

template<class Number>
class Settings {
private:
	static int _max_value;
	static MulOfPolynomialTerms<Number> *_is_positive;
public:
	static void max_value(int, const ZP &);
	static int max_value() { return _max_value; }
	static Number is_positive(const Number &n) { return _is_positive->compute(n); }
};

template<class Number>
MulOfPolynomialTerms<Number> *Settings<Number>::_is_positive = NULL;

template<class Number>
int Settings<Number>::_max_value;

template<class Number>
void Settings<Number>::max_value(int m, const ZP &zp) {
	_max_value = m;

	MulOfPolynomialTerms<int> p;
	std::vector<int> good;
	for (int i = 0; i < m; ++i)
		good.push_back(i);
	p.semiCharacteristicTerm(0, good, zp);
	if (_is_positive != NULL)
		delete _is_positive;
	_is_positive = new MulOfPolynomialTerms<Number>;
	_is_positive->encrypt(p);
}

#endif

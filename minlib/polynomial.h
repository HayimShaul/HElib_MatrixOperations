#ifndef __POLYNOMIAL__
#define __POLYNOMIAL__

#include <algorithm>
#include "assert.h"

#include "zp.h"
#include "converter.h"
#include "cache.h"

//class PlainNumberSuite {
//public:
//	typedef PlainBit Bit;
//	typedef int Number;
//};

//class EncryptedNumberSuite {
//public:
//	typedef EncryptedBit Bit;
//	typedef Ctxt Number;
//};

// Number should be either int or Ctxt
template<class Number>
class Polynomial {
private:
	std::vector<Number> _coef;

	void set(const Number &n, const char *x = NULL);
public:
	Polynomial(const Polynomial<Number> &p) : _coef(p._coef) {}
	Polynomial(const std::vector<Number> &coef) : _coef(coef) {}

	// init one linear factor of a polynomial: Polynomial("x", -4) or Polynomial(4, "-x")
	// Note: x can be either "x", "-x" or "+x"
	Polynomial(const char *x, const Number &n) { set(n, x); }
	Polynomial(const Number &n, const char *x = NULL) { set(n, x); }

//	Polynomial operator*(const Polynomial<Number> &p) const { Polynomial<Number> q(*this); q *= p; return q; }
//	void operator*=(const Polynomial<Number> &p);

	Number compute(const Number &x) const;
};

// Number should be either int or Ctxt
// holds a term of the form (x_i - n), (x_i + n) or (n)
template<class Number>
class PolynomialTerm {
private:
	Number _n;
	int _x_index;

	std::string _signature;
public:
	const Number &n() const { return _n; }
	int x_index() const { return _x_index; }
public:
	struct X {
		int _index;
		X() : _index(0) {}
		X(int i) : _index(i) {}
	};

	PolynomialTerm() : _n(0), _x_index(-1) {}
	PolynomialTerm(const PolynomialTerm<Number> &p) : _n(p._n), _x_index(p._x_index) {}
	PolynomialTerm(const X &x, const Number &n) : _n(n), _x_index(x._index) {}
	~PolynomialTerm() {}

	Number compute(const Number &x) const;
	Number compute(const std::vector<Number> &x) const;

	std::string getSignature() const;

	void encrypt(const PolynomialTerm<int> &p);

	std::string stringify() const;
};

template<class Number>
class MulOfPolynomialTerms {
private:
	std::vector< PolynomialTerm<Number> > _terms;
	Number _factor;
	bool _factor_is_1;

	std::string _signature;
public:
	const std::vector< PolynomialTerm<Number> > &terms() const { return _terms; }
	const Number &factor() const { return _factor; }
	bool factor_is_1() const { return _factor_is_1; }
public:
	MulOfPolynomialTerms() : _factor_is_1(true) {}
	MulOfPolynomialTerms(const MulOfPolynomialTerms<Number> &p) { *this = p; }

	MulOfPolynomialTerms<Number> &operator=(const MulOfPolynomialTerms<Number> &p);

	MulOfPolynomialTerms<Number> operator*(const PolynomialTerm<Number> &p) const { MulOfPolynomialTerms a(*this); a *= p; return a; }
	void operator*=(const PolynomialTerm<Number> &p) { _terms.push_back(p); }

	MulOfPolynomialTerms<Number> operator*(const Number &p) const { MulOfPolynomialTerms a(*this); a *= p; return a; }
	void operator*=(const Number &p);

	MulOfPolynomialTerms<Number> operator*(const MulOfPolynomialTerms<Number> &p) const { MulOfPolynomialTerms a(*this); a *= p; return a; }
	void operator*=(const MulOfPolynomialTerms<Number> &p);

	void semiCharacteristicTerm(const typename PolynomialTerm<Number>::X &x, const std::vector<int> &ch, const ZP &zp);
	void semiCharacteristicTerm(const typename PolynomialTerm<Number>::X &x, int ch, const ZP &zp);
	void characteristicTerm(const typename PolynomialTerm<Number>::X &x, int ch, const ZP &zp);
	void characteristicTerm(const std::vector<int> &ch, const ZP &zp);

	Number compute(const Number &x) const;
	Number compute(const std::vector<Number> &n, Cache<Number> *cache = NULL, int start = 0, int end = -1, const ZP *zp = NULL) const;

	std::string getSignature() const;
	std::string getSignature(int start, int end) const;

	void minPolynomial(int unknowns, int min_index, const ZP &zp);

	void encrypt(const MulOfPolynomialTerms<int> &p);

	std::string stringify() const;
};


template<class Number>
class SumOfMulOfPolynomialTerms {
private:
	std::vector< MulOfPolynomialTerms<Number> > _terms;
	std::string _signature;
public:
	const std::vector< MulOfPolynomialTerms<Number> > &terms() const { return _terms; }
public:
	SumOfMulOfPolynomialTerms() {}
	SumOfMulOfPolynomialTerms(const SumOfMulOfPolynomialTerms<Number> &s) : _terms(s._terms) {}

	SumOfMulOfPolynomialTerms<Number> &operator=(const SumOfMulOfPolynomialTerms<Number> &p) { _terms = p._terms; return *this; }

	SumOfMulOfPolynomialTerms &operator+(const SumOfMulOfPolynomialTerms<Number> &s) const { SumOfMulOfPolynomialTerms<Number> n(*this); n += s; return n; }
	void operator+=(const SumOfMulOfPolynomialTerms<Number> &s) { _terms.insert(_terms.end(), s._terms.begin(), s._terms.end()); }

	void operator+=(const MulOfPolynomialTerms<Number> &s) { _terms.push_back(s); }

	Number compute(const std::vector<Number> &x, Cache<Number> *cache = NULL, int begin = 0, int end = -1, const ZP *zp = NULL) const;

	std::string getSignature(int start, int end) const;

	void constructMinPolynomial(int unknowns, const ZP &zp);
//	void constructMinIndexPolynomial(int unknowns, const ZP &zp);

	void encrypt(const SumOfMulOfPolynomialTerms<int> &p);

	std::string stringify() const;
};


template<class Number>
Number PolynomialTerm<Number>::compute(const Number &x) const {
	return x + _n;
}

template<class Number>
Number PolynomialTerm<Number>::compute(const std::vector<Number> &x) const {
	if (_x_index == -1)
		return _n;
	return x[_x_index] + _n;
}

template<class Number>
std::string PolynomialTerm<Number>::getSignature() const {
	if (_signature != "")
		return _signature;
	return stringify();
}

template<>
inline void PolynomialTerm<ZP>::encrypt(const PolynomialTerm<int> &p) {
	_n = Converter<ZP>::fromInt(p.n());
	_x_index = p.x_index();
	_signature = p.stringify();
}

template<>
inline void PolynomialTerm<EncryptedNumber>::encrypt(const PolynomialTerm<int> &p) {
	_n = Converter<EncryptedNumber>::fromInt(p.n());
	_x_index = p.x_index();
	_signature = p.stringify();
}

template<class Number>
std::string PolynomialTerm<Number>::stringify() const {
	std::ostringstream out;

	if (x_index() == -1)
		out << Converter<Number>::toInt(n());
	else {
		int nn = Converter<Number>::toInt(n());
		out << "(x_" << x_index() << ((nn<0)?" - ":" + ") << abs(nn) << ")";
	}
	return out.str();
}


template<class Number>
void Polynomial<Number>::set(const Number &n, const char *x) {
	_coef.resize(2);
	_coef[0] = n;

	if (x == NULL) {
		_coef.resize(1);
	} else if (x[1] == '-') {
		_coef[1] = -1;
	} else if ((x[1] == '+') || (x[0] == 'x')) {
		_coef[1] = 1;
	} else {
		fprintf(stderr, "unknown term %s\n", x);
		exit(1);
	}
}

template<class Number>
Number Polynomial<Number>::compute(const Number &x) const {
	int i;
	std::vector<Number> powers;

	powers.resize(_coef.size());

//	powers[0] = 1; - we don't actually use it
	powers[1] = x;

	for (i = 2; i < _coef.size(); ++i) {
		if ((i & 1) == 0)
			powers[i] = powers[i/2] * powers[i/2];
		else
			powers[i] = powers[i - 1] * x;
	}

	std::vector<Number> monomials;
	monomials.resize(_coef.size());
	monomials[0] = _coef[0];
	for (i = 1; i < _coef.size(); ++i) {
		monomials[i] = _coef[i] * powers[i];
	}

	return addArray(monomials);
}

//template<class Number>
//void Polynomial<Number>::operator*=(const Polynomial<Number> &p) {
//	std::vector<Number> old_coef = _coef;
//	_coef.resize(old_coef.size() * p._coef.size());
//
//	for (int coef_i = 0; coef_i < _coef.size(); ++coef_i)
//		_coef[coef_i] = 0;
//
//	for (int old_i = 0; old_i < old_coef.size(); ++old_i) {
//		for (int p_i = 0; p_i < p._coef.size(); ++p_i) {
//			_coef[old_i*p_i] += old_coef[old_i] * p._coef[p_i];
//		}
//	}
//}


template<class Number>
MulOfPolynomialTerms<Number> &MulOfPolynomialTerms<Number>::operator=(const MulOfPolynomialTerms<Number> &p) {
	_factor_is_1 = p._factor_is_1;
	_factor = p._factor;
	_terms = p._terms;
	return *this;
}

// This method is really called when building the polynomial when Number==int  the polynomial will be encrypted only later
template<class Number>
void MulOfPolynomialTerms<Number>::operator*=(const Number &p) {
	if (_factor_is_1) {
		if (p != 1) {
			_factor = p;
			_factor_is_1 = false;
		}
	} else {
		_factor *= p;
		if (_factor == 1)
			_factor_is_1 = true;
	}
}

template<class Number>
void MulOfPolynomialTerms<Number>::operator*=(const MulOfPolynomialTerms<Number> &p) {
	if (!p._factor_is_1)
		(*this) *= p._factor;

	_terms.insert(_terms.end(), p._terms.begin(), p._terms.end());
}

template<class Number>
std::string MulOfPolynomialTerms<Number>::getSignature() const {
	if (_signature != "")
		return _signature;
	return stringify();
}

template<class Number>
std::string MulOfPolynomialTerms<Number>::getSignature(int start, int end) const {
	std::string res;
	for (int i = start; i < end; ++i)
		res += _terms[i].getSignature() + std::string(" * ");
	return res;
}

template<class Number>
Number MulOfPolynomialTerms<Number>::compute(const Number &n) const {
	std::vector<Number> _n;
	_n.push_back(n);
	return compute(_n);
}

template<class Number>
Number MulOfPolynomialTerms<Number>::compute(const std::vector<Number> &n, Cache<Number> *cache, int start, int end, const ZP *zp) const {
	Cache<Number> _cache;

	if (cache == NULL)
		cache = &_cache;

	bool mulFactor = false;
	if (end == -1) {
		end = _terms.size();
		mulFactor = !_factor_is_1;
	}

	std::string sig;
	Number res;

	if (mulFactor)
		sig = getSignature();
	else
		sig = getSignature(start, end);

	if (cache->read(sig, res))
		return res;

	if (mulFactor) {
		res = _factor * compute(n, cache, start, end, zp);
	} else if (start == end - 1) {
		res = _terms[start].compute(n);
	} else if (start == end - 2) {
		res = compute(n, cache, start, start+1, zp) * compute(n, cache, start+1, start+2, zp);
	} else {
		int mid = (start + end) / 2;
		res = compute(n, cache, start, mid, zp) * compute(n, cache, mid, end, zp);
	}

	if (zp != NULL) {
//		std::cerr << "mod(" << Converter<Number>::toInt(res);
		res = Converter<Number>::fromInt(zp->mod(Converter<Number>::toInt(res)));
//		std::cerr << ")=" << Converter<Number>::toInt(res) << std::endl;
	}

	cache->write(sig, res);
	return res;
}

template<class Number>
void MulOfPolynomialTerms<Number>::semiCharacteristicTerm(const typename PolynomialTerm<Number>::X &x, const std::vector<int> &ch, const ZP &zp) {
	_terms.resize(0);
	for (int i = 0; i < zp.p(); ++i) {
		if (std::find(ch.begin(), ch.end(), i) != ch.end())
			(*this) *= PolynomialTerm<Number>(x, -i);
	}
}

template<class Number>
void MulOfPolynomialTerms<Number>::semiCharacteristicTerm(const typename PolynomialTerm<Number>::X &x, int ch, const ZP &zp) {
	_terms.resize(0);
	for (int i = 0; i < zp.p(); ++i) {
		if (ch != i)
			(*this) *= PolynomialTerm<Number>(x, -i);
	}
}

template<class Number>
void MulOfPolynomialTerms<Number>::characteristicTerm(const typename PolynomialTerm<Number>::X &x, int ch, const ZP &zp) {
	semiCharacteristicTerm(x, ch, zp);
	Number f = (*this).compute(ch);
	(*this) *= Converter<Number>::fromInt(zp.inv(Converter<Number>::toInt(f)));
}

template<class Number>
void MulOfPolynomialTerms<Number>::characteristicTerm(const std::vector<int> &ch, const ZP &zp) {
	_terms.resize(0);

	for (int i = 0; i < ch.size(); ++i) {
		MulOfPolynomialTerms<Number> xi;
	 	xi.semiCharacteristicTerm(typename PolynomialTerm<Number>::X(i), ch[i], zp);
		(*this) *= xi;
	}
//	std::cerr << "computing " << (*this) << std::endl;
	Number f = (*this).compute(Converter<Number>::fromInt(ch), NULL, 0, -1, &zp);

//	std::cerr << "_factor_is_1=" << _factor_is_1 << " _factor=" << _factor << "  f=" << Converter<Number>::toInt(f) << std::endl;
	(*this) *= Converter<Number>::fromInt(zp.inv(Converter<Number>::toInt(f)));
}

template<>
inline void MulOfPolynomialTerms<int>::encrypt(const MulOfPolynomialTerms<int> &p) { *this = p; }

template<>
inline void MulOfPolynomialTerms<EncryptedNumber>::encrypt(const MulOfPolynomialTerms<int> &p) {
	_terms.resize(p.terms().size());
	for (int i = 0; i < _terms.size(); ++i) {
		_terms[i].encrypt(p.terms()[i]);
	}
	_factor = Converter<EncryptedNumber>::fromInt(p.factor());
	_factor_is_1 = p.factor_is_1();
	_signature = p.stringify();
}

template<>
inline void MulOfPolynomialTerms<ZP>::encrypt(const MulOfPolynomialTerms<int> &p) {
	_terms.resize(p.terms().size());
	for (int i = 0; i < _terms.size(); ++i) {
		_terms[i].encrypt(p.terms()[i]);
	}
	_factor = Converter<ZP>::fromInt(p.factor());
	_factor_is_1 = p.factor_is_1();
	_signature = p.stringify();
}

template<class Number>
std::string MulOfPolynomialTerms<Number>::stringify() const {
	std::ostringstream out;
	out << terms()[0];
	for (int i = 1; i < terms().size(); ++i) {
		out << "*" << terms()[i];
	}
	if (!_factor_is_1)
		out << "*" << Converter<Number>::toInt(_factor);
	return out.str();
}

template<>
inline void SumOfMulOfPolynomialTerms<int>::encrypt(const SumOfMulOfPolynomialTerms<int> &p) { *this = p; }

template<>
inline void SumOfMulOfPolynomialTerms<EncryptedNumber>::encrypt(const SumOfMulOfPolynomialTerms<int> &p) {
	_terms.resize(p.terms().size());
	for (int i = 0; i < _terms.size(); ++i) {
		_terms[i].encrypt(p.terms()[i]);
	}
	_signature = p.stringify();
}

template<>
inline void SumOfMulOfPolynomialTerms<ZP>::encrypt(const SumOfMulOfPolynomialTerms<int> &p) {
	_terms.resize(p.terms().size());
	for (int i = 0; i < _terms.size(); ++i) {
		_terms[i].encrypt(p.terms()[i]);
	}
	_signature = p.stringify();
}


template<class Number>
std::string SumOfMulOfPolynomialTerms<Number>::getSignature(int begin, int end) const {
	std::string res;
	for (int i = begin; i < end; ++i)
		res += _terms[i].getSignature() + std::string(" + ");
	return res;
}

template<class Number>
Number SumOfMulOfPolynomialTerms<Number>::compute(const std::vector<Number> &x, Cache<Number> *cache, int start, int end, const ZP *zp) const {
	Cache<Number> _cache;

	if (cache == NULL)
		cache = &_cache;

//	std::cerr << "cache is of size " << cache->size() << std::endl;

	if (end == -1)
		end = _terms.size();

	std::string sig = getSignature(start, end);
	Number res;

	if (cache->read(sig, res))
		return res;

	if (start == end - 1) {
		res = _terms[start].compute(x, cache, 0, -1, zp);
//		if (Converter<Number>::toInt(res) != 0) {
//			std::cerr << "computed:" << std::endl << _terms[start] << std::endl
//					<< "and got: " << Converter<Number>::toInt(res) << std::endl;
//		}
	} else if (start == end - 2) {
		res = compute(x, cache, start, start+1, zp) + compute(x, cache, start+1, start+2, zp);
	} else {
		int mid = (start + end) / 2;
		res = compute(x, cache, start, mid, zp) + compute(x, cache, mid, end, zp);
	}

	cache->write(sig, res);
	return res;
}

inline bool isMin(std::vector<int> &x, int min) {
	int i;
	for (i = 0; i < min; ++i)
		if (x[i] <= x[min])
			return false;
	for (i = min + 1; i < x.size(); ++i)
		if (x[i] < x[min])
			return false;
	return true;
}

template<class Number>
void SumOfMulOfPolynomialTerms<Number>::constructMinPolynomial(int unknowns, const ZP &zp) {
	std::vector<int> x;
	x.resize(unknowns);

	for (int i = 0; i < unknowns; ++i)
		x[i] = 0;

	while (!zp.isEnd(x)) {
		std::vector<int>::iterator min_it = std::min_element(x.begin(), x.end());
		int min = *min_it;
		
		if (min != 0) {
			MulOfPolynomialTerms<Number> t;
			t.characteristicTerm(x, zp);
			if (min != 1)
				t *= min;
			(*this) += t;
		}
		zp.advance(x);
	}
}

//template<class Number>
//void SumOfMulOfPolynomialTerms<Number>::constructMinIndexPolynomial(int unknowns, const ZP &zp) {
//	std::vector<int> x;
//	x.resize(unknowns);
//
//	for (int i = 0; i < unknowns; ++i)
//		x[i] = 0;
//
//	while (!zp.isEnd(x)) {
//		if (isMin(x, min_index)) {
//			MulOfPolynomialTerms<Number> t;
//			t.characteristicTerm(x, zp);
//			(*this) += t;
//		}
//		zp.advance(x);
//	}
//}

template<class Number>
std::string SumOfMulOfPolynomialTerms<Number>::stringify() const {
	std::ostringstream out;
	out << terms()[0];
	for (int i = 1; i < terms().size(); ++i) {
		out << " + " << terms()[i];
	}
	return out.str();
}

#endif

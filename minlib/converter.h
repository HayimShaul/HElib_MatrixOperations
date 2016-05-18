#ifndef __CONVERTER__
#define __CONVERTER__

#include "encrypted_number.h"
#include "zp.h"

template<class Number>
class Converter {
public:
	static Number fromInt(int);
	static int toInt(const Number &);

	static std::vector<Number> fromInt(const std::vector<int> &v);
	static std::vector<int> toInt(const std::vector<Number> &v);
};


// Int Converter
template<>
inline int Converter<int>::fromInt(int i) { return i; }
template<>
inline int Converter<int>::toInt(const int &i) { return i; }

template<>
inline std::vector<int> Converter<int>::fromInt(const std::vector<int> &v) { return v; }

template<>
inline std::vector<int> Converter<int>::toInt(const std::vector<int> &v) { return v; }

// ZP Converter
template<>
inline ZP Converter<ZP>::fromInt(int i) { return ZP(i); }
template<>
inline int Converter<ZP>::toInt(const ZP &i) { return i.v(); }

template<>
inline std::vector<ZP> Converter<ZP>::fromInt(const std::vector<int> &v) {
	std::vector<ZP> r;
	r.resize(v.size());
	for (int i = 0; i < v.size(); ++i) {
		r[i] = v[i];
	}
	return r;
}

template<>
inline std::vector<int> Converter<ZP>::toInt(const std::vector<ZP> &v) {
	std::vector<int> r;
	r.resize(v.size());
	for (int i = 0; i < v.size(); ++i) {
		r[i] = Converter<ZP>::toInt(v[i]);
	}
	return r;
}


// EncryptedNumber Converter
template<>
inline EncryptedNumber Converter<EncryptedNumber>::fromInt(int i) { return EncryptedNumber(i); }

template<>
inline int Converter<EncryptedNumber>::toInt(const EncryptedNumber &i) { return Keys::decrypt(i.number()); }

template<>
inline std::vector<EncryptedNumber> Converter<EncryptedNumber>::fromInt(const std::vector<int> &v) {
	std::vector<EncryptedNumber> r;
	r.resize(v.size());
	for (int i = 0; i < v.size(); ++i) {
		r[i] = v[i];
	}
	return r;
}

template<>
inline std::vector<int> Converter<EncryptedNumber>::toInt(const std::vector<EncryptedNumber> &v) {
	std::vector<int> r;
	r.resize(v.size());
	for (int i = 0; i < v.size(); ++i) {
		r[i] = Converter<EncryptedNumber>::toInt(v[i]);
	}
	return r;
}

#endif


#ifndef __CTXT_WRAPPER__
#define __CTXT_WRAPPER__

#include <algorithm>
#include <iostream>

#include <NTL/ZZ.h>
#include <FHE.h>
#include <EncryptedArray.h>
#include "keys.h"

class EncryptedNumber {
private:
	Ctxt _ctxt;

	int _depth;
	
public:
	EncryptedNumber() : _ctxt(Keys::publicKey()), _depth(0) {}
	EncryptedNumber(const Ctxt &c) : _ctxt(c), _depth(0) {}
	EncryptedNumber(const EncryptedNumber &c) : _ctxt(c._ctxt), _depth(c._depth) {}
	EncryptedNumber(int c) : _ctxt(Keys::encrypt(c)), _depth(0) {}

	EncryptedNumber &operator=(const EncryptedNumber &c) { _ctxt = c._ctxt; _depth = c._depth; return *this; }
	EncryptedNumber &operator=(const Ctxt &c) { _ctxt = c; _depth = 0; return *this; }
	EncryptedNumber &operator=(int c) { _ctxt = Keys::encrypt(c); _depth = 0; return *this; }

	Ctxt &number() { return _ctxt; }
	const Ctxt &number() const { return _ctxt; }

	void operator*=(const EncryptedNumber &e) {
		clock_t start = clock();
		_ctxt *= e._ctxt;
		std::cerr << "mul(depth " << _depth << ", depth " << e._depth << ") took " << (clock() - start) << std::endl;
		_depth = std::max(_depth, e._depth) + 1;
	}
	EncryptedNumber operator*(const EncryptedNumber &e) const { EncryptedNumber a(*this); a *= e; return a; }

	void operator+=(const EncryptedNumber &e) {
		clock_t start = clock();
		_ctxt += e._ctxt;
		std::cerr << "add(depth " << _depth << ", depth " << e._depth << ") took " << (clock() - start) << std::endl;
		_depth = std::max(_depth, e._depth);
	}
	EncryptedNumber operator+(const EncryptedNumber &e) const { EncryptedNumber a(*this); a += e; return a; }

	void operator-=(const EncryptedNumber &e) {
		clock_t start = clock();
		_ctxt -= e._ctxt;
		std::cerr << "sub(depth " << _depth << ", depth " << e._depth << ") took " << (clock() - start) << std::endl;
		_depth = std::max(_depth, e._depth);
	}
	EncryptedNumber operator-(const EncryptedNumber &e) const { EncryptedNumber a(*this); a -= e; return a; }
};

#endif

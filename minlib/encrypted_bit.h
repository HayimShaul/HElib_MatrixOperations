#ifndef ___ENCRYPTED_BIT___
#define ___ENCRYPTED_BIT___

// Author: Hayim Shaul 2016.

// EncryptedBit.
// implements an encrypted bit. all operations are performed on encrypted bits, and return encrypted bits

#include <NTL/ZZ.h>
#include <FHE.h>
#include <EncryptedArray.h>
#include "keys.h"

#include "encrypted_number.h"

class EncryptedBit {
private:
	Ctxt _bit;

//	void set(int i) { _bit = Keys::encrypt(i); }
	void set(int i) { Keys::encrypt(_bit, i); }
public:

	EncryptedBit() : _bit(Keys::publicKey())  { set(0); }
	EncryptedBit(int b) : _bit(Keys::publicKey()) { set(!!b); }
	EncryptedBit(const EncryptedBit &b) : _bit(Keys::publicKey()) { operator=(b); }

	EncryptedBit &operator=(const EncryptedNumber &b) { _bit = b.number(); return *this; }

	EncryptedBit &operator=(const EncryptedBit &b) { _bit = b._bit; return *this; }
	void operator+=(const EncryptedBit &b) { _bit += b._bit; }
	void operator-=(const EncryptedBit &b) { operator-(b); }
	void operator*=(const EncryptedBit &b) { _bit *= b._bit; }

	EncryptedBit operator+(const EncryptedBit &b) const { EncryptedBit c(*this); c+=b; return c; }
	EncryptedBit operator-(const EncryptedBit &b) const { return operator-(b); }
	EncryptedBit operator*(const EncryptedBit &b) const { EncryptedBit c(*this); c*=b; return c; }

	void neg() { _bit.addConstant(ZZ(1)); }
	EncryptedBit operator!() const { EncryptedBit c(*this); c.neg(); return c; }

	int get() const { return Keys::decrypt(_bit); }
};

#endif

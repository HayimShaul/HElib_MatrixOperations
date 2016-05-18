#include <assert.h>

#include "encrypted_bit.h"

//void EncryptedBit::set(int b) {
//	NewPlaintextArray myPlain(*_ea);
//	_ea->decode(myPlain, ZZX(b));
//	_ea->encrypt(_bit, *_publicKey, myPlain);
//}
//
//
//int EncryptedBit::get() const {
//	NewPlaintextArray myDecrypt(*_ea);
//	_ea->decrypt(_bit, *_secretKey, myDecrypt);
//	ZZX myOutput;
//	_ea->encode(myOutput, myDecrypt);
//	if (myOutput == ZZX(0))
//		return 0;
//	else if (myOutput == ZZX(1))
//  		return 1;
//	else
//		assert(0);
//}

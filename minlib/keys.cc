#include "keys.h"

EncryptedArray *Keys::_ea = NULL;
FHEPubKey *Keys::_publicKey = NULL;
FHESecKey *Keys::_secretKey = NULL;
FHEcontext *Keys::_context = NULL;

void Keys::initKeys(long s, long R, long p, long r, long d, long c, long k, long w, long L, long chosen_m, const Vec<long> &gens, const Vec<long> &ords) {
//  long R=1;
//  long p=2;
//  long r=1;
//  long d=1;
//  long c=2;
//  long k=80;
//  long L=0;
//  long s=0;
//  long repeat=1;
//  long chosen_m=0;
//  Vec<long> gens;
//  Vec<long> ords;
//
//  long seed=0;

	if (L==0) { // determine L based on R,r
		L = 3*R+3;
		if (p>2 || r>1) { // add some more primes for each round
			long addPerRound = 2*ceil(log((double)p)*r*3)/(log(2.0)*FHE_p2Size) +1;
			L += R * addPerRound;
		}
	}
	long m = FindM(k, L, c, p, d, s, chosen_m, true);

	vector<long> gens1, ords1;
	convert(gens1, gens);
	convert(ords1, ords);

	_context = new FHEcontext(m, p, r, gens1, ords1);
	buildModChain(*_context, L, c);

	_secretKey = new FHESecKey(*_context);
	_secretKey->GenSecKey(w); // A Hamming-weight-w secret key

	ZZX G;

	if (d == 0)
		G = _context->alMod.getFactorsOverZZ()[0];
	else
		G = makeIrredPoly(p, d); 

	cerr << "G = " << G << "\n";
	cerr << "generating key-switching matrices... ";
	addSome1DMatrices(*_secretKey); // compute key-switching matrices that we need
	cerr << "done\n";

//	FHESecKey *shadowPubKey = new FHESecKey(*_context);
//	shadowPubKey->GenSecKey(w);
//	addSome1DMatrices(*shadowPubKey); // compute key-switching matrices that we need
//	_publicKey = shadowPubKey;

	_publicKey = _secretKey;


	cerr << "computing masks and tables for rotation...";
	_ea = new EncryptedArray(*_context, G);
	cerr << "done\n";
}

Ctxt Keys::encrypt(int b) {
	Ctxt r(*_publicKey);
	NewPlaintextArray myPlain(*_ea);
	_ea->decode(myPlain, ZZX(b));
	_ea->encrypt(r, *_publicKey, myPlain);
	return r;
}

void Keys::encrypt(Ctxt &r, int b) {
	NewPlaintextArray myPlain(*_ea);
	_ea->decode(myPlain, ZZX(b));
	_ea->encrypt(r, *_publicKey, myPlain);
}


long Keys::decrypt(const Ctxt &b) {
	NewPlaintextArray myDecrypt(*_ea);
	_ea->decrypt(b, *_secretKey, myDecrypt);
	ZZX myOutput;
	_ea->encode(myOutput, myDecrypt);

	if (myOutput == ZZX(0))
		return 0;

	return to_long(myOutput[0]);
}


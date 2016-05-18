#ifndef __KEYS__
#define __KEYS__

#include <NTL/ZZ.h>
#include <FHE.h>
#include <EncryptedArray.h>

class Keys {
private:
	static FHEPubKey *_publicKey;
	static FHESecKey *_secretKey;
	static EncryptedArray *_ea;
	static FHEcontext *_context;
public:
	static void initKeys(long s, long R, long p, long r, long d, long c, long k, long w, long L, long m, const Vec<long> &gens, const Vec<long> &ords);
	static void setKeys(FHEPubKey *pub, FHESecKey *sec, EncryptedArray *ea, FHEcontext *ctx)
		{ _publicKey = pub; _secretKey = sec; _ea = ea; _context = ctx; }

	static FHEPubKey &publicKey() { return *_publicKey; }
	static Ctxt encrypt(int b);
	static void encrypt(Ctxt &c, int b);
	static long decrypt(const Ctxt &b);
};

#endif

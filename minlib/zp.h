#ifndef ___ZP____
#define ___ZP____

class ZP {
private:
	static long long _prev_p;
	long long _p;
	long long _val;

	long long find_inv(long long x, long long p) const;
public:
	ZP() : _p(_prev_p) {}
	ZP(long long v) : _p(_prev_p) { _val = mod(v); }
	ZP(long long v, long long p) : _p(p) { _val = mod(v); }

	static void set_global_p(long long p) { _prev_p = p; }
	void set_p(long long p) { _p = p; _prev_p = _p; }
	long long p() const { return _p; }
	long long inv(long long x) const { return (p() <= 3) ? x : power(x, p() - 2); }
	long long mod(long long x) const { return ((x % p()) + p()) % p(); }
	long long power(long long x, long long e) const;

	void advance(std::vector<int> &x) const;
	bool isEnd(std::vector<int> &x) const { return x[0] == p(); }

	long long v() const { return _val; }
	void operator*=(const ZP &z) { assert(_p == z._p); _val = mod(_val * z._val);  }
	ZP operator*(const ZP &z) const { assert(_p == z._p); return ZP(_val * z._val, _p); }
	ZP operator-(const ZP &z) const { assert(_p == z._p); return ZP(_val - z._val, _p); }
	ZP operator+(const ZP &z) const { assert(_p == z._p); return ZP(_val + z._val, _p); }

	bool operator>(const ZP &z) const { assert(_p == z._p); return _val > z._val; }
	bool operator!=(const ZP &z) const { assert(_p == z._p); return _val != z._val; }
};

inline long long ZP::power(long long x, long long e) const {
//	std::cerr << "coputing x^" << e << std::endl;
	if (e == 2) {
		long long res = mod(x*x);
//		std::cerr << x << "^2" << " = " << res << std::endl;
		return res;
	}
	if (e == 1)
		return x;

	if (e & 1) {
		long long res = mod(x*power(x, e-1));
//		std::cerr << x << "^" << e << " = " << res << std::endl;
		return res;
	} else {
		long long m = power(x, e/2);
		long long res = mod(m*m);
//		std::cerr << x << "^" << e << " = " << res << std::endl;
		return res;
	}
}

//inline int ZP::power(int x, int e) const {
//	if (e == 2)
//		return mod(x*x);
//	if (e == 1)
//		return x;
//	return (e & 1) ? mod(x*power(x, e-1)) : power(power(x, e/2), 2);
//}

inline void ZP::advance(std::vector<int> &x) const {
	int i = x.size() - 1;

	while (i >= 0) {
		if (x[i] < _p - 1) {
			++x[i];
			return;
		}
		x[i] = 0;
		--i;
	}
	x[0] = _p;
}




#endif

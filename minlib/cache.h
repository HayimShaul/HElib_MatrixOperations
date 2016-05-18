#ifndef ___COMP_CACHE___
#define ___COMP_CACHE___

#include <map>
#include "converter.h"

template<class Number>
class Cache {
private:
	std::map<std::string, Number> cache;
public:
	bool read(const std::string &sig, Number &n) {
return false;
//		std::cout << "Trying to read " << sig << std::endl;
		typename std::map<std::string, Number>::iterator it = cache.find(sig);
		if (it == cache.end()) {
//			std::cout << "Cache miss" << std::endl;
//			std::cout << "Not available" << std::endl;
			return false;
		}
//		std::cout << "success! It is " << Converter<Number>::toInt(it->second) << std::endl;
//		std::cout << "Cache hit" << std::endl;
		n = it->second;
		return true;
	}
	void write(const std::string &sig, const Number &n) {
//		std::cout << "Adding cache[" << sig << "] = " << Converter<Number>::toInt(n) << std::endl;
		cache[sig] = n;
	}

	int size() const { return cache.size(); }
};

#endif


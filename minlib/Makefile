CCFLAGS = -g

NTLINCDIR = -I../ntl-9.6.2/include
NTLLIBDIR = -L../ntl-9.6.2/src

FHEINCDIR = -I../HElib-master/src
FHELIBDIR = -L../HElib-master/src

LIBS = $(FHELIBDIR) -lfhe $(NTLLIBDIR) -lntl
INCS = $(NTLINCDIR) $(FHEINCDIR)

FHE_EXT_LIB_HEADERS = plain_bit.h unsigned_word.h polynomial.h

all: min3

#all: test min fast_min min2 min3

test: test.o keys.o polynomial.o
	g++ $(LDFLAGS) -o $@ $^ $(LIBS)

min: min_interface.o keys.o polynomial.o
	g++ $(LDFLAGS) -o $@ $^ $(LIBS)

min2: min_interface2.o keys.o polynomial.o
	g++ $(LDFLAGS) -o $@ $^ $(LIBS)

min3: min_interface3.o keys.o polynomial.o
	g++ $(LDFLAGS) -o $@ $^ $(LIBS)

fast_min: fast_min.o keys.o polynomial.o
	g++ $(LDFLAGS) -o $@ $^ $(LIBS)


%.o: %.cc $(FHE_EXT_LIB_HEADERS)
	g++ $(CCFLAGS) -c  $(INCS) $<

clean:
	rm -f test min min2 min3 fast_min *.o


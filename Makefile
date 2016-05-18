CCFLAGS = -g

NTLINCDIR = -I../ntl-9.6.2/include
NTLLIBDIR = -L../ntl-9.6.2/src

FHEINCDIR = -I../HElib-master/src
FHELIBDIR = -L../HElib-master/src

AVIVINCDIR = -I..
AVIVLIBDIR = -L..

LIBS = $(FHELIBDIR) -lfhe $(NTLLIBDIR) -lntl
INCS = $(NTLINCDIR) $(FHEINCDIR)


k_means: k_means.o
	g++ $(LDFLAGS) -o $@ $^ $(LIBS)

%.o: %.cc $(FHE_EXT_LIB_HEADERS)
	g++ $(CCFLAGS) -c  $(INCS) $<

%.o: %.cpp $(FHE_EXT_LIB_HEADERS)
	g++ $(CCFLAGS) -c  $(INCS) $<


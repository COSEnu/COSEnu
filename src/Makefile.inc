CXX=g++

INCLUDE=
LIBS=

OPT=-O3 -std=c++0x -fopenmp
## -fp-model precise
#-unroll-aggressive -parallel -opt-prefetch

.cpp.o:
	$(CXX) $(OPT) $(INCLUDE) $< -c -o $@

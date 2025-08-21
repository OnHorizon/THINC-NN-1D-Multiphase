CXX=g++
CXXFLAGS=-O3 -Wall -Wextra -MD
SRC := $(wildcard *.cc)
LIBS = -lgsl -lcblas -latlas

OBJ=${SRC:.cc=.o}

%.o: %.cc
	${CXX} ${CXXFLAGS} -o  $@  -c $< ${LIBS}
a.out: ${OBJ}
	${CXX} ${CXXFLAGS} -o $@  ${OBJ} ${LIBS}
clean:
	@rm -f *.o *.d a.out
distclean:
	@rm -f *.o *.d *.csv a.out
scrub:
	@rm -f *.csv

-include $(OBJ:.o=.d)

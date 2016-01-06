# Makefile for TVDN library
# Gilles Puy

# ======== LINKS ========
MLAB 		= /Applications/MATLAB_R2010b.app
MLABINC		= ${MLAB}/extern/include
MLABLIB		= ${MLAB}/extern/lib

# ======== COMPILER ========
CC = cc -Wall -pedantic -no-cpp-precomp
C_OPTIMIZE_SWITCH = -O3 ## For GCC
INC = -I${MLABINC}
LIB = -L${MLABLIB}
CFLAGS = ${C_OPTIMIZE_SWITCH} ${INC} ${LIB}

# ======== MAKE RULES ========
all: mdwt_r.o midwt_r.o

clean:
	/bin/rm mdwt_r.o mdiwt_r.o

%.o : %.c
	${CC} -c ${CFLAGS} $<		


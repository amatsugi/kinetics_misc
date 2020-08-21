CC = gcc
F77 = gfortran
CCFLAG = -O2 -Wall
F77FLAG = -O2 -Wall

RM = rm -f

LLIBS = -lminpack -lopenblas  -lgfortran -lm
#LLIBS = -lminpack -llapack -lblas -lgfortran -lm


all: thfit

thfit: thfit.o thflib.o thflib.h
	$(CC) ${CCFLAG} thfit.o thflib.o $(LLIBS) -o $@

thfit.o : thfit.c thflib.h
	$(CC) ${CCFLAG} -c thfit.c

thflib.o : thflib.c thflib.h
	$(CC) ${CCFLAG} -c thflib.c

clean:
	$(RM) *.o


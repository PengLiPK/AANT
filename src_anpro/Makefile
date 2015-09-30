fc=gfortran


all: cacuns norm rdsac stksac


cacuns : cacuns.o
	$(fc) -o cacuns cacuns.o
cacuns.o : cacuns.f90
	$(fc) -c cacuns.f90
norm : norm.o
	$(fc) -o norm norm.o
norm.o : norm.f90
	$(fc) -c norm.f90
rdsac : rdsac.o
	$(fc) -o rdsac rdsac.o
rdsac.o : rdsac.f90
	$(fc) -c rdsac.f90
stksac : stksac.o
	$(fc) -o stksac stksac.o
stksac.o : stksac.f90
	$(fc) -c stksac.f90

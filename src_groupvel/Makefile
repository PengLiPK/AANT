fc=gfortran


all: ftan


ftan : ftan.o
	$(fc) -o ftan ftan.o
ftan.o : ftan.f90
	$(fc) -c ftan.f90

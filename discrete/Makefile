DEBUG =
# DEBUG = -g -fbacktrace -fbounds-check 
PROFILE =
PROFILE = -pg
PARALLEL =
# OMP
#PARALLEL =  -fopenmp -DOMP
#
# Binary format
#BINFMT = -fconvert=big-endian
BINFMT =
#
#OPTIONS = -O3 -132  $(PARALLEL)   -fpp -cm -unroll -align -align commons  -ip 
OPTIONS = -O3 $(PARALLEL) -fimplicit-none -funroll-loops -mtune=native -march=native -Wall -ftree-vectorizer-verbose=1 -ffast-math  
#
#LOPTIONS=   $(PARALLEL)  -Vaxlib  -static-libcxa
LOPTIONS=   $(PARALLEL) 
#OBJS = IO.o random.o energy.o colisio.o
OBJS=
MODS = commLine.o geometry.o geometryDP.o stepPotentials.o intList.o paramSet.o IO.o random.o energy.o colisio.o ANM.o least_sq_fit.o rmsd.o dims.o
#F77= ifort
F77= /sw/bin/gfortran
FC=$(F77)
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
AS=as
CPP= cpp
CFLAGS = -c $(DEBUG) $(PROFILE) 
LFLAGS = $(PROFILE)

discrete: ${MODS} $(OBJS) main.o
	$(F77) $(LFLAGS) -o discrete $(LOPTIONS) $(OBJS) $(MODS) main.o
commLine.o:
	$(F77) -ffree-form $(OPTIONS) $(PC32) $(CFLAGS) commLine.f
geometry.o:
	$(F77) -ffree-form $(OPTIONS) $(PC32) $(CFLAGS) geometry.f
geometryDP.o: geometry.o
	$(F77) -ffree-form $(OPTIONS) $(PC32) $(CFLAGS) geometryDP.f
stepPotentials.o:
	$(F77) -ffree-form $(OPTIONS) $(PC32) $(CFLAGS) stepPotentials.f
intList.o: stepPotentials.o
	$(F77) -ffree-form $(OPTIONS) $(PC32) $(CFLAGS) intList.f
paramSet.o:
	$(F77) -ffree-form $(OPTIONS) $(PC32) $(CFLAGS) paramSet.f
IO.o: commLine.o geometryDP.o
	$(F77) -ffree-form $(OPTIONS) $(PC32) $(CFLAGS) IO.f
random.o:
	$(F77) -ffree-form $(OPTIONS) $(PC32) $(CFLAGS) random.f
energy.o: stepPotentials.o geometryDP.o intList.o geometry.o random.o
	$(F77) -ffree-form $(OPTIONS) $(PC32) $(CFLAGS) energy.f
colisio.o: geometry.o geometryDP.o intList.o stepPotentials.o 
	$(F77) -ffree-form $(OPTIONS) $(PC32) $(CFLAGS) colisio.f
ANM.o: geometry.o geometryDP.o rmsd.o
	$(F77) -ffree-form $(OPTIONS) $(PC32) $(CFLAGS) ANM.f
least_sq_fit.o:
	$(F77) -ffree-form $(OPTIONS) $(PC32) $(CFLAGS) least_sq_fit.f
rmsd.o: geometryDP.o
	$(F77) -ffree-form $(OPTIONS) $(PC32) $(CFLAGS) rmsd.f
dims.o: geometryDP.o geometry.o random.o ANM.o paramSet.o 
	$(F77) -ffree-form $(OPTIONS) $(PC32) $(CFLAGS) dims.f
main.o: $(MODS)
	$(F77) -ffree-form $(OPTIONS) $(PC32) $(CFLAGS) main.f
clean: 
	rm *.o *.mod discrete

IFLAG = -I
MPILIB = /opt/openmpi/lib


# The compiler

#When using Intel:
FC     = ifort
FFLAGS = -qopenmp -O3

#If you have MPI installed (adjust to your MPI compiler)
#FC = ifort

#FFLAGS = -O2 -cpp -DMPI -I /opt/openmpi/include
#LDFLAGS = -L. -L$(MPILIB)  -lpthread
#LINKFLAGS = 
#F90FLAGS = $(FFLAGS) $(LINKFLAGS) $(INCLUDE) $(LDFLAGS)

#If you dont have MPI installed
#FC = gfortran 
#-fbounds-check -fbacktrace -g (flaggs for debugging) 
#FFLAGS = -fopenmp -O2  -cpp -ffixed-line-length-none -fbounds-check -fbacktrace -g 
#FFLAGS = -fopenmp -O2  -cpp -ffixed-line-length-none 
#LDFLAGS = -L. -lpthread
#LINKFLAGS = 
#F90FLAGS = $(FFLAGS) $(LINKFLAGS) $(INCLUDE) $(LDFLAGS)



# "make" builds all
default: NGFish

Fisher_BS.o: 

OBJFILES = Fisher_Komatsu.o

%.o: %.f90
	$(FC) $(FFLAGS) -c $*.f90 -o $*.o

NGFish:  $(OBJFILES)
	$(FC) -o NGFish $(OBJFILES) $(FFLAGS)

clean: cleanNGFish

cleancflat:
	rm -f *.o *.mod NGFish

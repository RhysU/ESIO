## C++ directives ##

#NOTE: USER MUST LOAD HDF5
# i.e. module load hdf5

# choose compiler
#CC=icc
#CC=icpc #intel
#CC=g++  #gnu
CC=mpicc # mpi compiler

CFLAGS  = #-03 #-g
LDFLAGS =
SOURCES = esio.c
#SOURCES += posix_writer.c

## fortran directives ##

# choose compiler
#FF= gfort  # gnu compiler
#FF=ifort  # intel compiler
FF=mpif90 # mpi compiler

FFLAGS   = #-132
SOURCES += esiof.f90 header.f90 main.f90 

## executable name ##
EXEC=esio_chunk

## building object list ##
OBJECTS = $(addsuffix .o,$(basename $(SOURCES)))

## append include and library paths -- choose what is being used here ##
#INCLUDE= -I$(TACC_HDF5_INC)
#LIBS= -L$(TACC_HDF5_LIB) -lhdf5 -lz
INCLUDE=
LIBS= -lhdf5

## build the executable ##
$(EXEC): $(OBJECTS)
	$(FF) $(LDFLAGS) $(OBJECTS) -o $(EXEC) $(INCLUDE) $(LIBS)

# build cpp object files
%.o: %.cpp
	@echo building $< 
	$(CC) -c $(CFLAGS) $(INCLUDE) $< -o $@

# build c object files
%.o: %.c
	@echo building $<
	$(CC) -c $(CFLAGS) $(INCLUDE) $< -o $@

# build fortran 90 object files
%.o: %.f90
	@echo building $< 
	$(FF) -c $(FFLAGS) $(INCLUDE) $<

# build fortran 77 object files
%.o: %.f
	@echo building $< 
	$(FF) -c $(FFLAGS) $(INCLUDE) $<

# clean directive 'make clean'
clean:
	- /bin/rm $(EXEC) *.o *.mod *~ \#* 
	@echo 'files cleaned'

cleanfield:
	cd field
	- /bin/rm $(EXEC) *.o *.mod *~ \#* *.h5
	@echo 'fields cleaned'

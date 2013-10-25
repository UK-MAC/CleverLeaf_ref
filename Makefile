GIT_VERSION := $(shell git rev-parse --short HEAD)
HOST_NAME := $(shell uname -n)

HDF_DIR=/home/dab/opt/hdf5/1.8.11/intel-13.1.1.163/impi-4.1.0.24
HDF_INC=-I$(HDF_DIR)/include
HDF_LIB=-L$(HDF_DIR)/lib -lhdf5

BOOST_DIR=/home/dab/opt/boost/1.52.0/intel-13.1.1.163
BOOST_INC=-I$(BOOST_DIR)/include

SAMRAI_DIR=/home/dab/opt/SAMRAI/3.6.3/intel-13.1.1.163/impi-4.1.0.24/opt
SAMRAI_INC=-I$(SAMRAI_DIR)/include
SAMRAI_LDIR=$(SAMRAI_DIR)/lib
SAMRAI_LIB=$(SAMRAI_LDIR)/libSAMRAI_appu.a $(SAMRAI_LDIR)/libSAMRAI_algs.a $(SAMRAI_LDIR)/libSAMRAI_solv.a $(SAMRAI_LDIR)/libSAMRAI_geom.a $(SAMRAI_LDIR)/libSAMRAI_mesh.a $(SAMRAI_LDIR)/libSAMRAI_math.a $(SAMRAI_LDIR)/libSAMRAI_pdat.a $(SAMRAI_LDIR)/libSAMRAI_xfer.a $(SAMRAI_LDIR)/libSAMRAI_hier.a $(SAMRAI_LDIR)/libSAMRAI_tbox.a

LAPACK_DIR=/home/dab/opt/lapack/3.4.2/intel-13.1.1.163
BLAS_DIR=/home/dab/opt/lapack/3.4.2/intel-13.1.1.163
MATH_INC=-I$(LAPACK_DIR)/include
MATH_LIB=-L$(LAPACK_DIR)/lib -llapack -lblas

CXX=mpiicpc
F90=mpiifort

CPPFLAGS=-O3 -fp-model source -fp-model strict -prec-div -prec-sqrt -lz $(BOOST_INC) $(HDF_INC) $(SAMRAI_INC) $(MATH_INC) -DVERSION=\"$(GIT_VERSION)\" -DHOST_NAME=\"$(HOST_NAME)\"
FFLAGS=-O3 -warn all -fp-model strict -fp-model source -prec-div -prec-sqrt -module obj/
LDFLAGS=-lz $(SAMRAI_LIB) $(HDF_LIB) $(MATH_LIB) -lstdc++ -nofor_main

CPP_FILES := $(wildcard src/*.C)
F90_FILES := $(wildcard src/fortran/*.f90)
OBJ_FILES := $(addprefix obj/,$(notdir $(F90_FILES:.f90=.o) $(CPP_FILES:.C=.o)))

ref: cleverleaf

openmp: CPPFLAGS+=-openmp
openmp: LDFLAGS+=-openmp
openmp: FFLAGS+=-openmp
openmp: ref

cleverleaf: obj $(OBJ_FILES)
	$(F90) $^ $(LDFLAGS) -o $@

obj/%.o: src/%.C
	$(CXX) $(CPPFLAGS) -c -o $@ $<

obj/%.o: src/fortran/%.f90
	$(F90) $(FFLAGS) -c -o $@ $<

obj: 
	mkdir -p obj

clean:
	rm -rf obj/*.o obj/*.mod cleverleaf doc/dox

doc:
	( cat doc/Doxyfile ; echo "PROJECT_NUMBER=$(GIT_VERSION)") | doxygen -
	rsync doc/dox/html/* /shared/general/docs/Cleverleaf/ -r

test: cleverleaf
	mpirun -n 1 ./cleverleaf test/cleverleaf_test.in

.PHONY: clean, doc

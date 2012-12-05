HDF_DIR=/home/dab/opt/hdf5/1.8.7/intel-12/ompi-1.4.3
HDF_INC=-I$(HDF_DIR)/include
HDF_LIB=-L$(HDF_DIR)/lib -lhdf5

SAMRAI_DIR=/home/dab/opt/SAMRAI/3.4.1/intel-12/ompi-1.4.3
SAMRAI_INC=-I$(SAMRAI_DIR)/include
SAMRAI_LDIR=-I$(SAMRAI_DIR)/lib
SAMRAI_LIB=$(SAMRAI_LDIR)/libSAMRAI_appu.a $(SAMRAI_LDIR)/libSAMRAI_algs.a $(SAMRAI_LDIR)/libSAMRAI_solv.a $(SAMRAI_LDIR)/libSAMRAI_geom.a $(SAMRAI_LDIR)/libSAMRAI_mesh.a $(SAMRAI_LDIR)/libSAMRAI_math.a $(SAMRAI_LDIR)/libSAMRAI_pdat.a $(SAMRAI_LDIR)/libSAMRAI_xfer.a $(SAMRAI_LDIR)/libSAMRAI_hier.a $(SAMRAI_LDIR)/libSAMRAI_tbox.a

LAPACK_DIR=/home/dab/opt/lapack/3.3.1/intel-12/serial
BLAS_DIR=/home/dab/opt/lapack/3.3.1/intel-12/serial
MATH_INC=-I$(LAPACK_DIR)/include
MATH_LIB=-L$(LAPACK_DIR)/lib -llapack -lblas

CXX=mpiCC

CPPFLAGS=-g -O3 -ffloat-store -lz $(HDF_INC) $(SAMRAI_INC) $(MATH_INC)
LDFLAGS=-g -lz $(SAMRAI_LIB) $(HDF_LIB) $(MATH_LIB) -lstdc++ -lgfortran

CPP_FILES := $(wildcard src/*.C)
OBJ_FILES := $(addprefix obj/,$(notdir $(CPP_FILES:.C=.o)))

all: obj cleverleaf

cleverleaf: $(OBJ_FILES)
	$(CXX) $^ $(LDFLAGS) -o $@

obj/%.o: src/%.C
	$(CXX) $(CPPFLAGS) -c -o $@ $<

obj: 
	mkdir -p obj

clean:
	rm -f obj/*.o cleverleaf

docs:
	doxygen doc/Doxyfile
	rsync doc/dox/html/* /shared/general/docs/Cleverleaf/ -r

cleandoc:
	rm -rf doc/dox

.PHONY: clean, doc

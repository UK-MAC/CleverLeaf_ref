HYPRE_DIR=/home/dab/opt/hypre/2.7.0b/gnu/ompi-1.4.3
HYPRE_INC=-I$(HYPRE_DIR)/include
HYPRE_LIB=-L$(HYPRE_DIR)/lib -lHYPRE

HDF_DIR=/home/dab/opt/hdf5/1.8.7/gnu/ompi-1.4.3
HDF_INC=-I$(HDF_DIR)/include
HDF_LIB=-L$(HDF_DIR)/lib -lhdf5

SAMRAI_DIR=/home/dab/opt/SAMRAI/3.3.3/gnu/ompi-1.4.3
SAMRAI_INC=-I$(SAMRAI_DIR)/include
SAMRAI_LIB=$(SAMRAI_DIR)/lib/libSAMRAI_appu.a $(SAMRAI_DIR)/lib/libSAMRAI_algs.a $(SAMRAI_DIR)/lib/libSAMRAI_solv.a $(SAMRAI_DIR)/lib/libSAMRAI_geom.a $(SAMRAI_DIR)/lib/libSAMRAI_mesh.a $(SAMRAI_DIR)/lib/libSAMRAI_math.a $(SAMRAI_DIR)/lib/libSAMRAI_pdat.a $(SAMRAI_DIR)/lib/libSAMRAI_xfer.a $(SAMRAI_DIR)/lib/libSAMRAI_hier.a $(SAMRAI_DIR)/lib/libSAMRAI_tbox.a

LAPACK_DIR=/home/dab/opt/lapack/3.3.1/gnu/serial
BLAS_DIR=/home/dab/opt/lapack/3.3.1/gnu/serial
MATH_INC=-I$(LAPACK_DIR)/include
MATH_LIB=-L$(LAPACK_DIR)/lib -llapack -lblas

CXX=mpiCC

CPPFLAGS=-g -DDEBUG -lz $(HYPRE_INC) $(HDF_INC) $(SAMRAI_INC) $(MATH_INC)
LDFLAGS=-g -lz $(SAMRAI_LIB) $(HYPRE_LIB) $(HDF_LIB) $(MATH_LIB)

all: cleverleaf

cleverleaf: LagrangianEulerianPatchStrategy.o LagrangianEulerianIntegrator.o Cleverleaf.o main.o 
	$(CXX) $^ $(LDFLAGS) -o $@

%.o: src/%.C
	$(CXX) $(CPPFLAGS) -c $<

clean:
	rm -f *.o cleverleaf

.PHONY: clean

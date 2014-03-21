COMPILER=INTEL

HDF_DIR=/home/dab/opt/hdf5/1.8.11/intel-13.1.1.163/impi-4.1.0.24
HDF_INC=-I$(HDF_DIR)/include
HDF_LIB=-L$(HDF_DIR)/lib -lhdf5

BOOST_DIR=/home/dab/opt/boost/1.52.0/intel-13.1.1.163
BOOST_INC=-I$(BOOST_DIR)/include

SAMRAI_DIR=/home/dab/opt/SAMRAI/3.7.3/intel-13.1.1.163/impi-4.1.0.24/opt
SAMRAI_INC=-I$(SAMRAI_DIR)/include
SAMRAI_LDIR=$(SAMRAI_DIR)/lib
SAMRAI_LIB=$(SAMRAI_LDIR)/libSAMRAI_appu.a $(SAMRAI_LDIR)/libSAMRAI_algs.a $(SAMRAI_LDIR)/libSAMRAI_solv.a $(SAMRAI_LDIR)/libSAMRAI_geom.a $(SAMRAI_LDIR)/libSAMRAI_mesh.a $(SAMRAI_LDIR)/libSAMRAI_math.a $(SAMRAI_LDIR)/libSAMRAI_pdat.a $(SAMRAI_LDIR)/libSAMRAI_xfer.a $(SAMRAI_LDIR)/libSAMRAI_hier.a $(SAMRAI_LDIR)/libSAMRAI_tbox.a

LAPACK_DIR=/home/dab/opt/lapack/3.4.2/intel-13.1.1.163
BLAS_DIR=/home/dab/opt/lapack/3.4.2/intel-13.1.1.163
MATH_INC=-I$(LAPACK_DIR)/include
MATH_LIB=-L$(LAPACK_DIR)/lib -llapack -lblas

CXX=mpiicpc
F90=mpiifort

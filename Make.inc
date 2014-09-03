COMPILER=INTEL

HDF_DIR=$(HOME)/opt/libs/hdf5/1.8.12
HDF_INC=-I$(HDF_DIR)/include
HDF_LIB=-L$(HDF_DIR)/lib -lhdf5

BOOST_DIR=$(HOME)/opt/libs/boost/1.52.0
BOOST_INC=-I$(BOOST_DIR)/include

SAMRAI_DIR=$(HOME)/opt/libs/SAMRAI/3.7.3
SAMRAI_INC=-I$(SAMRAI_DIR)/include
SAMRAI_LDIR=$(SAMRAI_DIR)/lib
SAMRAI_LIB=$(SAMRAI_LDIR)/libSAMRAI_appu.a $(SAMRAI_LDIR)/libSAMRAI_algs.a $(SAMRAI_LDIR)/libSAMRAI_solv.a $(SAMRAI_LDIR)/libSAMRAI_geom.a $(SAMRAI_LDIR)/libSAMRAI_mesh.a $(SAMRAI_LDIR)/libSAMRAI_math.a $(SAMRAI_LDIR)/libSAMRAI_pdat.a $(SAMRAI_LDIR)/libSAMRAI_xfer.a $(SAMRAI_LDIR)/libSAMRAI_hier.a $(SAMRAI_LDIR)/libSAMRAI_tbox.a


CXX=mpiicpc-13.1.163
F90=mpiifort-13.1.163

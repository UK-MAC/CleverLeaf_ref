# 
# Copyright 2013 David Beckingsale.
# 
# This file is part of CleverLeaf.
# 
# CleverLeaf is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
# 
# CleverLeaf is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
# more details.
# 
# You should have received a copy of the GNU General Public License along with
# CleverLeaf. If not, see http://www.gnu.org/licenses/.
#  
include Make.inc

GIT_VERSION := $(shell git rev-parse --short HEAD)
HOST_NAME := $(shell uname -n)

CXXFLAGS_INTEL=-O3 -ipo -fp-model source -fp-model strict -prec-div -prec-sqrt
FFLAGS_INTEL=-O3 -ipo -fpe0 -warn all -fp-model strict -fp-model source -prec-div -prec-sqrt -module obj/
# CXXFLAGS_INTEL=-O0 -g -debug all -fp-model source -fp-model strict -prec-div -prec-sqrt
# FFLAGS_INTEL=-O0 -g -fpe0 -warn all -debug all -ftrapuv -check uninit -fp-model strict -fp-model source -prec-div -prec-sqrt -module obj/
LDFLAGS_INTEL=-ipo -nofor_main
#LDFLAGS_INTEL=-nofor_main -O0 -g -debug all
OMP_INTEL=-openmp

CXXFLAGS_GNU=-O3 -march=native -funroll-loops -ffloat-store
FFLAGS_GNU=-O3 -march=native -funroll-loops -ffloat-store
LDFLAGS_GNU=
OMP_GNU=-fopenmp

CXXFLAGS=$(CXXFLAGS_$(COMPILER)) \
				 -lz $(BOOST_INC) $(HDF_INC) $(SAMRAI_INC) $(MATH_INC) \
				 -DVERSION=\"$(GIT_VERSION)\" -DHOST_NAME=\"$(HOST_NAME)\"
FFLAGS=$(FFLAGS_$(COMPILER)) 
LDFLAGS=$(LDFLAGS_$(COMPILER)) -lz $(SAMRAI_LIB) $(HDF_LIB) $(MATH_LIB) -lstdc++

CPP_FILES := $(wildcard src/*.C)
F90_FILES := $(wildcard src/fortran/*.f90)
OBJ_FILES := $(addprefix obj/,$(notdir $(F90_FILES:.f90=.o) $(CPP_FILES:.C=.o)))

ref: obj cleverleaf

openmp: CXXFLAGS+=$(OMP_$(COMPILER))
openmp: LDFLAGS+=$(OMP_$(COMPILER))
openmp: FFLAGS+=$(OMP_$(COMPILER))
openmp: ref

cleverleaf: $(OBJ_FILES)
	$(F90) $^ $(LDFLAGS) -o $@

obj/%.o: src/%.C
	$(CXX) $(CXXFLAGS) -c -o $@ $<

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

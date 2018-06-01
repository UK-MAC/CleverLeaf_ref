//////////////////////////////////////////////////////////////////////////////
// Crown Copyright 2014 AWE, Copyright 2014 David Beckingsale.
//
// This file is part of CleverLeaf.
//
// CleverLeaf is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.
//
// CleverLeaf is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with
// CleverLeaf. If not, see http://www.gnu.org/licenses/.
//////////////////////////////////////////////////////////////////////////////
#ifndef CLEVERLEAF_PDAT_CLEVERSIDEVARIABLE_C_
#define CLEVERLEAF_PDAT_CLEVERSIDEVARIABLE_C_

#include "CleverSideVariable.h"

#include "CleverSideDataFactory.h"


namespace clever {
namespace pdat {

template<typename TYPE>
CleverSideVariable<TYPE>::CleverSideVariable(
    const SAMRAI::tbox::Dimension& dim,
    const std::string& name,
    int depth):
  SAMRAI::hier::Variable(name,
      std::make_shared<CleverSideDataFactory<TYPE> >(depth,
        SAMRAI::hier::IntVector::getZero(dim)))
{
}

template<typename TYPE>
CleverSideVariable<TYPE>::~CleverSideVariable(){}

template<typename TYPE>
bool CleverSideVariable<TYPE>::fineBoundaryRepresentsVariable() const
{
  return true;
}

template<typename TYPE>
bool CleverSideVariable<TYPE>::dataLivesOnPatchBorder() const
{
  return true;
}

template<typename TYPE>
int CleverSideVariable<TYPE>::getDepth() const
{
  std::shared_ptr<CleverSideDataFactory<TYPE> > clever_side_data_factory(
      SHARED_PTR_CAST(CleverSideDataFactory<TYPE> ,
        getPatchDataFactory()));

  return clever_side_data_factory->getDepth();
}

}
}

#endif // CLEVERLEAF_PDAT_CLEVERSIDEVARIABLE_C_

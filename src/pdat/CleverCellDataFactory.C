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
#ifndef CLEVERLEAF_PDAT_CLEVERCELLDATAFACTORY_C_
#define CLEVERLEAF_PDAT_CLEVERCELLDATAFACTORY_C_

#include "CleverCellDataFactory.h"

#include "SAMRAI/tbox/MemoryUtilities.h"
#include "SAMRAI/pdat/CellGeometry.h"

#include "CleverCellData.h"

#include "macros.h"

namespace clever {
namespace pdat {

template<typename TYPE>
CleverCellDataFactory<TYPE>::CleverCellDataFactory(
    int depth,
    const SAMRAI::hier::IntVector& ghosts):
  SAMRAI::hier::PatchDataFactory(ghosts),
  d_depth(depth)
{
}

template<typename TYPE>
CleverCellDataFactory<TYPE>::~CleverCellDataFactory(){}

template<typename TYPE>
std::shared_ptr<SAMRAI::hier::PatchDataFactory> 
CleverCellDataFactory<TYPE>::cloneFactory(
    const SAMRAI::hier::IntVector& ghosts)
{
  return std::make_shared<CleverCellDataFactory>(d_depth, ghosts);
}

template<typename TYPE>
std::shared_ptr<SAMRAI::hier::PatchData> CleverCellDataFactory<TYPE>::allocate(
    const SAMRAI::hier::Patch& patch) const
{
  return std::make_shared<CleverCellData<TYPE> >(
      patch.getBox(),
      d_depth,
      d_ghosts);
}

template<typename TYPE>
std::shared_ptr<SAMRAI::hier::BoxGeometry>
CleverCellDataFactory<TYPE>::getBoxGeometry(const SAMRAI::hier::Box& box) const
{
  return std::make_shared<SAMRAI::pdat::CellGeometry>(box, d_ghosts);
}

template<typename TYPE>
size_t CleverCellDataFactory<TYPE>::getSizeOfMemory(const SAMRAI::hier::Box& box) const
{
  const size_t obj = 
    SAMRAI::tbox::MemoryUtilities::align(sizeof(CleverCellData<TYPE>));
  const size_t data = 
    CleverCellData<TYPE>::getSizeOfData(box, d_depth, d_ghosts);

  return obj + data;
}

template<typename TYPE>
bool CleverCellDataFactory<TYPE>::fineBoundaryRepresentsVariable() const
{
  return true;
}

template<typename TYPE>
bool CleverCellDataFactory<TYPE>::dataLivesOnPatchBorder() const
{
  return false;
}

template<typename TYPE>
bool CleverCellDataFactory<TYPE>::validCopyTo(
    const std::shared_ptr<SAMRAI::hier::PatchDataFactory>& dst_pdf) const
{
  bool valid_copy = false;

  std::shared_ptr<CleverCellDataFactory> cell_data_factory(
      SHARED_PTR_CAST(CleverCellDataFactory,
        dst_pdf));

  if(cell_data_factory) {
    valid_copy = true;
  }

  return valid_copy;
}

template<typename TYPE>
int CleverCellDataFactory<TYPE>::getDepth()
{
  return d_depth;
}

}
}

#endif // CLEVERLEAF_PDAT_CLEVERCELLDATAFACTORY_C_

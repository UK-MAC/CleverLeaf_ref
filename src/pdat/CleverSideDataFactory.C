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
#ifndef CLEVERLEAF_PDAT_CLEVERSIDEDATAFACTORY_C_
#define CLEVERLEAF_PDAT_CLEVERSIDEDATAFACTORY_C_

#include "CleverSideDataFactory.h"

#include "SAMRAI/tbox/MemoryUtilities.h"
#include "SAMRAI/pdat/SideGeometry.h"

#include "macros.h"

#include "CleverSideData.h"

namespace clever {
namespace pdat {

template<typename TYPE>
CleverSideDataFactory<TYPE>::CleverSideDataFactory(
    int depth,
    const SAMRAI::hier::IntVector& ghosts):
  SAMRAI::hier::PatchDataFactory(ghosts),
  d_depth(depth)
{
}

template<typename TYPE>
CleverSideDataFactory<TYPE>::~CleverSideDataFactory(){}

template<typename TYPE>
std::shared_ptr<SAMRAI::hier::PatchDataFactory> 
CleverSideDataFactory<TYPE>::cloneFactory(
    const SAMRAI::hier::IntVector& ghosts)
{
  return std::make_shared<CleverSideDataFactory>(d_depth, ghosts);
}

template<typename TYPE>
std::shared_ptr<SAMRAI::hier::PatchData> CleverSideDataFactory<TYPE>::allocate(
    const SAMRAI::hier::Patch& patch) const
{
  return std::make_shared<CleverSideData<TYPE> >(
      patch.getBox(),
      d_depth,
      d_ghosts);
}

template<typename TYPE>
std::shared_ptr<SAMRAI::hier::BoxGeometry>
CleverSideDataFactory<TYPE>::getBoxGeometry(const SAMRAI::hier::Box& box) const
{
  return std::make_shared<SAMRAI::pdat::SideGeometry>(box, d_ghosts);
}

template<typename TYPE>
size_t CleverSideDataFactory<TYPE>::getSizeOfMemory(const SAMRAI::hier::Box& box) const
{
  const size_t obj = 
    SAMRAI::tbox::MemoryUtilities::align(sizeof(CleverSideData<TYPE>));
  const size_t data = 
    CleverSideData<TYPE>::getSizeOfData(box, d_depth, d_ghosts);

  return obj + data;
}

template<typename TYPE>
bool CleverSideDataFactory<TYPE>::fineBoundaryRepresentsVariable() const
{
  return true;
}

template<typename TYPE>
bool CleverSideDataFactory<TYPE>::dataLivesOnPatchBorder() const
{
  return true;
}

template<typename TYPE>
bool CleverSideDataFactory<TYPE>::validCopyTo(
    const std::shared_ptr<SAMRAI::hier::PatchDataFactory>& dst_pdf) const
{
  bool valid_copy = false;

  std::shared_ptr<CleverSideDataFactory> side_data_factory(
      SHARED_PTR_CAST(CleverSideDataFactory,
        dst_pdf));

  if(side_data_factory) {
    valid_copy = true;
  }

  return valid_copy;
}

template<typename TYPE>
int CleverSideDataFactory<TYPE>::getDepth()
{
  return d_depth;
}

}
}

#endif // CLEVERLEAF_PDAT_CLEVERSIDEDATAFACTORY_C_

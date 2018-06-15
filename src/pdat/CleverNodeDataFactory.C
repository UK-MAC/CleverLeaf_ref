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
#ifndef CLEVERLEAF_PDAT_CLEVERNODEDATAFACTORY_C_
#define CLEVERLEAF_PDAT_CLEVERNODEDATAFACTORY_C_

#include "CleverNodeDataFactory.h"

#include "SAMRAI/tbox/MemoryUtilities.h"
#include "SAMRAI/pdat/NodeGeometry.h"

#include "CleverNodeData.h"

#include "macros.h"

namespace clever {
namespace pdat {

template<typename TYPE>
CleverNodeDataFactory<TYPE>::CleverNodeDataFactory(
    int depth,
    const SAMRAI::hier::IntVector& ghosts):
  SAMRAI::hier::PatchDataFactory(ghosts),
  d_depth(depth)
{
}

template<typename TYPE>
CleverNodeDataFactory<TYPE>::~CleverNodeDataFactory(){}

template<typename TYPE>
std::shared_ptr<SAMRAI::hier::PatchDataFactory> 
CleverNodeDataFactory<TYPE>::cloneFactory(
    const SAMRAI::hier::IntVector& ghosts)
{
  return std::make_shared<CleverNodeDataFactory>(d_depth, ghosts);
}

template<typename TYPE>
std::shared_ptr<SAMRAI::hier::PatchData> CleverNodeDataFactory<TYPE>::allocate(
    const SAMRAI::hier::Patch& patch) const
{
  return std::make_shared<CleverNodeData<TYPE> >(
      patch.getBox(),
      d_depth,
      d_ghosts);
}

template<typename TYPE>
std::shared_ptr<SAMRAI::hier::BoxGeometry>
CleverNodeDataFactory<TYPE>::getBoxGeometry(const SAMRAI::hier::Box& box) const
{
  return std::make_shared<SAMRAI::pdat::NodeGeometry>(box, d_ghosts);
}

template<typename TYPE>
size_t CleverNodeDataFactory<TYPE>::getSizeOfMemory(const SAMRAI::hier::Box& box) const
{
  const size_t obj = 
    SAMRAI::tbox::MemoryUtilities::align(sizeof(CleverNodeData<TYPE>));
  const size_t data = 
    CleverNodeData<TYPE>::getSizeOfData(box, d_depth, d_ghosts);

  return obj + data;
}

template<typename TYPE>
bool CleverNodeDataFactory<TYPE>::fineBoundaryRepresentsVariable() const
{
  return true;
}

template<typename TYPE>
bool CleverNodeDataFactory<TYPE>::dataLivesOnPatchBorder() const
{
  return true;
}

template<typename TYPE>
bool CleverNodeDataFactory<TYPE>::validCopyTo(
    const std::shared_ptr<SAMRAI::hier::PatchDataFactory>& dst_pdf) const
{
  bool valid_copy = false;

  std::shared_ptr<CleverNodeDataFactory> node_data_factory(
      SHARED_PTR_CAST(CleverNodeDataFactory,
        dst_pdf));

  if(node_data_factory) {
    valid_copy = true;
  }

  return valid_copy;
}

template<typename TYPE>
int CleverNodeDataFactory<TYPE>::getDepth()
{
  return d_depth;
}

}
}

#endif // CLEVERLEAF_PDAT_CLEVERNODEDATAFACTORY_C_

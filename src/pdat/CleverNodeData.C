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
#ifndef CLEVERLEAF_PDAT_CLEVERNODEDATA_C_
#define CLEVERLEAF_PDAT_CLEVERNODEDATA_C_

#include "CleverNodeData.h"

#include "SAMRAI/pdat/NodeOverlap.h"
#include "SAMRAI/pdat/NodeGeometry.h"

namespace clever {
namespace pdat {

template<typename TYPE>
size_t CleverNodeData<TYPE>::getSizeOfData(
    const SAMRAI::hier::Box& box,
    int depth,
    const SAMRAI::hier::IntVector& ghosts)
{
  const SAMRAI::hier::Box node_box = 
    SAMRAI::pdat::NodeGeometry::toNodeBox(box);
  const SAMRAI::hier::Box ghost_box = SAMRAI::hier::Box::grow(node_box, ghosts);

  return AlignedArrayData<TYPE, 64>::getSizeOfData(ghost_box, depth);
}

template<typename TYPE>
CleverNodeData<TYPE>::CleverNodeData(
    const SAMRAI::hier::Box& domain,
    int depth,
    const SAMRAI::hier::IntVector& ghosts):
  SAMRAI::hier::PatchData(domain,ghosts),
  d_depth(depth)
{
  const SAMRAI::hier::Box node_box = 
    SAMRAI::pdat::NodeGeometry::toNodeBox(getGhostBox());

  d_array_data.reset(new AlignedArrayData<TYPE, 64>(node_box, depth));
}

template<typename TYPE>
CleverNodeData<TYPE>::~CleverNodeData(){}

template<typename TYPE>
TYPE* CleverNodeData<TYPE>::getPointer(int depth)
{
  return d_array_data->getPointer(depth);
}

template<typename TYPE>
void CleverNodeData<TYPE>::fill(const TYPE& value)
{
  d_array_data->fill(value);
}

template<typename TYPE>
void CleverNodeData<TYPE>::fillAll(const TYPE& value)
{
  d_array_data->fillAll(value);
}

template<typename TYPE>
void CleverNodeData<TYPE>::copy(const SAMRAI::hier::PatchData& src)
{
  const CleverNodeData* node_data_source = 
    dynamic_cast<const CleverNodeData<TYPE> *>(&src);

  if (node_data_source == 0) {
    src.copy2(*this);
  } else {
    const SAMRAI::hier::Box box = 
      d_array_data->getBox() * node_data_source->d_array_data->getBox();
    if (!box.empty()) {
      d_array_data->copy(*(node_data_source->d_array_data), box);
    }
  }
}

template<typename TYPE>
void CleverNodeData<TYPE>::copy2(SAMRAI::hier::PatchData& dst) const
{
  CleverNodeData* node_data_destination = 
   static_cast<CleverNodeData<TYPE> *>(&dst); 

  const SAMRAI::hier::Box box = 
    d_array_data->getBox()*node_data_destination->d_array_data->getBox();

  if (!box.empty()) {
    node_data_destination->d_array_data->copy(*d_array_data, box);
  }
}

template<typename TYPE>
void CleverNodeData<TYPE>::copy(
    const SAMRAI::hier::PatchData& src,
    const SAMRAI::hier::BoxOverlap& overlap)
{
  const CleverNodeData* node_data_source = 
    dynamic_cast<const CleverNodeData<TYPE> *>(&src);

  const SAMRAI::pdat::NodeOverlap* node_overlap = 
    dynamic_cast<const SAMRAI::pdat::NodeOverlap *>(&overlap);

  if ((node_data_source == 0) || (node_overlap == 0)) {
    src.copy2(*this, overlap);
  } else {
    if (node_overlap->getTransformation().getRotation() ==
        SAMRAI::hier::Transformation::NO_ROTATE) {
      d_array_data->copy(*(node_data_source->d_array_data),
          node_overlap->getDestinationBoxContainer(),
          node_overlap->getTransformation());
    } else {
      std::cerr << "CleverNodeData: rotated boxes not supported in copy" 
        << std::endl;
      std::abort();
    }
  }
}

template<typename TYPE>
void CleverNodeData<TYPE>::copy2(
    SAMRAI::hier::PatchData& dst, 
    const SAMRAI::hier::BoxOverlap& overlap) const
{
  CleverNodeData* node_data_destination =
    static_cast<CleverNodeData<TYPE> *>(&dst);

  const SAMRAI::pdat::NodeOverlap* node_overlap =
    static_cast<const SAMRAI::pdat::NodeOverlap *>(&overlap);

  if (node_overlap->getTransformation().getRotation() ==
      SAMRAI::hier::Transformation::NO_ROTATE) {
    node_data_destination->d_array_data->copy(*d_array_data,
        node_overlap->getDestinationBoxContainer(),
        node_overlap->getTransformation());
  } else {
      std::cerr << "CleverNodeData: rotated boxes not supported in copy2"
        << std::endl;
      std::abort();
  }
}

template<typename TYPE>
bool CleverNodeData<TYPE>::canEstimateStreamSizeFromBox() const
{
    return true;
}

template<typename TYPE>
size_t CleverNodeData<TYPE>::getDataStreamSize(
    const SAMRAI::hier::BoxOverlap& overlap) const
{
  const SAMRAI::pdat::NodeOverlap* node_overlap = 
    static_cast<const SAMRAI::pdat::NodeOverlap *>(&overlap);

  return d_array_data->getDataStreamSize(
      node_overlap->getDestinationBoxContainer());
}

template<typename TYPE>
void CleverNodeData<TYPE>::packStream(
    SAMRAI::tbox::MessageStream& stream,
        const SAMRAI::hier::BoxOverlap& overlap) const
{
  const SAMRAI::pdat::NodeOverlap* node_overlap = 
    static_cast<const SAMRAI::pdat::NodeOverlap *>(&overlap);

  if(node_overlap->getTransformation().getRotation()
      == SAMRAI::hier::Transformation::NO_ROTATE)
  {
    d_array_data->packStream(
        stream,
        node_overlap->getDestinationBoxContainer(),
        node_overlap->getTransformation());
  } else {
    std::cerr << "CleverNodeData: rotated boxes not supported in packStream!" 
      << std::endl;
    std::abort();
  }
}

template<typename TYPE>
void CleverNodeData<TYPE>::unpackStream(SAMRAI::tbox::MessageStream& stream,
    const SAMRAI::hier::BoxOverlap& overlap)
{
  const SAMRAI::pdat::NodeOverlap* node_overlap = 
    static_cast<const SAMRAI::pdat::NodeOverlap *>(&overlap);

  d_array_data->unpackStream(
      stream,
      node_overlap->getDestinationBoxContainer(),
      node_overlap->getSourceOffset());
}

template<typename TYPE>
int CleverNodeData<TYPE>::getDepth() const
{
  return d_depth;
}

}
}

#endif // CLEVERLEAF_PDAT_CLEVERNODEDATA_C_

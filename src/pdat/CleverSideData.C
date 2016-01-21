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
#ifndef CLEVERLEAF_PDAT_CLEVERSIDEDEDATA_C_
#define CLEVERLEAF_PDAT_CLEVERSIDEDEDATA_C_

#include "CleverSideData.h"

#include "SAMRAI/pdat/SideOverlap.h"
#include "SAMRAI/pdat/SideGeometry.h"

namespace clever {
namespace pdat {

template<typename TYPE>
size_t CleverSideData<TYPE>::getSizeOfData(
    const SAMRAI::hier::Box& box,
    int depth,
    const SAMRAI::hier::IntVector& ghosts)
{
  size_t size = 0;
  const SAMRAI::hier::Box ghost_box = SAMRAI::hier::Box::grow(box, ghosts);

  for (int d = 0; d < box.getDim().getValue(); ++d) {
      const SAMRAI::hier::Box side_box = 
        SAMRAI::pdat::SideGeometry::toSideBox(ghost_box, d);
      size += AlignedArrayData<TYPE, 64>::getSizeOfData(side_box, depth);
  }

  return size;
}

template<typename TYPE>
CleverSideData<TYPE>::CleverSideData(
    const SAMRAI::hier::Box& domain,
    int depth,
    const SAMRAI::hier::IntVector& ghosts):
  SAMRAI::hier::PatchData(domain,ghosts),
  d_depth(depth)
{
  for (int dimension = 0; dimension < getDim().getValue(); dimension++) {
    const SAMRAI::hier::Box side_box = 
      SAMRAI::pdat::SideGeometry::toSideBox(getGhostBox(), dimension);

    d_array_data[dimension].reset(new AlignedArrayData<TYPE, 64>(side_box, depth));
  }
}

template<typename TYPE>
CleverSideData<TYPE>::~CleverSideData(){}

template<typename TYPE>
TYPE* CleverSideData<TYPE>::getPointer(int axis, int depth)
{
  return d_array_data[axis]->getPointer(depth);
}

template<typename TYPE>
void CleverSideData<TYPE>::fill(const TYPE& value)
{
  for (int dimension = 0; dimension < getDim().getValue(); dimension++) {
    d_array_data[dimension]->fill(value);
  }
}

template<typename TYPE>
void CleverSideData<TYPE>::fillAll(const TYPE& value)
{
  for (int dimension = 0; dimension < getDim().getValue(); dimension++) {
    d_array_data[dimension]->fillAll(value);
  }
}

template<typename TYPE>
void CleverSideData<TYPE>::copy(const SAMRAI::hier::PatchData& src)
{
  const CleverSideData* side_data_source = 
    dynamic_cast<const CleverSideData<TYPE> *>(&src);

  if (side_data_source == 0) {
    src.copy2(*this);
  } else {
    for (int dimension = 0; dimension < getDim().getValue(); dimension++) {
      const SAMRAI::hier::Box box = 
        d_array_data[dimension]->getBox() * side_data_source->d_array_data[dimension]->getBox();
      if (!box.empty()) {
        d_array_data[dimension]->copy(*(side_data_source->d_array_data[dimension]), box);
      }
    }
  }
}

template<typename TYPE>
void CleverSideData<TYPE>::copy2(SAMRAI::hier::PatchData& dst) const
{
  CleverSideData* side_data_destination = 
   static_cast<CleverSideData<TYPE> *>(&dst); 

  for (int dimension = 0; dimension < getDim().getValue(); dimension++) {
    const SAMRAI::hier::Box box = 
      d_array_data[dimension]->getBox()*side_data_destination->d_array_data[dimension]->getBox();

    if (!box.empty()) {
      side_data_destination->d_array_data[dimension]->copy(*d_array_data[dimension], box);
    }
  }
}

template<typename TYPE>
void CleverSideData<TYPE>::copy(
    const SAMRAI::hier::PatchData& src,
    const SAMRAI::hier::BoxOverlap& overlap)
{
  const CleverSideData* side_data_source = 
    dynamic_cast<const CleverSideData<TYPE> *>(&src);

  const SAMRAI::pdat::SideOverlap* side_overlap = 
    dynamic_cast<const SAMRAI::pdat::SideOverlap *>(&overlap);

  if ((side_data_source == 0) || (side_overlap == 0)) {
    src.copy2(*this, overlap);
  } else {
    if (side_overlap->getTransformation().getRotation() ==
        SAMRAI::hier::Transformation::NO_ROTATE) {
      for (int dimension = 0; dimension < getDim().getValue(); dimension++) {
        d_array_data[dimension]->copy(*(side_data_source->d_array_data[dimension]),
            side_overlap->getDestinationBoxContainer(dimension),
            side_overlap->getTransformation());
      }
    } else {
      std::cerr << "CleverSideData: rotated boxes not supported in copy" 
        << std::endl;
      std::abort();
    }
  }
}

template<typename TYPE>
void CleverSideData<TYPE>::copy2(
    SAMRAI::hier::PatchData& dst, 
    const SAMRAI::hier::BoxOverlap& overlap) const
{
  CleverSideData* side_data_destination =
    static_cast<CleverSideData<TYPE> *>(&dst);

  const SAMRAI::pdat::SideOverlap* side_overlap =
    static_cast<const SAMRAI::pdat::SideOverlap *>(&overlap);

  if (side_overlap->getTransformation().getRotation() ==
      SAMRAI::hier::Transformation::NO_ROTATE) {
    for (int dimension = 0; dimension < getDim().getValue(); dimension++) {

      const SAMRAI::hier::Transformation& transformation = side_overlap->getTransformation();

      side_data_destination->d_array_data[dimension]->copy(*d_array_data[dimension],
          side_overlap->getDestinationBoxContainer(dimension),
          transformation);
    }
  } else {
    std::cerr << "CleverSideData: rotated boxes not supported in copy2"
      << std::endl;
    std::abort();
  }
}

template<typename TYPE>
bool CleverSideData<TYPE>::canEstimateStreamSizeFromBox() const
{
    return true;
}

template<typename TYPE>
size_t CleverSideData<TYPE>::getDataStreamSize(
    const SAMRAI::hier::BoxOverlap& overlap) const
{
  const SAMRAI::pdat::SideOverlap* side_overlap = 
    static_cast<const SAMRAI::pdat::SideOverlap *>(&overlap);

  int size = 0;
  for (int dimension = 0; dimension < getDim().getValue(); dimension++) {
    size += d_array_data[dimension]->getDataStreamSize(
        side_overlap->getDestinationBoxContainer(dimension));
  }

  return size;
}

template<typename TYPE>
void CleverSideData<TYPE>::packStream(
    SAMRAI::tbox::MessageStream& stream,
        const SAMRAI::hier::BoxOverlap& overlap) const
{
  const SAMRAI::pdat::SideOverlap* side_overlap = 
    static_cast<const SAMRAI::pdat::SideOverlap *>(&overlap);

  if(side_overlap->getTransformation().getRotation()
      == SAMRAI::hier::Transformation::NO_ROTATE)
  {
    const SAMRAI::hier::Transformation& transformation = 
      side_overlap->getTransformation();

    for (int dimension = 0; dimension < getDim().getValue(); dimension++) {
      const SAMRAI::hier::BoxContainer& boxes = 
        side_overlap->getDestinationBoxContainer(dimension);

      d_array_data[dimension]->packStream(stream, boxes, transformation);
    }
  } else {
    std::cerr << "CleverSideData: rotated boxes not supported in packStream!" 
      << std::endl;
    std::abort();
  }
}

template<typename TYPE>
void CleverSideData<TYPE>::unpackStream(SAMRAI::tbox::MessageStream& stream,
    const SAMRAI::hier::BoxOverlap& overlap)
{
  const SAMRAI::pdat::SideOverlap* side_overlap = 
    static_cast<const SAMRAI::pdat::SideOverlap *>(&overlap);

  const SAMRAI::hier::IntVector& offset = side_overlap->getSourceOffset();

  for (int dimension = 0; dimension < getDim().getValue(); dimension++) {
    const SAMRAI::hier::BoxContainer& boxes = 
      side_overlap->getDestinationBoxContainer(dimension);

    if (!boxes.isEmpty()) {
      d_array_data[dimension]->unpackStream(stream, boxes, offset);
    }
  }
}

template<typename TYPE>
int CleverSideData<TYPE>::getDepth() const
{
  return d_depth;
}

}
}

#endif // CLEVERLEAF_PDAT_CLEVERSIDEDEDATA_C_

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
#ifndef CLEVERLEAF_PDAT_ALIGNEDARRAYDATA_C_
#define CLEVERLEAF_PDAT_ALIGNEDARRAYDATA_C_

#include "AlignedArrayData.h"
#include "AlignedArrayDataOperationUtilities.h"

#include "SAMRAI/tbox/MemoryUtilities.h"

#include "stdlib.h"

namespace clever {
namespace pdat {

template<typename TYPE, int ALIGNMENT>
size_t AlignedArrayData<TYPE, ALIGNMENT>::getSizeOfData(
    const SAMRAI::hier::Box& box,
    int depth)
{
  return SAMRAI::tbox::MemoryUtilities::align(
      box.size() * depth * sizeof(TYPE));
}

template<typename TYPE, int ALIGNMENT>
AlignedArrayData<TYPE, ALIGNMENT>::AlignedArrayData(
    const SAMRAI::hier::Box& box, int depth):
  d_depth(depth),
  d_offset(box.size()),
  d_box(box)
{
  posix_memalign((void **)&d_array, ALIGNMENT, d_depth*d_offset*sizeof(TYPE));
}

template<typename TYPE, int ALIGNMENT>
AlignedArrayData<TYPE, ALIGNMENT>::~AlignedArrayData()
{
  free(d_array);
}

template<typename TYPE, int ALIGNMENT>
void AlignedArrayData<TYPE, ALIGNMENT>::packStream(
    SAMRAI::tbox::MessageStream& stream,
    const SAMRAI::hier::BoxContainer& destination_boxes,
    const SAMRAI::hier::Transformation& transformation) const
{
   const int size = d_depth * destination_boxes.getTotalSizeOfBoxes();
   TYPE* buffer;
   posix_memalign((void **)&buffer, ALIGNMENT, size*sizeof(TYPE));

   int pointer = 0;
   for (SAMRAI::hier::BoxContainer::const_iterator box = 
       destination_boxes.begin();
       box != destination_boxes.end();
       ++box)
   {
     SAMRAI::hier::Box pack_box(*box);
     transformation.inverseTransform(pack_box);

     packBuffer(&buffer[pointer], pack_box);

     pointer += d_depth*pack_box.size();
   }

   stream.pack(&buffer[0], size);
   free(buffer);
}

template<typename TYPE, int ALIGNMENT>
void AlignedArrayData<TYPE, ALIGNMENT>::packBuffer(
   TYPE* buffer,
   const SAMRAI::hier::Box& box) const
{
   AlignedArrayDataOperationUtilities<TYPE, ALIGNMENT>::
     doArrayDataBufferOperationOnBox<false>(*this,
         buffer,
         box);
}


template<typename TYPE, int ALIGNMENT>
void AlignedArrayData<TYPE, ALIGNMENT>::unpackStream(
    SAMRAI::tbox::MessageStream& stream,
    const SAMRAI::hier::BoxContainer& destination_boxes,
    const SAMRAI::hier::IntVector& source_offset) const
{
  const int size = d_depth*destination_boxes.getTotalSizeOfBoxes();
  TYPE* buffer;
  posix_memalign((void **)&buffer, ALIGNMENT, size*sizeof(TYPE));

  stream.unpack(&buffer[0], size);

  int pointer = 0;
  for (SAMRAI::hier::BoxContainer::const_iterator box = 
      destination_boxes.begin();
      box != destination_boxes.end();
      ++box)
  {
    unpackBuffer(&buffer[pointer], *box);
    pointer += d_depth*box->size();
  }

  free(buffer);
}

template<typename TYPE, int ALIGNMENT>
void AlignedArrayData<TYPE, ALIGNMENT>::unpackBuffer(
   TYPE* buffer,
   const SAMRAI::hier::Box& box) const
{
   AlignedArrayDataOperationUtilities<TYPE, ALIGNMENT>::
     doArrayDataBufferOperationOnBox<true>(*this,
         buffer,
         box);
}

template<typename TYPE, int ALIGNMENT>
void AlignedArrayData<TYPE, ALIGNMENT>::copy(
    const AlignedArrayData<TYPE, ALIGNMENT>& src,
    const SAMRAI::hier::Box& box)
{
  if ((d_depth == src.d_depth) &&
      (d_box.isSpatiallyEqual(src.d_box)) &&
      (box.isSpatiallyEqual(d_box))) {

    const int n = d_offset * d_depth;
    std::copy(&src.d_array[0], &src.d_array[n], d_array);
  } else {
    const SAMRAI::hier::Box copybox = box * d_box * src.d_box;

    if (!copybox.empty()) {
      const int dst_start_depth = 0;
      const int src_start_depth = 0;
      const int num_depth = (d_depth < src.d_depth ? d_depth : src.d_depth);
      const SAMRAI::hier::IntVector src_shift(box.getDim(), 0);

      AlignedArrayDataOperationUtilities<TYPE, ALIGNMENT>::
        doArrayDataOperationOnBox(*this,
            src,
            copybox,
            src_shift,
            dst_start_depth,
            src_start_depth,
            num_depth);
    }

  }
}

template<typename TYPE, int ALIGNMENT>
void AlignedArrayData<TYPE, ALIGNMENT>::copy(
    const AlignedArrayData<TYPE, ALIGNMENT>& src,
    const SAMRAI::hier::BoxContainer& boxes,
    const SAMRAI::hier::Transformation& transformation)
{
   for (SAMRAI::hier::BoxContainer::const_iterator b = boxes.begin();
       b != boxes.end(); ++b) {

     SAMRAI::hier::Box box(*b);

     if (transformation.getRotation() == SAMRAI::hier::Transformation::NO_ROTATE
         && transformation.getOffset() == SAMRAI::hier::IntVector::getZero(box.getDim())
         && transformation.getBeginBlock() == transformation.getEndBlock()) {

       copy(src, box);

     } else {
       SAMRAI::hier::Box transformed_src(src.d_box);
       transformation.transform(transformed_src);
       const SAMRAI::hier::Box copybox(
           box * d_box * transformed_src);

       if (!copybox.empty()) {
         const int dst_start_depth = 0;
         const int src_start_depth = 0;
         const int num_depth = (d_depth < src.d_depth ? d_depth : src.d_depth);

         AlignedArrayDataOperationUtilities<TYPE, ALIGNMENT>::
           doArrayDataOperationOnBox(*this,
               src,
               copybox,
               transformation.getOffset(),
               dst_start_depth,
               src_start_depth,
               num_depth);
       }
     }
   }
}

template<typename TYPE, int ALIGNMENT>
int AlignedArrayData<TYPE, ALIGNMENT>::getDataStreamSize(
    const SAMRAI::hier::BoxContainer& boxes) const
{
  const int number_of_elements = boxes.getTotalSizeOfBoxes();

  return SAMRAI::tbox::MessageStream::getSizeof<TYPE>(
      d_depth*number_of_elements);
}

template<typename TYPE, int ALIGNMENT>
TYPE* AlignedArrayData<TYPE, ALIGNMENT>::getPointer(int depth) const
{
  return &d_array[depth*d_offset];
}

template<typename TYPE, int ALIGNMENT>
const SAMRAI::hier::Box& AlignedArrayData<TYPE, ALIGNMENT>::getBox() const
{
  return d_box;
}

template<typename TYPE, int ALIGNMENT>
const SAMRAI::tbox::Dimension& AlignedArrayData<TYPE, ALIGNMENT>::getDim() const
{
  return d_box.getDim();
}

template<typename TYPE, int ALIGNMENT>
int AlignedArrayData<TYPE, ALIGNMENT>::getOffset() const
{
  return d_offset;
}

template<typename TYPE, int ALIGNMENT>
int AlignedArrayData<TYPE, ALIGNMENT>::getDepth() const
{
  return d_depth;
}

template<typename TYPE, int ALIGNMENT>
void AlignedArrayData<TYPE, ALIGNMENT>::fill(const TYPE& value) 
{
  fillAll(value);
}

template<typename TYPE, int ALIGNMENT>
void AlignedArrayData<TYPE, ALIGNMENT>::fillAll(const TYPE& value)
{
   if (!d_box.empty()) {
      TYPE* ptr = &d_array[0];
      const int n = d_depth * d_offset;
      for (int i = 0; i < n; ++i) {
         ptr[i] = value;
      }
   }
}

}
}

#endif // CLEVERLEAF_PDAT_ALIGNEDARRAYDATA_C_

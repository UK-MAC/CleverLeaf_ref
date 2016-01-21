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
#ifndef CLEVERLEAF_PDAT_ALIGNEDARRAYDATA_H_
#define CLEVERLEAF_PDAT_ALIGNEDARRAYDATA_H_

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/BoxContainer.h"
#include "SAMRAI/hier/Transformation.h"
#include "SAMRAI/tbox/MessageStream.h"

namespace clever {
namespace pdat {

template<typename TYPE, int ALIGNMENT>
class AlignedArrayData
{
  public:
    static size_t getSizeOfData(
        const SAMRAI::hier::Box& box,
        int depth);

    AlignedArrayData(const SAMRAI::hier::Box& box, int depth);

    ~AlignedArrayData();

    void packStream(
        SAMRAI::tbox::MessageStream& stream,
        const SAMRAI::hier::BoxContainer& destination_boxes,
        const SAMRAI::hier::Transformation& transformation) const;

    void unpackStream(
        SAMRAI::tbox::MessageStream& stream,
        const SAMRAI::hier::BoxContainer& destination_boxes,
        const SAMRAI::hier::IntVector& source_offset) const;

    void copy(
        const AlignedArrayData<TYPE, ALIGNMENT>& source,
        const SAMRAI::hier::Box& box);

    void copy(
        const AlignedArrayData<TYPE, ALIGNMENT>& src,
        const SAMRAI::hier::BoxContainer& boxes,
        const SAMRAI::hier::Transformation& transformation);

    size_t getDataStreamSize(
        const SAMRAI::hier::BoxContainer& boxes) const;

    TYPE* getPointer(int depth = 0) const;

    const SAMRAI::hier::Box& getBox() const;
    const SAMRAI::tbox::Dimension& getDim() const;
    int getDepth() const;
    int getOffset() const;

    void fill(const TYPE& value);
    void fillAll(const TYPE& value);
  private:
    TYPE* d_array;

    SAMRAI::hier::Box d_box;
    int d_depth;
    int d_offset;

    void packBuffer(
        TYPE* buffer,
        const SAMRAI::hier::Box& box) const;

    void unpackBuffer(
        TYPE* buffer,
        const SAMRAI::hier::Box& box) const;

    void copy(
        const AlignedArrayData<TYPE, ALIGNMENT>& src,
        const SAMRAI::hier::Box& box,
        const SAMRAI::hier::Transformation& transformation);
};

}
}

#include "AlignedArrayData.C"

#endif // CLEVERLEAF_PDAT_ALIGNEDARRAYDATA_H_

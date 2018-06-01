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
#ifndef CLEVERLEAF_PDAT_CLEVERNODEDATA_H_
#define CLEVERLEAF_PDAT_CLEVERNODEDATA_H_

#include "AlignedArrayData.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/hier/BoxOverlap.h"
#include "SAMRAI/tbox/MessageStream.h"

namespace clever {
namespace pdat {

template<typename TYPE>
class CleverNodeData : public SAMRAI::hier::PatchData
{
  public:
    static size_t getSizeOfData(
        const SAMRAI::hier::Box& box,
        int depth,
        const SAMRAI::hier::IntVector& ghosts);

    CleverNodeData(
        const SAMRAI::hier::Box& box,
        int depth,
        const SAMRAI::hier::IntVector& ghosts);

    virtual ~CleverNodeData();

    void copy(const SAMRAI::hier::PatchData& src);

    void copy2(SAMRAI::hier::PatchData& dst) const;

    void copy(const SAMRAI::hier::PatchData& src,
        const SAMRAI::hier::BoxOverlap& overlap);

    void copy2(SAMRAI::hier::PatchData& dst,
        const SAMRAI::hier::BoxOverlap& overlap) const;

    bool canEstimateStreamSizeFromBox() const;

    size_t getDataStreamSize(const SAMRAI::hier::BoxOverlap& overlap) const;

    void packStream(SAMRAI::tbox::MessageStream& stream,
        const SAMRAI::hier::BoxOverlap& overlap) const;

    void unpackStream(SAMRAI::tbox::MessageStream& stream,
        const SAMRAI::hier::BoxOverlap& overlap);

    int getDepth() const;

    TYPE* getPointer(int depth=0);

    void fill(const TYPE& value);

    void fillAll(const TYPE& value);
  private:
    std::shared_ptr<AlignedArrayData<TYPE, 64> > d_array_data;
    int d_depth;
};

}
}

#include "CleverNodeData.C"

#endif // CLEVERLEAF_PDAT_CLEVERNODEDATA_H_

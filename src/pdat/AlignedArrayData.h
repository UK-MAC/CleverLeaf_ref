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

    int getDataStreamSize(
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

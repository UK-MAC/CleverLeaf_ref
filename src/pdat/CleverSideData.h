#ifndef CLEVERLEAF_PDAT_CLEVERSIDEDEDATA_H_
#define CLEVERLEAF_PDAT_CLEVERSIDEDEDATA_H_

#include "AlignedArrayData.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/hier/BoxOverlap.h"
#include "SAMRAI/tbox/MessageStream.h"

namespace clever {
namespace pdat {

template<typename TYPE>
class CleverSideData : public SAMRAI::hier::PatchData
{
  public:
    static size_t getSizeOfData(
        const SAMRAI::hier::Box& box,
        int depth,
        const SAMRAI::hier::IntVector& ghosts);

    CleverSideData(
        const SAMRAI::hier::Box& box,
        int depth,
        const SAMRAI::hier::IntVector& ghosts);

    virtual ~CleverSideData();

    void copy(const SAMRAI::hier::PatchData& src);

    void copy2(SAMRAI::hier::PatchData& dst) const;

    void copy(const SAMRAI::hier::PatchData& src,
        const SAMRAI::hier::BoxOverlap& overlap);

    void copy2(SAMRAI::hier::PatchData& dst,
        const SAMRAI::hier::BoxOverlap& overlap) const;

    bool canEstimateStreamSizeFromBox() const;

    int getDataStreamSize(const SAMRAI::hier::BoxOverlap& overlap) const;

    void packStream(SAMRAI::tbox::MessageStream& stream,
        const SAMRAI::hier::BoxOverlap& overlap) const;

    void unpackStream(SAMRAI::tbox::MessageStream& stream,
        const SAMRAI::hier::BoxOverlap& overlap);

    int getDepth() const;

    TYPE* getPointer(int axis, int depth=0);

    void fill(const TYPE& value);

    void fillAll(const TYPE& value);
  private:
    boost::shared_ptr<AlignedArrayData<TYPE, 64> > d_array_data[3];
    int d_depth;
};

}
}

#include "CleverSideData.C"

#endif // CLEVERLEAF_PDAT_CLEVERSIDEDEDATA_H_

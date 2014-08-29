#ifndef CLEVERLEAF_PDAT_CLEVERCELLDATAFACTORY_H_
#define CLEVERLEAF_PDAT_CLEVERCELLDATAFACTORY_H_

#include "SAMRAI/hier/PatchDataFactory.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/BoxGeometry.h"
#include "SAMRAI/hier/Box.h"

namespace clever {
namespace pdat {

template<typename TYPE>
class CleverCellDataFactory : public SAMRAI::hier::PatchDataFactory
{
  public:
    CleverCellDataFactory(
        int depth,
        const SAMRAI::hier::IntVector& ghosts);

    virtual ~CleverCellDataFactory();

    boost::shared_ptr<SAMRAI::hier::PatchDataFactory> cloneFactory(
          const SAMRAI::hier::IntVector& ghosts);

    boost::shared_ptr<SAMRAI::hier::PatchData>
      allocate(
          const SAMRAI::hier::Patch& patch) const;

    boost::shared_ptr<SAMRAI::hier::BoxGeometry> getBoxGeometry(
          const SAMRAI::hier::Box& box) const;

    size_t getSizeOfMemory( const SAMRAI::hier::Box& box) const;

    bool fineBoundaryRepresentsVariable() const;

    bool dataLivesOnPatchBorder() const;

    bool validCopyTo(
          const boost::shared_ptr<SAMRAI::hier::PatchDataFactory>& dst_pdf)
      const;

    int getDepth();
  private:
    int d_depth;
};

}
}

#include "CleverCellDataFactory.C"

#endif // CLEVERLEAF_PDAT_CLEVERCELLDATAFACTORY_H_

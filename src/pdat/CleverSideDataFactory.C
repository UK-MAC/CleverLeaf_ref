#ifndef CLEVERLEAF_PDAT_CLEVERSIDEDATAFACTORY_C_
#define CLEVERLEAF_PDAT_CLEVERSIDEDATAFACTORY_C_

#include "CleverSideDataFactory.h"

#include "SAMRAI/tbox/MemoryUtilities.h"
#include "SAMRAI/pdat/SideGeometry.h"
#include <boost/make_shared.hpp>

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
boost::shared_ptr<SAMRAI::hier::PatchDataFactory> 
CleverSideDataFactory<TYPE>::cloneFactory(
    const SAMRAI::hier::IntVector& ghosts)
{
  return boost::make_shared<CleverSideDataFactory>(d_depth, ghosts);
}

template<typename TYPE>
boost::shared_ptr<SAMRAI::hier::PatchData> CleverSideDataFactory<TYPE>::allocate(
    const SAMRAI::hier::Patch& patch) const
{
  return boost::make_shared<CleverSideData<TYPE> >(
      patch.getBox(),
      d_depth,
      d_ghosts);
}

template<typename TYPE>
boost::shared_ptr<SAMRAI::hier::BoxGeometry>
CleverSideDataFactory<TYPE>::getBoxGeometry(const SAMRAI::hier::Box& box) const
{
  return boost::make_shared<SAMRAI::pdat::SideGeometry>(box, d_ghosts);
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
    const boost::shared_ptr<SAMRAI::hier::PatchDataFactory>& dst_pdf) const
{
  bool valid_copy = false;

  boost::shared_ptr<CleverSideDataFactory> side_data_factory(
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

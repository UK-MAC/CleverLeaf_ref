#ifndef CLEVERLEAF_PDAT_CLEVERCELLDATAFACTORY_C_
#define CLEVERLEAF_PDAT_CLEVERCELLDATAFACTORY_C_

#include "CleverCellDataFactory.h"

#include "SAMRAI/tbox/MemoryUtilities.h"
#include "SAMRAI/pdat/CellGeometry.h"
#include <boost/make_shared.hpp>

#include "CleverCellData.h"

namespace clever {
namespace pdat {

template<typename TYPE>
CleverCellDataFactory<TYPE>::CleverCellDataFactory(
    int depth,
    const SAMRAI::hier::IntVector& ghosts):
  SAMRAI::hier::PatchDataFactory(ghosts),
  d_depth(depth)
{
}

template<typename TYPE>
CleverCellDataFactory<TYPE>::~CleverCellDataFactory(){}

template<typename TYPE>
boost::shared_ptr<SAMRAI::hier::PatchDataFactory> 
CleverCellDataFactory<TYPE>::cloneFactory(
    const SAMRAI::hier::IntVector& ghosts)
{
  return boost::make_shared<CleverCellDataFactory>(d_depth, ghosts);
}

template<typename TYPE>
boost::shared_ptr<SAMRAI::hier::PatchData> CleverCellDataFactory<TYPE>::allocate(
    const SAMRAI::hier::Patch& patch) const
{
  return boost::make_shared<CleverCellData<TYPE> >(
      patch.getBox(),
      d_depth,
      d_ghosts);
}

template<typename TYPE>
boost::shared_ptr<SAMRAI::hier::BoxGeometry>
CleverCellDataFactory<TYPE>::getBoxGeometry(const SAMRAI::hier::Box& box) const
{
  return boost::make_shared<SAMRAI::pdat::CellGeometry>(box, d_ghosts);
}

template<typename TYPE>
size_t CleverCellDataFactory<TYPE>::getSizeOfMemory(const SAMRAI::hier::Box& box) const
{
  const size_t obj = 
    SAMRAI::tbox::MemoryUtilities::align(sizeof(CleverCellData<TYPE>));
  const size_t data = 
    CleverCellData<TYPE>::getSizeOfData(box, d_depth, d_ghosts);

  return obj + data;
}

template<typename TYPE>
bool CleverCellDataFactory<TYPE>::fineBoundaryRepresentsVariable() const
{
  return true;
}

template<typename TYPE>
bool CleverCellDataFactory<TYPE>::dataLivesOnPatchBorder() const
{
  return false;
}

template<typename TYPE>
bool CleverCellDataFactory<TYPE>::validCopyTo(
    const boost::shared_ptr<SAMRAI::hier::PatchDataFactory>& dst_pdf) const
{
  bool valid_copy = false;

  boost::shared_ptr<CleverCellDataFactory> cell_data_factory(
      dst_pdf,
      boost::detail::dynamic_cast_tag());

  if(cell_data_factory) {
    valid_copy = true;
  }

  return valid_copy;
}

template<typename TYPE>
int CleverCellDataFactory<TYPE>::getDepth()
{
  return d_depth;
}

}
}

#endif // CLEVERLEAF_PDAT_CLEVERCELLDATAFACTORY_C_

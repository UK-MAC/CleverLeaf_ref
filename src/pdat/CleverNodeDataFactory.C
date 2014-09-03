#ifndef CLEVERLEAF_PDAT_CLEVERNODEDATAFACTORY_C_
#define CLEVERLEAF_PDAT_CLEVERNODEDATAFACTORY_C_

#include "CleverNodeDataFactory.h"

#include "SAMRAI/tbox/MemoryUtilities.h"
#include "SAMRAI/pdat/NodeGeometry.h"
#include <boost/make_shared.hpp>

#include "CleverNodeData.h"

namespace clever {
namespace pdat {

template<typename TYPE>
CleverNodeDataFactory<TYPE>::CleverNodeDataFactory(
    int depth,
    const SAMRAI::hier::IntVector& ghosts):
  SAMRAI::hier::PatchDataFactory(ghosts),
  d_depth(depth)
{
}

template<typename TYPE>
CleverNodeDataFactory<TYPE>::~CleverNodeDataFactory(){}

template<typename TYPE>
boost::shared_ptr<SAMRAI::hier::PatchDataFactory> 
CleverNodeDataFactory<TYPE>::cloneFactory(
    const SAMRAI::hier::IntVector& ghosts)
{
  return boost::make_shared<CleverNodeDataFactory>(d_depth, ghosts);
}

template<typename TYPE>
boost::shared_ptr<SAMRAI::hier::PatchData> CleverNodeDataFactory<TYPE>::allocate(
    const SAMRAI::hier::Patch& patch) const
{
  return boost::make_shared<CleverNodeData<TYPE> >(
      patch.getBox(),
      d_depth,
      d_ghosts);
}

template<typename TYPE>
boost::shared_ptr<SAMRAI::hier::BoxGeometry>
CleverNodeDataFactory<TYPE>::getBoxGeometry(const SAMRAI::hier::Box& box) const
{
  return boost::make_shared<SAMRAI::pdat::NodeGeometry>(box, d_ghosts);
}

template<typename TYPE>
size_t CleverNodeDataFactory<TYPE>::getSizeOfMemory(const SAMRAI::hier::Box& box) const
{
  const size_t obj = 
    SAMRAI::tbox::MemoryUtilities::align(sizeof(CleverNodeData<TYPE>));
  const size_t data = 
    CleverNodeData<TYPE>::getSizeOfData(box, d_depth, d_ghosts);

  return obj + data;
}

template<typename TYPE>
bool CleverNodeDataFactory<TYPE>::fineBoundaryRepresentsVariable() const
{
  return true;
}

template<typename TYPE>
bool CleverNodeDataFactory<TYPE>::dataLivesOnPatchBorder() const
{
  return true;
}

template<typename TYPE>
bool CleverNodeDataFactory<TYPE>::validCopyTo(
    const boost::shared_ptr<SAMRAI::hier::PatchDataFactory>& dst_pdf) const
{
  bool valid_copy = false;

  boost::shared_ptr<CleverNodeDataFactory> node_data_factory(
      dst_pdf,
      boost::detail::dynamic_cast_tag());

  if(node_data_factory) {
    valid_copy = true;
  }

  return valid_copy;
}

template<typename TYPE>
int CleverNodeDataFactory<TYPE>::getDepth()
{
  return d_depth;
}

}
}

#endif // CLEVERLEAF_PDAT_CLEVERNODEDATAFACTORY_C_

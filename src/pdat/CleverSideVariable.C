#ifndef CLEVERLEAF_PDAT_CLEVERSIDEVARIABLE_C_
#define CLEVERLEAF_PDAT_CLEVERSIDEVARIABLE_C_

#include "CleverSideVariable.h"

#include "CleverSideDataFactory.h"

#include <boost/make_shared.hpp>

namespace clever {
namespace pdat {

template<typename TYPE>
CleverSideVariable<TYPE>::CleverSideVariable(
    const SAMRAI::tbox::Dimension& dim,
    const std::string& name,
    int depth):
  SAMRAI::hier::Variable(name,
      boost::make_shared<CleverSideDataFactory<TYPE> >(depth,
        SAMRAI::hier::IntVector::getZero(dim)))
{
}

template<typename TYPE>
CleverSideVariable<TYPE>::~CleverSideVariable(){}

template<typename TYPE>
bool CleverSideVariable<TYPE>::fineBoundaryRepresentsVariable() const
{
  return true;
}

template<typename TYPE>
bool CleverSideVariable<TYPE>::dataLivesOnPatchBorder() const
{
  return true;
}

template<typename TYPE>
int CleverSideVariable<TYPE>::getDepth() const
{
  boost::shared_ptr<CleverSideDataFactory<TYPE> > clever_side_data_factory(
      SHARED_PTR_CAST(CleverSideDataFactory<TYPE> ,
        getPatchDataFactory()));

  return clever_side_data_factory->getDepth();
}

}
}

#endif // CLEVERLEAF_PDAT_CLEVERSIDEVARIABLE_C_

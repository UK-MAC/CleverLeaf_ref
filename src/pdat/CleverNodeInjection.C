#ifndef CLEVERLEAF_PDAT_CLEVERNODEINJECTION_C_
#define CLEVERLEAF_PDAT_CLEVERNODEINJECTION_C_

#include "CleverNodeInjection.h"

#include "CleverNodeData.h"

#define F90_FUNC(name,NAME) name ## _

extern "C" {
  void F90_FUNC(conavgclevernodedoub2d, CONAVGCLEVERNODEDOUB2D)
    (const int&, const int&, const int&, const int&,
      const int&, const int&, const int&, const int&,
      const int&, const int&, const int&, const int&,
      const int *, const double *, double *);
}

namespace clever {
namespace pdat {

template<typename TYPE>
CleverNodeInjection<TYPE>::CleverNodeInjection():
  SAMRAI::hier::CoarsenOperator("CONSTANT_COARSEN")
{
}

template<typename TYPE>
CleverNodeInjection<TYPE>::~CleverNodeInjection()
{
}

template<typename TYPE>
int CleverNodeInjection<TYPE>::getOperatorPriority() const
{
  return 0;
}

template<typename TYPE>
SAMRAI::hier::IntVector CleverNodeInjection<TYPE>::getStencilWidth(
        const SAMRAI::tbox::Dimension &dim) const
{
  return SAMRAI::hier::IntVector::getZero(dim);
}

template<typename TYPE>
void CleverNodeInjection<TYPE>::coarsen(
    SAMRAI::hier::Patch& coarse,
    const SAMRAI::hier::Patch& fine,
    const int dst_component,
    const int src_component,
    const SAMRAI::hier::Box& coarse_box,
    const SAMRAI::hier::IntVector& ratio) const
{
  boost::shared_ptr<CleverNodeData<TYPE> > fine_data(
      fine.getPatchData(src_component),
      boost::detail::dynamic_cast_tag());

  boost::shared_ptr<CleverNodeData<TYPE> > coarse_data(
      coarse.getPatchData(dst_component),
      boost::detail::dynamic_cast_tag());

  const SAMRAI::hier::Box coarse_ghost_box =
    SAMRAI::pdat::NodeGeometry::toNodeBox(coarse_data->getGhostBox());
  const SAMRAI::hier::Box fine_ghost_box = 
    SAMRAI::pdat::NodeGeometry::toNodeBox(fine_data->getGhostBox());

   const hier::Index filo = fine_data->getGhostBox().lower();
   const hier::Index fihi = fine_data->getGhostBox().upper();
   const hier::Index cilo = coarse_data->getGhostBox().lower();
   const hier::Index cihi = coarse_data->getGhostBox().upper();

   const hier::Index ifirstc = coarse_box.lower();
   const hier::Index ilastc = coarse_box.upper();

  for(int depth = 0; depth < coarse_data->getDepth(); depth++) {
    if (fine.getDim() == SAMRAI::tbox::Dimension(2)) {
         F90_FUNC(conavgclevernodedoub2d, CONAVGCLEVERNODEDOUB2D)
           (ifirstc(0), ifirstc(1), ilastc(0), ilastc(1),
            filo(0), filo(1), fihi(0), fihi(1),
            cilo(0), cilo(1), cihi(0), cihi(1),
            &ratio[0],
            fine_data->getPointer(depth),
            coarse_data->getPointer(depth));
    } else {
      TBOX_ERROR("CleverNodeInjection error...\n"
          << "dim != 2 not supported." << std::endl);
    }
  }
}
}
}

#endif

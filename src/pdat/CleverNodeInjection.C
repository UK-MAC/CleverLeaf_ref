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
      SHARED_PTR_CAST(CleverNodeData<TYPE> ,
        fine.getPatchData(src_component)));

  boost::shared_ptr<CleverNodeData<TYPE> > coarse_data(
      SHARED_PTR_CAST(CleverNodeData<TYPE> ,
        coarse.getPatchData(dst_component)));

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

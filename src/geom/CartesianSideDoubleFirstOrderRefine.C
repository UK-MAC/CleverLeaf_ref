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
#include "CartesianSideDoubleFirstOrderRefine.h"

#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideOverlap.h"

#include "macros.h"

#define F90_FUNC(name,NAME) name ## _

extern "C" {
  void F90_FUNC(cartsidedoubfirstorderrefine0,
      CARTSIDEDOUBFIRSTORDERREFINE0)
    (const int&, const int&, const int&, const int&, const int&, const int&,
     const int&,const int&, const int&, const int&, const int&, const int&,
     const int&,const int&, const int&, const int&, const int*, const double*,
     const double*);

  void F90_FUNC(cartsidedoubfirstorderrefine1,
      CARTSIDEDOUBFIRSTORDERREFINE1)
    (const int&, const int&, const int&, const int&, const int&, const int&,
     const int&,const int&, const int&, const int&, const int&, const int&,
     const int&,const int&, const int&, const int&, const int*, const double*,
     const double*);
}

CartesianSideDoubleFirstOrderRefine::CartesianSideDoubleFirstOrderRefine():
  SAMRAI::hier::RefineOperator("FIRST_ORDER_REFINE")
{
}

CartesianSideDoubleFirstOrderRefine::~CartesianSideDoubleFirstOrderRefine()
{
}

int CartesianSideDoubleFirstOrderRefine::getOperatorPriority() const
{
  return 1;
}

SAMRAI::hier::IntVector CartesianSideDoubleFirstOrderRefine::getStencilWidth(
    const SAMRAI::tbox::Dimension &dim) const
{
  return SAMRAI::hier::IntVector::getZero(dim);
}

void CartesianSideDoubleFirstOrderRefine::refine(
    SAMRAI::hier::Patch& fine,
    const SAMRAI::hier::Patch& coarse,
    const int dst_component,
    const int src_component,
    const SAMRAI::hier::BoxOverlap& fine_overlap,
    const SAMRAI::hier::IntVector& ratio) const
{
  const SAMRAI::tbox::Dimension& dim(fine.getDim());

  boost::shared_ptr<SAMRAI::pdat::SideData<double> > cdata(
      SHARED_PTR_CAST(SAMRAI::pdat::SideData<double>,
        coarse.getPatchData(src_component)));

  boost::shared_ptr<SAMRAI::pdat::SideData<double> > fdata(
      SHARED_PTR_CAST(SAMRAI::pdat::SideData<double>,
        fine.getPatchData(dst_component)));

  const SAMRAI::pdat::SideOverlap* t_overlap =
    static_cast<const SAMRAI::pdat::SideOverlap *>(&fine_overlap);

  const SAMRAI::hier::Box cgbox(cdata->getGhostBox());

  const SAMRAI::hier::Index cilo = cgbox.lower();
  const SAMRAI::hier::Index cihi = cgbox.upper();
  const SAMRAI::hier::Index filo = fdata->getGhostBox().lower();
  const SAMRAI::hier::Index fihi = fdata->getGhostBox().upper();

  for (int axis = 0; axis < dim.getValue(); axis++) {
    const SAMRAI::hier::BoxContainer& boxes = 
      t_overlap->getDestinationBoxContainer(axis);

    for (SAMRAI::hier::BoxContainer::const_iterator b = boxes.begin();
        b != boxes.end(); ++b) {

      SAMRAI::hier::Box fine_box(*b);

      fine_box.upper(axis) -= 1;

      const SAMRAI::hier::Box coarse_box = SAMRAI::hier::Box::coarsen(fine_box, ratio);
      const SAMRAI::hier::Index ifirstc = coarse_box.lower();
      const SAMRAI::hier::Index ilastc = coarse_box.upper();
      const SAMRAI::hier::Index ifirstf = fine_box.lower();
      const SAMRAI::hier::Index ilastf = fine_box.upper();

      for (int d = 0; d < fdata->getDepth(); d++) {
        if (dim == SAMRAI::tbox::Dimension(2)) {
          if (axis == 0) {
            F90_FUNC(cartsidedoubfirstorderrefine0,
                CARTSIDEDOUBFIRSTORDERREFINE0)
              (ifirstc(0), ifirstc(1), ilastc(0), ilastc(1),
               ifirstf(0), ifirstf(1), ilastf(0), ilastf(1),
               filo(0), filo(1), fihi(0), fihi(1),
               cilo(0), cilo(1), cihi(0), cihi(1),
               &ratio[0],
               fdata->getPointer(0, d),
               cdata->getPointer(0, d));
          }
          if (axis == 1) {
            F90_FUNC(cartsidedoubfirstorderrefine1,
                CARTSIDEDOUBFIRSTORDERREFINE1)
              (ifirstc(0), ifirstc(1), ilastc(0), ilastc(1),
               ifirstf(0), ifirstf(1), ilastf(0), ilastf(1),
               filo(0), filo(1), fihi(0), fihi(1),
               cilo(0), cilo(1), cihi(0), cihi(1),
               &ratio[0],
               fdata->getPointer(1, d),
               cdata->getPointer(1, d));
          }
        } else {
          TBOX_ERROR( "CartesianSideFirstOrderRefine::refine dimension != 2 not supported"
              << std::endl);
        }
      }
    }
  }
}

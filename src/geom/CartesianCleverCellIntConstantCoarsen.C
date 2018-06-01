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
#include "CartesianCleverCellIntConstantCoarsen.h"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/VariableDatabase.h"

#include "pdat/CleverCellData.h"
#include "pdat/CleverCellVariable.h"

#include "macros.h"

#define F90_FUNC(name,NAME) name ## _

extern "C" {
  void F90_FUNC(cartesian_cell_constant_int_coarsen,
      CARTESIAN_CELL_CONSTANT_INT_COARSEN)
    (const int&,const int&,const int&,const int&,
     const int&,const int&,const int&,const int&,
     const int&,const int&,const int&,const int&,
     const int*, int*, const int*);
}

namespace clever {
namespace geom {

SAMRAI::tbox::StartupShutdownManager::Handler 
CartesianCleverCellIntConstantCoarsen::s_initialize_handler(
    CartesianCleverCellIntConstantCoarsen::initializeCallback,
    0,
    0,
    CartesianCleverCellIntConstantCoarsen::finalizeCallback,
    SAMRAI::tbox::StartupShutdownManager::priorityTimers);

std::shared_ptr<SAMRAI::tbox::Timer> CartesianCleverCellIntConstantCoarsen::t_coarsen;

CartesianCleverCellIntConstantCoarsen::CartesianCleverCellIntConstantCoarsen(
    const SAMRAI::tbox::Dimension& dim):
  SAMRAI::hier::CoarsenOperator("CONSTANT_INDICATOR_COARSEN")
{
}

CartesianCleverCellIntConstantCoarsen::~CartesianCleverCellIntConstantCoarsen()
{
}

bool CartesianCleverCellIntConstantCoarsen::findCoarsenOperator(
    const std::shared_ptr<SAMRAI::hier::Variable>& var,
    const std::string& op_name) const
{
  const std::shared_ptr<clever::pdat::CleverCellVariable<int> > cast_var(
      SHARED_PTR_CAST(clever::pdat::CleverCellVariable<int>, var));

  if (cast_var && (op_name == getOperatorName())) {
    return true;
  } else {
    return false;
  }
}

int CartesianCleverCellIntConstantCoarsen::getOperatorPriority() const
{
  return 0;
}

SAMRAI::hier::IntVector CartesianCleverCellIntConstantCoarsen::getStencilWidth(
    const SAMRAI::tbox::Dimension &dim) const 
{
  return SAMRAI::hier::IntVector::getZero(dim);
}

void CartesianCleverCellIntConstantCoarsen::coarsen(
    SAMRAI::hier::Patch& coarse,
    const SAMRAI::hier::Patch& fine,
    const int dst_component,
    const int src_component,
    const SAMRAI::hier::Box& coarse_box,
    const SAMRAI::hier::IntVector& ratio) const
{
  t_coarsen->start();
  const SAMRAI::tbox::Dimension& dim(fine.getDim());

  std::shared_ptr<clever::pdat::CleverCellData<int> > fdata(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<int>,
        fine.getPatchData(src_component)));
  std::shared_ptr<clever::pdat::CleverCellData<int> > cdata(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<int>,
        coarse.getPatchData(dst_component)));

  const SAMRAI::hier::Index filo = fdata->getGhostBox().lower();
  const SAMRAI::hier::Index fihi = fdata->getGhostBox().upper();
  const SAMRAI::hier::Index cilo = cdata->getGhostBox().lower();
  const SAMRAI::hier::Index cihi = cdata->getGhostBox().upper();

  const std::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> fgeom(
      SHARED_PTR_CAST(SAMRAI::geom::CartesianPatchGeometry,
        fine.getPatchGeometry()));
  const std::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> cgeom(
      SHARED_PTR_CAST(SAMRAI::geom::CartesianPatchGeometry,
        coarse.getPatchGeometry()));

  const SAMRAI::hier::Index ifirstc = coarse_box.lower();
  const SAMRAI::hier::Index ilastc = coarse_box.upper();

  for (int d = 0; d < cdata->getDepth(); d++) {
    if ((dim == SAMRAI::tbox::Dimension(2))) {

      F90_FUNC(cartesian_cell_constant_int_coarsen,
          CARTESIAN_CELL_CONSTANT_INT_COARSEN)
        (ifirstc(0),ifirstc(1),ilastc(0),ilastc(1),
         filo(0),filo(1),fihi(0),fihi(1),
         cilo(0),cilo(1),cihi(0),cihi(1),
         &ratio[0], fdata->getPointer(d), cdata->getPointer(d));

    } else {
      TBOX_ERROR("CartesianCleverCellIntConstantCoarsen error...\n"
          << "dim != 2 not supported." << std::endl);
    }
  }
  t_coarsen->stop();
}

void CartesianCleverCellIntConstantCoarsen::initializeCallback()
{
  t_coarsen = SAMRAI::tbox::TimerManager::getManager()->getTimer(
      "CartesianCleverCellIntConstantCoarsen::t_coarsen");
}

void CartesianCleverCellIntConstantCoarsen::finalizeCallback()
{
  t_coarsen.reset();
}

}
}

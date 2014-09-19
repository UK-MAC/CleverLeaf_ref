/*
 * Copyright 2013 David Beckingsale.
 * 
 * This file is part of CleverLeaf.
 * 
 * CleverLeaf is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 * 
 * CleverLeaf is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * CleverLeaf. If not, see http://www.gnu.org/licenses/.
 */ 
#include "CartesianCleverCellDoubleVolumeWeightedAverage.h"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Utilities.h"

#include "pdat/CleverCellData.h"
#include "pdat/CleverCellVariable.h"

#include "macros.h"

#define F90_FUNC(name,NAME) name ## _

extern "C" {
  void F90_FUNC(cartesian_cell_double_volume_weighted_coarsen,
      CARTESIAN_CELL_DOUBLE_VOLUME_WEIGHTED_COARSEN)
    (const int&, const int&, const int&, const int&, const int&, const int&,
     const int&, const int&, const int&, const int&, const int&, const int&,
     const int*, const double*, double*, const double*, const double*);
}

namespace clever {
namespace geom {

SAMRAI::tbox::StartupShutdownManager::Handler 
CartesianCleverCellDoubleVolumeWeightedAverage::s_initialize_handler(
    CartesianCleverCellDoubleVolumeWeightedAverage::initializeCallback,
    0,
    0,
    CartesianCleverCellDoubleVolumeWeightedAverage::finalizeCallback,
    SAMRAI::tbox::StartupShutdownManager::priorityTimers);

boost::shared_ptr<SAMRAI::tbox::Timer>
CartesianCleverCellDoubleVolumeWeightedAverage::t_coarsen;

CartesianCleverCellDoubleVolumeWeightedAverage::CartesianCleverCellDoubleVolumeWeightedAverage(
    const SAMRAI::tbox::Dimension& dim):
  SAMRAI::hier::CoarsenOperator("VOLUME_WEIGHTED_COARSEN")
{
}

CartesianCleverCellDoubleVolumeWeightedAverage::~CartesianCleverCellDoubleVolumeWeightedAverage()
{
}

bool CartesianCleverCellDoubleVolumeWeightedAverage::findCoarsenOperator(
    const boost::shared_ptr<SAMRAI::hier::Variable>& var,
    const std::string& op_name) const
{
  const boost::shared_ptr<clever::pdat::CleverCellVariable<double> > cast_var(
      SHARED_PTR_CAST(clever::pdat::CleverCellVariable<double>,
        var));

  if (cast_var && (op_name == getOperatorName())) {
    return true;
  } else {
    return false;
  }
}

int CartesianCleverCellDoubleVolumeWeightedAverage::getOperatorPriority() const
{
  return 1;
}

SAMRAI::hier::IntVector CartesianCleverCellDoubleVolumeWeightedAverage::getStencilWidth(
    const SAMRAI::tbox::Dimension &dim) const 
{
  return SAMRAI::hier::IntVector::getZero(dim);
}

void CartesianCleverCellDoubleVolumeWeightedAverage::coarsen(
    SAMRAI::hier::Patch& coarse,
    const SAMRAI::hier::Patch& fine,
    const int dst_component,
    const int src_component,
    const SAMRAI::hier::Box& coarse_box,
    const SAMRAI::hier::IntVector& ratio) const
{
  t_coarsen->start();

  const SAMRAI::tbox::Dimension& dim(fine.getDim());

  boost::shared_ptr<clever::pdat::CleverCellData<double> > fdata(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double>,
        fine.getPatchData(src_component)));
  boost::shared_ptr<clever::pdat::CleverCellData<double> > cdata(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double>,
        coarse.getPatchData(dst_component)));

  const SAMRAI::hier::Index filo = fdata->getGhostBox().lower();
  const SAMRAI::hier::Index fihi = fdata->getGhostBox().upper();
  const SAMRAI::hier::Index cilo = cdata->getGhostBox().lower();
  const SAMRAI::hier::Index cihi = cdata->getGhostBox().upper();

  const boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> fgeom(
      SHARED_PTR_CAST(SAMRAI::geom::CartesianPatchGeometry,
        fine.getPatchGeometry()));
  const boost::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> cgeom(
      SHARED_PTR_CAST(SAMRAI::geom::CartesianPatchGeometry,
        coarse.getPatchGeometry()));

  const SAMRAI::hier::Index ifirstc = coarse_box.lower();
  const SAMRAI::hier::Index ilastc = coarse_box.upper();

  const double* fdx = fgeom->getDx();
  const double* cdx = cgeom->getDx();

  for (int d = 0; d < cdata->getDepth(); d++) {
    double* farray = fdata->getPointer(d);
    double* carray = cdata->getPointer(d);

    if ((dim == SAMRAI::tbox::Dimension(2))) {
      double Vf = fdx[0]*fdx[1];
      double Vc = cdx[0]*cdx[1];

      F90_FUNC(cartesian_cell_double_volume_weighted_coarsen,
          CARTESIAN_CELL_DOUBLE_VOLUME_WEIGHTED_COARSEN)
        (ifirstc(0),
         ifirstc(1),
         ilastc(0),
         ilastc(1),
         filo(0),
         filo(1),
         fihi(0),
         fihi(1),
         cilo(0),
         cilo(1),
         cihi(0),
         cihi(1),
         &ratio[0],
         fdata->getPointer(d),
         cdata->getPointer(d),
         &Vf,
         &Vc);
    } else {
      TBOX_ERROR("CartesianCleverCellDoubleVolumeWeightedAverage error...\n"
          << "dim != 2 not supported." << std::endl);
    }
  }
  t_coarsen->stop();
}

void CartesianCleverCellDoubleVolumeWeightedAverage::initializeCallback()
{
  t_coarsen = SAMRAI::tbox::TimerManager::getManager()->getTimer(
      "CartesianCleverCellDoubleVolumeWeightedAverage::t_coarsen");
}

void CartesianCleverCellDoubleVolumeWeightedAverage::finalizeCallback()
{
  t_coarsen.reset();
}

}
}

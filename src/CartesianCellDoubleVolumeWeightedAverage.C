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
#include "CartesianCellDoubleVolumeWeightedAverage.h"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Utilities.h"

#define F90_FUNC(name,NAME) name ## _

extern "C" {
  void F90_FUNC(cartesian_cell_double_volume_weighted_coarsen,
      CARTESIAN_CELL_DOUBLE_VOLUME_WEIGHTED_COARSEN)
    (const int&, const int&, const int&, const int&, const int&, const int&,
     const int&, const int&, const int&, const int&, const int&, const int&,
     const int*, const double*, double*, const double*, const double*);
}

tbox::StartupShutdownManager::Handler 
CartesianCellDoubleVolumeWeightedAverage::s_initialize_handler(
    CartesianCellDoubleVolumeWeightedAverage::initializeCallback,
    0,
    0,
    CartesianCellDoubleVolumeWeightedAverage::finalizeCallback,
    tbox::StartupShutdownManager::priorityTimers);

boost::shared_ptr<tbox::Timer>
CartesianCellDoubleVolumeWeightedAverage::t_coarsen;

CartesianCellDoubleVolumeWeightedAverage::CartesianCellDoubleVolumeWeightedAverage(
    const tbox::Dimension& dim):
  hier::CoarsenOperator("VOLUME_WEIGHTED_COARSEN")
{
}

CartesianCellDoubleVolumeWeightedAverage::~CartesianCellDoubleVolumeWeightedAverage()
{
}

bool CartesianCellDoubleVolumeWeightedAverage::findCoarsenOperator(
    const boost::shared_ptr<SAMRAI::hier::Variable>& var,
    const std::string& op_name) const
{
  const boost::shared_ptr<pdat::CellVariable<double> > cast_var(
      var, boost::detail::dynamic_cast_tag());

  if (cast_var && (op_name == getOperatorName())) {
    return true;
  } else {
    return false;
  }
}

int CartesianCellDoubleVolumeWeightedAverage::getOperatorPriority() const
{
  return 1;
}

hier::IntVector CartesianCellDoubleVolumeWeightedAverage::getStencilWidth(
    const tbox::Dimension &dim) const 
{
  return hier::IntVector::getZero(dim);
}

void CartesianCellDoubleVolumeWeightedAverage::coarsen(
    hier::Patch& coarse,
    const hier::Patch& fine,
    const int dst_component,
    const int src_component,
    const hier::Box& coarse_box,
    const hier::IntVector& ratio) const
{
  t_coarsen->start();

  const tbox::Dimension& dim(fine.getDim());

  boost::shared_ptr<pdat::CellData<double> > fdata(
      fine.getPatchData(src_component), boost::detail::dynamic_cast_tag());
  boost::shared_ptr<pdat::CellData<double> > cdata(
      coarse.getPatchData(dst_component), boost::detail::dynamic_cast_tag());

  const hier::Index filo = fdata->getGhostBox().lower();
  const hier::Index fihi = fdata->getGhostBox().upper();
  const hier::Index cilo = cdata->getGhostBox().lower();
  const hier::Index cihi = cdata->getGhostBox().upper();

  const boost::shared_ptr<geom::CartesianPatchGeometry> fgeom(
      fine.getPatchGeometry(), boost::detail::dynamic_cast_tag());
  const boost::shared_ptr<geom::CartesianPatchGeometry> cgeom(
      coarse.getPatchGeometry(), boost::detail::dynamic_cast_tag());

  const hier::Index ifirstc = coarse_box.lower();
  const hier::Index ilastc = coarse_box.upper();

  const double* fdx = fgeom->getDx();
  const double* cdx = cgeom->getDx();

  for (int d = 0; d < cdata->getDepth(); d++) {
    double* farray = fdata->getPointer(d);
    double* carray = cdata->getPointer(d);

    if ((dim == tbox::Dimension(2))) {
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
      TBOX_ERROR("CartesianCellDoubleVolumeWeightedAverage error...\n"
          << "dim != 2 not supported." << std::endl);
    }
  }
  t_coarsen->stop();
}

void CartesianCellDoubleVolumeWeightedAverage::initializeCallback()
{
  t_coarsen = tbox::TimerManager::getManager()->getTimer(
      "CartesianCellDoubleVolumeWeightedAverage::t_coarsen");
}

void CartesianCellDoubleVolumeWeightedAverage::finalizeCallback()
{
  t_coarsen.reset();
}

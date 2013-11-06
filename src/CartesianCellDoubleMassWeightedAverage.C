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
#include "CartesianCellDoubleMassWeightedAverage.h"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/VariableDatabase.h"

#define F90_FUNC(name,NAME) name ## _

extern "C" {
  void F90_FUNC(cartesian_cell_double_mass_weighted_coarsen,
      CARTESIAN_CELL_DOUBLE_MASS_WEIGHTED_COARSEN)
    (const int&, const int&, const int&, const int&, const int&, const int&,
     const int&, const int&, const int&, const int&, const int&, const int&,
     const int*, const double*, double*, const double*, const double*,
     const double*, const double*);
}

tbox::StartupShutdownManager::Handler 
CartesianCellDoubleMassWeightedAverage::s_initialize_handler(
    CartesianCellDoubleMassWeightedAverage::initializeCallback,
    0,
    0,
    CartesianCellDoubleMassWeightedAverage::finalizeCallback,
    tbox::StartupShutdownManager::priorityTimers);

boost::shared_ptr<tbox::Timer>
CartesianCellDoubleMassWeightedAverage::t_coarsen;

CartesianCellDoubleMassWeightedAverage::CartesianCellDoubleMassWeightedAverage(
    const tbox::Dimension& dim):
  hier::CoarsenOperator("MASS_WEIGHTED_COARSEN")
{
}

CartesianCellDoubleMassWeightedAverage::~CartesianCellDoubleMassWeightedAverage()
{
}

bool CartesianCellDoubleMassWeightedAverage::findCoarsenOperator(
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

int CartesianCellDoubleMassWeightedAverage::getOperatorPriority() const
{
  /*
   * This priority is lower than the priority of the volume_weighted coarsen,
   * since we need to do that first.
   */
  return 2;
}

hier::IntVector
CartesianCellDoubleMassWeightedAverage::getStencilWidth(
    const tbox::Dimension &dim) const {
  return hier::IntVector::getZero(dim);
}

void CartesianCellDoubleMassWeightedAverage::coarsen(
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

  /*
   * We need to access the current density values in order to work out the mass
   * for our weighted coarsening.
   */
  hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

  int density_id = variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("density"), 
      variable_db->getContext("CURRENT"));

  boost::shared_ptr<pdat::CellData<double> > fmass(
      fine.getPatchData(density_id), boost::detail::dynamic_cast_tag());
  boost::shared_ptr<pdat::CellData<double> > cmass(
      coarse.getPatchData(density_id), boost::detail::dynamic_cast_tag());

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

    double* fmass_array = fmass->getPointer(d);
    double* cmass_array = cmass->getPointer(d);

    if ((dim == tbox::Dimension(2))) {
      double Vf = fdx[0]*fdx[1];
      double Vc = cdx[0]*cdx[1];

      F90_FUNC(cartesian_cell_double_mass_weighted_coarsen,
          CARTESIAN_CELL_DOUBLE_MASS_WEIGHTED_COARSEN)
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
         fmass->getPointer(d),
         cmass->getPointer(d),
         &Vf,
         &Vc);
    } else {
      TBOX_ERROR("CartesianCellDoubleMassWeightedAverage error...\n"
          << "dim != 2 not supported." << std::endl);
    }
  }
  t_coarsen->stop();
}

void CartesianCellDoubleMassWeightedAverage::initializeCallback()
{
  t_coarsen = tbox::TimerManager::getManager()->getTimer(
      "CartesianCellDoubleMassWeightedAverage::t_coarsen");
}

void CartesianCellDoubleMassWeightedAverage::finalizeCallback()
{
  t_coarsen.reset();
}

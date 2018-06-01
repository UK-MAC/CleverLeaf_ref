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
#include "CartesianCleverCellDoubleMassWeightedAverage.h"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/VariableDatabase.h"

#include "pdat/CleverCellData.h"
#include "pdat/CleverCellVariable.h"

#include "macros.h"

#define F90_FUNC(name,NAME) name ## _

extern "C" {
  void F90_FUNC(cartesian_cell_double_mass_weighted_coarsen,
      CARTESIAN_CELL_DOUBLE_MASS_WEIGHTED_COARSEN)
    (const int&, const int&, const int&, const int&, const int&, const int&,
     const int&, const int&, const int&, const int&, const int&, const int&,
     const int*, const double*, double*, const double*, const double*,
     const double*, const double*);
}

namespace clever {
namespace geom {

tbox::StartupShutdownManager::Handler 
CartesianCleverCellDoubleMassWeightedAverage::s_initialize_handler(
    CartesianCleverCellDoubleMassWeightedAverage::initializeCallback,
    0,
    0,
    CartesianCleverCellDoubleMassWeightedAverage::finalizeCallback,
    tbox::StartupShutdownManager::priorityTimers);

std::shared_ptr<tbox::Timer>
CartesianCleverCellDoubleMassWeightedAverage::t_coarsen;

CartesianCleverCellDoubleMassWeightedAverage::CartesianCleverCellDoubleMassWeightedAverage(
    const tbox::Dimension& dim):
  hier::CoarsenOperator("MASS_WEIGHTED_COARSEN")
{
}

CartesianCleverCellDoubleMassWeightedAverage::~CartesianCleverCellDoubleMassWeightedAverage()
{
}

bool CartesianCleverCellDoubleMassWeightedAverage::findCoarsenOperator(
    const std::shared_ptr<SAMRAI::hier::Variable>& var,
    const std::string& op_name) const
{
  const std::shared_ptr<pdat::CleverCellVariable<double> > cast_var(
      SHARED_PTR_CAST(pdat::CleverCellVariable<double>,
        var));

  if (cast_var && (op_name == getOperatorName())) {
    return true;
  } else {
    return false;
  }
}

int CartesianCleverCellDoubleMassWeightedAverage::getOperatorPriority() const
{
  /*
   * This priority is lower than the priority of the volume_weighted coarsen,
   * since we need to do that first.
   */
  return 2;
}

hier::IntVector
CartesianCleverCellDoubleMassWeightedAverage::getStencilWidth(
    const tbox::Dimension &dim) const {
  return hier::IntVector::getZero(dim);
}

void CartesianCleverCellDoubleMassWeightedAverage::coarsen(
    hier::Patch& coarse,
    const hier::Patch& fine,
    const int dst_component,
    const int src_component,
    const hier::Box& coarse_box,
    const hier::IntVector& ratio) const
{
  t_coarsen->start();
  const tbox::Dimension& dim(fine.getDim());

  std::shared_ptr<clever::pdat::CleverCellData<double> > fdata(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double>,
        fine.getPatchData(src_component)));
  std::shared_ptr<clever::pdat::CleverCellData<double> > cdata(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double>,
        coarse.getPatchData(dst_component)));

  /*
   * We need to access the current density values in order to work out the mass
   * for our weighted coarsening.
   */
  hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

  int density_id = variable_db->mapVariableAndContextToIndex(
      variable_db->getVariable("density"), 
      variable_db->getContext("CURRENT"));

  std::shared_ptr<clever::pdat::CleverCellData<double> > fmass(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double>,
        fine.getPatchData(density_id)));
  std::shared_ptr<clever::pdat::CleverCellData<double> > cmass(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<double>,
        coarse.getPatchData(density_id)));

  const hier::Index filo = fdata->getGhostBox().lower();
  const hier::Index fihi = fdata->getGhostBox().upper();
  const hier::Index cilo = cdata->getGhostBox().lower();
  const hier::Index cihi = cdata->getGhostBox().upper();

  const std::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> fgeom(
      SHARED_PTR_CAST(SAMRAI::geom::CartesianPatchGeometry,
        fine.getPatchGeometry()));
  const std::shared_ptr<SAMRAI::geom::CartesianPatchGeometry> cgeom(
      SHARED_PTR_CAST(SAMRAI::geom::CartesianPatchGeometry,
        coarse.getPatchGeometry()));

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
      TBOX_ERROR("CartesianCleverCellDoubleMassWeightedAverage error...\n"
          << "dim != 2 not supported." << std::endl);
    }
  }
  t_coarsen->stop();
}

void CartesianCleverCellDoubleMassWeightedAverage::initializeCallback()
{
  t_coarsen = tbox::TimerManager::getManager()->getTimer(
      "CartesianCleverCellDoubleMassWeightedAverage::t_coarsen");
}

void CartesianCleverCellDoubleMassWeightedAverage::finalizeCallback()
{
  t_coarsen.reset();
}

}
}

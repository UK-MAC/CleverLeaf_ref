#include "CartesianCellIntConstantCoarsen.h"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/VariableDatabase.h"

#define F90_FUNC(name,NAME) name ## _

extern "C" {
  void F90_FUNC(cartesian_cell_constant_int_coarsen,
      CARTESIAN_CELL_CONSTANT_INT_COARSEN)
    (const int&,const int&,const int&,const int&,
     const int&,const int&,const int&,const int&,
     const int&,const int&,const int&,const int&,
     const int*, int*, const int*);
}

tbox::StartupShutdownManager::Handler 
CartesianCellIntConstantCoarsen::s_initialize_handler(
    CartesianCellIntConstantCoarsen::initializeCallback,
    0,
    0,
    CartesianCellIntConstantCoarsen::finalizeCallback,
    tbox::StartupShutdownManager::priorityTimers);

boost::shared_ptr<tbox::Timer> CartesianCellIntConstantCoarsen::t_coarsen;

CartesianCellIntConstantCoarsen::CartesianCellIntConstantCoarsen(
    const tbox::Dimension& dim):
  hier::CoarsenOperator("CONSTANT_INDICATOR_COARSEN")
{
}

CartesianCellIntConstantCoarsen::~CartesianCellIntConstantCoarsen()
{
}

bool CartesianCellIntConstantCoarsen::findCoarsenOperator(
    const boost::shared_ptr<SAMRAI::hier::Variable>& var,
    const std::string& op_name) const
{
  const boost::shared_ptr<pdat::CellVariable<int> > cast_var(
      var, boost::detail::dynamic_cast_tag());

  if (cast_var && (op_name == getOperatorName())) {
    return true;
  } else {
    return false;
  }
}

int CartesianCellIntConstantCoarsen::getOperatorPriority() const
{
  return 0;
}

hier::IntVector CartesianCellIntConstantCoarsen::getStencilWidth(
    const tbox::Dimension &dim) const 
{
  return hier::IntVector::getZero(dim);
}

void CartesianCellIntConstantCoarsen::coarsen(
    hier::Patch& coarse,
    const hier::Patch& fine,
    const int dst_component,
    const int src_component,
    const hier::Box& coarse_box,
    const hier::IntVector& ratio) const
{
  t_coarsen->start();
  const tbox::Dimension& dim(fine.getDim());

  boost::shared_ptr<pdat::CellData<int> > fdata(
      fine.getPatchData(src_component), boost::detail::dynamic_cast_tag());
  boost::shared_ptr<pdat::CellData<int> > cdata(
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

  for (int d = 0; d < cdata->getDepth(); d++) {
    if ((dim == tbox::Dimension(2))) {

      F90_FUNC(cartesian_cell_constant_int_coarsen,
          CARTESIAN_CELL_CONSTANT_INT_COARSEN)
        (ifirstc(0),ifirstc(1),ilastc(0),ilastc(1),
         filo(0),filo(1),fihi(0),fihi(1),
         cilo(0),cilo(1),cihi(0),cihi(1),
         &ratio[0], fdata->getPointer(d), cdata->getPointer(d));

    } else {
      TBOX_ERROR("CartesianCellIntConstantCoarsen error...\n"
          << "dim != 2 not supported." << std::endl);
    }
  }
  t_coarsen->stop();
}

void CartesianCellIntConstantCoarsen::initializeCallback()
{
  t_coarsen = tbox::TimerManager::getManager()->getTimer(
      "CartesianCellIntConstantCoarsen::t_coarsen");
}

void CartesianCellIntConstantCoarsen::finalizeCallback()
{
  t_coarsen.reset();
}

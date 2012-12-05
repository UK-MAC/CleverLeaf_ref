#include "CartesianCellDoubleMassWeightedAverage.h"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Utilities.h"
#include "SAMRAI/hier/VariableDatabase.h"

#define POLY2(i, j, imin, jmin, nx) ((i - imin) + (j - jmin) * (nx))

CartesianCellDoubleMassWeightedAverage::CartesianCellDoubleMassWeightedAverage(
   const tbox::Dimension& dim):
   hier::CoarsenOperator(dim, "MASS_WEIGHTED_COARSEN")
{
}

CartesianCellDoubleMassWeightedAverage::~CartesianCellDoubleMassWeightedAverage()
{
}

bool CartesianCellDoubleMassWeightedAverage::findCoarsenOperator(
   const boost::shared_ptr<SAMRAI::hier::Variable>& var,
   const std::string& op_name) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *var);

   const boost::shared_ptr<pdat::CellVariable<double> > cast_var(var, boost::detail::dynamic_cast_tag());

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
CartesianCellDoubleMassWeightedAverage::getStencilWidth() const {
   return hier::IntVector::getZero(getDim());
}

void CartesianCellDoubleMassWeightedAverage::coarsen(
   hier::Patch& coarse,
   const hier::Patch& fine,
   const int dst_component,
   const int src_component,
   const hier::Box& coarse_box,
   const hier::IntVector& ratio) const
{
   const tbox::Dimension& dim(getDim());

   boost::shared_ptr<pdat::CellData<double> > fdata(fine.getPatchData(src_component), boost::detail::dynamic_cast_tag());
   boost::shared_ptr<pdat::CellData<double> > cdata(coarse.getPatchData(dst_component), boost::detail::dynamic_cast_tag());

   /*
    * We need the density to get the mass for this averaging...
    */

   hier::VariableDatabase* variable_db = hier::VariableDatabase::getDatabase();

   int density_id = variable_db->mapVariableAndContextToIndex(variable_db->getVariable("density"), variable_db->getContext("CURRENT"));

   boost::shared_ptr<pdat::CellData<double> > fmass(fine.getPatchData(density_id), boost::detail::dynamic_cast_tag());
   boost::shared_ptr<pdat::CellData<double> > cmass(coarse.getPatchData(density_id), boost::detail::dynamic_cast_tag());

   const hier::Index filo = fdata->getGhostBox().lower();
   const hier::Index fihi = fdata->getGhostBox().upper();
   const hier::Index cilo = cdata->getGhostBox().lower();
   const hier::Index cihi = cdata->getGhostBox().upper();

   const boost::shared_ptr<geom::CartesianPatchGeometry> fgeom(fine.getPatchGeometry(), boost::detail::dynamic_cast_tag());
   const boost::shared_ptr<geom::CartesianPatchGeometry> cgeom(coarse.getPatchGeometry(), boost::detail::dynamic_cast_tag());



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

          /*
           * cell volumes:
           */
          double Vf = fdx[0]*fdx[1];
          double Vc = cdx[0]*cdx[1];

          for(int k = ifirstc(1); k <= ilastc(1); k++) {
              for(int j = ifirstc(0); j <= ilastc(0); j++) {

                  /* Sigma eM */
                  double seM = 0.0;

                  for(int ry = 0; ry < ratio[1]; ry++) {
                      int kfine = k*ratio[1]+ry;

                      for(int rx = 0; rx < ratio[0]; rx++) {
                          int jfine = j*ratio[0]+rx;

                          /*
                           * energy*density*volume
                           */

                          seM += farray[POLY2(jfine,kfine,filo(0), filo(1), (fihi(0)-filo(0)+1))]*fmass_array[POLY2(jfine,kfine,filo(0), filo(1), (fihi(0)-filo(0)+1))]*Vf;

                      }
                  }

                  carray[POLY2(j,k,cilo(0), cilo(1), (cihi(0)-cilo(0)+1))] = seM/(cmass_array[POLY2(j,k,cilo(0), cilo(1), (cihi(0) - cilo(0) + 1))]*Vc);

//                  std::cout << "Updating energy coarse(" << j << "," << k << ") = " << carray[POLY2(j,k,cilo(0), cilo(1), (cihi(0)-cilo(0)+1))];

              }
          }

      } else {
         TBOX_ERROR("CartesianCellDoubleMassWeightedAverage error...\n"
            << "dim != 2 not supported." << std::endl);
      }
   }
}

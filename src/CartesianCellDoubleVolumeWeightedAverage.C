#include "CartesianCellDoubleVolumeWeightedAverage.h"

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/tbox/Utilities.h"

#define POLY2(i, j, imin, jmin, nx) ((i - imin) + (j - jmin) * (nx))

CartesianCellDoubleVolumeWeightedAverage::CartesianCellDoubleVolumeWeightedAverage(
   const tbox::Dimension& dim):
   hier::CoarsenOperator(dim, "VOLUME_WEIGHTED_COARSEN")
{
}

CartesianCellDoubleVolumeWeightedAverage::~CartesianCellDoubleVolumeWeightedAverage()
{
}

bool CartesianCellDoubleVolumeWeightedAverage::findCoarsenOperator(
   const tbox::Pointer<SAMRAI::hier::Variable>& var,
   const std::string& op_name) const
{
   TBOX_DIM_ASSERT_CHECK_ARGS2(*this, *var);

   const tbox::Pointer<pdat::CellVariable<double> > cast_var(var);

   if (!cast_var.isNull() && (op_name == getOperatorName())) {
      return true;
   } else {
      return false;
   }
}

int CartesianCellDoubleVolumeWeightedAverage::getOperatorPriority() const
{
   return 1;
}

hier::IntVector
CartesianCellDoubleVolumeWeightedAverage::getStencilWidth() const {
   return hier::IntVector::getZero(getDim());
}

void CartesianCellDoubleVolumeWeightedAverage::coarsen(
   hier::Patch& coarse,
   const hier::Patch& fine,
   const int dst_component,
   const int src_component,
   const hier::Box& coarse_box,
   const hier::IntVector& ratio) const
{
   const tbox::Dimension& dim(getDim());

   tbox::Pointer<pdat::CellData<double> > fdata = fine.getPatchData(src_component);
   tbox::Pointer<pdat::CellData<double> > cdata = coarse.getPatchData(dst_component);

   const hier::Index filo = fdata->getGhostBox().lower();
   const hier::Index fihi = fdata->getGhostBox().upper();
   const hier::Index cilo = cdata->getGhostBox().lower();
   const hier::Index cihi = cdata->getGhostBox().upper();

   const tbox::Pointer<geom::CartesianPatchGeometry> fgeom = fine.getPatchGeometry();
   const tbox::Pointer<geom::CartesianPatchGeometry> cgeom = coarse.getPatchGeometry();

   const hier::Index ifirstc = coarse_box.lower();
   const hier::Index ilastc = coarse_box.upper();

   const double* fdx = fgeom->getDx();
   const double* cdx = cgeom->getDx();

   for (int d = 0; d < cdata->getDepth(); d++) {

       double* farray = fdata->getPointer(d);
       double* carray = cdata->getPointer(d);

      if ((dim == tbox::Dimension(2))) {

          /*
           * cell volumes:
           */
          double Vf = fdx[0]*fdx[1];
          double Vc = cdx[0]*cdx[1];

          for(int k = ifirstc(1); k <= ilastc(1); k++) {
              for(int j = ifirstc(0); j <= ilastc(0); j++) {

                  /* Sigma pv */
                  double spv = 0.0;

                  for(int ry = 0; ry < ratio[1]; ry++) {
                      int kfine = k*ratio[1]+ry;

                      for(int rx = 0; rx < ratio[0]; rx++) {
                          int jfine = j*ratio[0]+rx;

 //                         std::cout << "jf,ky = " << jfine << "," << kfine << std::endl;
                          spv += farray[POLY2(jfine,kfine,filo(0), filo(1), (fihi(0)-filo(0)+1))]*Vf;
                      }
                  }

                  carray[POLY2(j,k,cilo(0), cilo(1), (cihi(0)-cilo(0)+1))] = spv/Vc;
//                  std::cout << "Updating coarse(" << j << "," << k << ") = " << carray[POLY2(j,k,cilo(0), cilo(1), (cihi(0)-cilo(0)+1))];

              }
          }

      } else {
         TBOX_ERROR("CartesianCellDoubleVolumeWeightedAverage error...\n"
            << "dim != 2 not supported." << std::endl);
      }
   }
}

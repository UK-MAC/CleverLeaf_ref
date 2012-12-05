#ifndef CARTESIANCELLMASSVOLUMEWEIGHTEDAVERAGE_H
#define CARTESIANCELLMASSVOLUMEWEIGHTEDAVERAGE_H

#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"

#include <string>

using namespace SAMRAI;

class CartesianCellDoubleMassWeightedAverage:
    public hier::CoarsenOperator
{
public:
   explicit CartesianCellDoubleMassWeightedAverage(
      const tbox::Dimension& dim);

   virtual ~CartesianCellDoubleMassWeightedAverage();

   bool
   findCoarsenOperator(
      const boost::shared_ptr<hier::Variable>& var,
      const std::string& op_name) const;

   int
   getOperatorPriority() const;

   hier::IntVector
   getStencilWidth() const;

   void
   coarsen(
      hier::Patch& coarse,
      const hier::Patch& fine,
      const int dst_component,
      const int src_component,
      const hier::Box& coarse_box,
      const hier::IntVector& ratio) const;
};

#endif

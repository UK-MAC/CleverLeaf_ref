#ifndef included_CartesianCellDoubleVolumeWeightedAverage
#define included_CartesianCellDoubleVolumeWeightedAverage

#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"

#include <string>

using namespace SAMRAI;

class CartesianCellDoubleVolumeWeightedAverage:
    public hier::CoarsenOperator
{
public:
   explicit CartesianCellDoubleVolumeWeightedAverage(
      const tbox::Dimension& dim);

   virtual ~CartesianCellDoubleVolumeWeightedAverage();

   bool
   findCoarsenOperator(
      const boost::shared_ptr<hier::Variable>& var,
      const std::string& op_name) const;

   int
   getOperatorPriority() const;

   hier::IntVector
   getStencilWidth( const tbox::Dimension &dim ) const;

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

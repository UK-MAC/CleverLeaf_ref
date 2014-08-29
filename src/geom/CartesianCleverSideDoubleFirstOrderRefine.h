#ifndef CLEVERLEAF_GEOM_CARTESIANCLEVERSIDEDOUBLEFIRSTORDERREFINE_H_
#define  CLEVERLEAF_GEOM_CARTESIANCLEVERSIDEDOUBLEFIRSTORDERREFINE_H_

#include "SAMRAI/hier/RefineOperator.h"

namespace clever {
namespace geom {

class CartesianCleverSideDoubleFirstOrderRefine 
  : public SAMRAI::hier::RefineOperator
{
  public:
    CartesianCleverSideDoubleFirstOrderRefine();
    virtual ~CartesianCleverSideDoubleFirstOrderRefine();

    int getOperatorPriority() const;

    SAMRAI::hier::IntVector getStencilWidth(
        const SAMRAI::tbox::Dimension &dim) const;

    void refine(
        SAMRAI::hier::Patch& fine,
        const SAMRAI::hier::Patch& coarse,
        const int dst_component,
        const int src_component,
        const SAMRAI::hier::BoxOverlap& fine_overlap,
        const SAMRAI::hier::IntVector& ratio) const;

  private:
    void refine(
        SAMRAI::hier::Patch& fine,
        const SAMRAI::hier::Patch& coarse,
        const int dst_component,
        const int src_component,
        const SAMRAI::hier::Box& fine_box,
        const SAMRAI::hier::IntVector& ratio) const;
};

}
}

#endif  // CLEVERLEAF_GEOM_CARTESIANCLEVERSIDEDOUBLEFIRSTORDERREFINE_H_

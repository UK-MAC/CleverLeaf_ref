#ifndef CLEVERLEAF_CARTESIANSIDEDOUBLEFIRSTORDERREFINE_H_
#define CLEVERLEAF_CARTESIANSIDEDOUBLEFIRSTORDERREFINE_H_

#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/BoxOverlap.h"
#include "SAMRAI/hier/IntVector.h"

class CartesianSideDoubleFirstOrderRefine:
  public SAMRAI::hier::RefineOperator
{
  public:
    CartesianSideDoubleFirstOrderRefine();

    virtual ~CartesianSideDoubleFirstOrderRefine();

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
};

#endif

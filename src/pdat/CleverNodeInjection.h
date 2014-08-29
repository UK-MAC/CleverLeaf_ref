#ifndef CLEVERLEAF_PDAT_CLEVERNODEINJECTION_H_
#define CLEVERLEAF_PDAT_CLEVERNODEINJECTION_H_

#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/Box.h"

namespace clever {
namespace pdat {

template<typename TYPE>
class CleverNodeInjection : public SAMRAI::hier::CoarsenOperator
{
  public:
    CleverNodeInjection();

    virtual ~CleverNodeInjection();

    int getOperatorPriority() const;

    SAMRAI::hier::IntVector getStencilWidth(
        const SAMRAI::tbox::Dimension &dim) const;

    void coarsen(
        SAMRAI::hier::Patch& coarse,
        const SAMRAI::hier::Patch& fine,
        const int dst_component,
        const int src_component,
        const SAMRAI::hier::Box& coarse_box,
        const SAMRAI::hier::IntVector& ratio) const;
};
}
}

#include "CleverNodeInjection.C"

#endif // CLEVERLEAF_PDAT_CLEVERNODEINJECTION_H_

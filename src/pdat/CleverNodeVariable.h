#ifndef CLEVERLEAF_PDAT_CLEVERNODEVARIABLE_H_
#define CLEVERLEAF_PDAT_CLEVERNODEVARIABLE_H_

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/appu/VisDerivedDataStrategy.h"

#include "boost/enable_shared_from_this.hpp"

namespace clever {
namespace pdat {

template<typename TYPE>
class CleverNodeVariable : public SAMRAI::hier::Variable,
  public SAMRAI::appu::VisDerivedDataStrategy,
  public boost::enable_shared_from_this<CleverNodeVariable<TYPE> >
{
  public:
    CleverNodeVariable(
        const SAMRAI::tbox::Dimension& dim,
        const std::string& name,
        int depth = 1);

    virtual ~CleverNodeVariable();

    bool fineBoundaryRepresentsVariable() const;
    
    bool dataLivesOnPatchBorder() const;

    int getDepth() const;

    bool packDerivedDataIntoDoubleBuffer(
        double* buffer,
        const SAMRAI::hier::Patch& patch,
        const SAMRAI::hier::Box& region,
        const std::string& variable_name,
        int depth_index) const;
};

}
}

#include "CleverNodeVariable.C"

#endif // CLEVERLEAF_PDAT_CLEVERNODEVARIABLE_H_

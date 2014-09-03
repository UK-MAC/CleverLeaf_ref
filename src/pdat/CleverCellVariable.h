#ifndef CLEVERLEAF_PDAT_CLEVERCELLVARIABLE_H_
#define CLEVERLEAF_PDAT_CLEVERCELLVARIABLE_H_

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/appu/VisDerivedDataStrategy.h"

#include "boost/enable_shared_from_this.hpp"

namespace clever {
namespace pdat {

template<typename TYPE>
class CleverCellVariable : public SAMRAI::hier::Variable,
  public SAMRAI::appu::VisDerivedDataStrategy,
  public boost::enable_shared_from_this<CleverCellVariable<TYPE> >
{
  public:
    CleverCellVariable(
        const SAMRAI::tbox::Dimension& dim,
        const std::string& name,
        int depth = 1);

    virtual ~CleverCellVariable();

    bool fineBoundaryRepresentsVariable() const;
    
    bool dataLivesOnPatchBorder() const;

    int getDepth() const;

    bool packDerivedDataIntoDoubleBuffer(
        double* buffer,
        const SAMRAI::hier::Patch& patch,
        const SAMRAI::hier::Box& region,
        const std::string& variable_name,
        int depth_index) const;
  private:
};

}
}

#include "CleverCellVariable.C"

#endif // CLEVERLEAF_PDAT_CLEVERCELLVARIABLE_H_

#ifndef CLEVERLEAF_PDAT_CLEVERSIDEVARIABLE_H_
#define CLEVERLEAF_PDAT_CLEVERSIDEVARIABLE_H_

#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/tbox/Dimension.h"

namespace clever {
namespace pdat {

template<typename TYPE>
class CleverSideVariable : public SAMRAI::hier::Variable
{
  public:
    CleverSideVariable(
        const SAMRAI::tbox::Dimension& dim,
        const std::string& name,
        int depth = 1);

    virtual ~CleverSideVariable();

    bool fineBoundaryRepresentsVariable() const;
    
    bool dataLivesOnPatchBorder() const;

    int getDepth() const;
};

}
}

#include "CleverSideVariable.C"

#endif // CLEVERLEAF_PDAT_CLEVERSIDEVARIABLE_H_

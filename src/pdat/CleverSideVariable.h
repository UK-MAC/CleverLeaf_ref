//////////////////////////////////////////////////////////////////////////////
// Crown Copyright 2014 AWE, Copyright 2014 David Beckingsale.
//
// This file is part of CleverLeaf.
//
// CleverLeaf is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.
//
// CleverLeaf is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with
// CleverLeaf. If not, see http://www.gnu.org/licenses/.
//////////////////////////////////////////////////////////////////////////////
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

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
#ifndef CLEVERLEAF_PDAT_CLEVERNODEVARIABLE_H_
#define CLEVERLEAF_PDAT_CLEVERNODEVARIABLE_H_

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/appu/VisDerivedDataStrategy.h"

//#include "boost/enable_shared_from_this.hpp"

namespace clever {
namespace pdat {

template<typename TYPE>
class CleverNodeVariable : public SAMRAI::hier::Variable,
  public SAMRAI::appu::VisDerivedDataStrategy,
  public std::enable_shared_from_this<CleverNodeVariable<TYPE> >
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
        int depth_index,
        double simulation_time) const;
};

}
}

#include "CleverNodeVariable.C"

#endif // CLEVERLEAF_PDAT_CLEVERNODEVARIABLE_H_

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

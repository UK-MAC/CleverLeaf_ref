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

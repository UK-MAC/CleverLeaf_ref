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
#ifndef CLEVERLEAF_PDAT_CLEVERCELLDATAFACTORY_H_
#define CLEVERLEAF_PDAT_CLEVERCELLDATAFACTORY_H_

#include "SAMRAI/hier/PatchDataFactory.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/BoxGeometry.h"
#include "SAMRAI/hier/Box.h"

namespace clever {
namespace pdat {

template<typename TYPE>
class CleverCellDataFactory : public SAMRAI::hier::PatchDataFactory
{
  public:
    CleverCellDataFactory(
        int depth,
        const SAMRAI::hier::IntVector& ghosts);

    virtual ~CleverCellDataFactory();

    std::shared_ptr<SAMRAI::hier::PatchDataFactory> cloneFactory(
          const SAMRAI::hier::IntVector& ghosts);

    std::shared_ptr<SAMRAI::hier::PatchData>
      allocate(
          const SAMRAI::hier::Patch& patch) const;

    std::shared_ptr<SAMRAI::hier::BoxGeometry> getBoxGeometry(
          const SAMRAI::hier::Box& box) const;

    size_t getSizeOfMemory( const SAMRAI::hier::Box& box) const;

    bool fineBoundaryRepresentsVariable() const;

    bool dataLivesOnPatchBorder() const;

    bool validCopyTo(
          const std::shared_ptr<SAMRAI::hier::PatchDataFactory>& dst_pdf)
      const;

    int getDepth();
  private:
    int d_depth;
};

}
}

#include "CleverCellDataFactory.C"

#endif // CLEVERLEAF_PDAT_CLEVERCELLDATAFACTORY_H_

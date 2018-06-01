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
#ifndef CLEVERLEAF_PDAT_CLEVERCELLVARIABLE_C_
#define CLEVERLEAF_PDAT_CLEVERCELLVARIABLE_C_

#include "CleverCellVariable.h"


#include "SAMRAI/hier/VariableDatabase.h"

#include "CleverCellDataFactory.h"

namespace clever {
namespace pdat {

template<typename TYPE>
CleverCellVariable<TYPE>::CleverCellVariable(
    const SAMRAI::tbox::Dimension& dim,
    const std::string& name,
    int depth):
  SAMRAI::hier::Variable(name,
      std::make_shared<CleverCellDataFactory<TYPE> >(depth,
        SAMRAI::hier::IntVector::getZero(dim)))
{
}

template<typename TYPE>
CleverCellVariable<TYPE>::~CleverCellVariable(){}

template<typename TYPE>
bool CleverCellVariable<TYPE>::fineBoundaryRepresentsVariable() const
{
  return true;
}

template<typename TYPE>
bool CleverCellVariable<TYPE>::dataLivesOnPatchBorder() const
{
  return false;
}

template<typename TYPE>
int CleverCellVariable<TYPE>::getDepth() const
{
  std::shared_ptr<CleverCellDataFactory<TYPE> > clever_cell_data_factory(
      SHARED_PTR_CAST(CleverCellDataFactory<TYPE> ,
        getPatchDataFactory()));
  return clever_cell_data_factory->getDepth();
}

template<typename TYPE>
bool CleverCellVariable<TYPE>::packDerivedDataIntoDoubleBuffer(
    double* buffer,
    const SAMRAI::hier::Patch& patch,
    const SAMRAI::hier::Box& region,
    const std::string& variable_name,
    int depth_index,
    double simulation_time) const
{
  const SAMRAI::tbox::Dimension& dim(patch.getDim());

  SAMRAI::hier::VariableDatabase* variable_db = 
    SAMRAI::hier::VariableDatabase::getDatabase();

  std::shared_ptr<clever::pdat::CleverCellData<TYPE> > cell_data(
      SHARED_PTR_CAST(clever::pdat::CleverCellData<TYPE> ,
        patch.getPatchData(std::const_pointer_cast<clever::pdat::CleverCellVariable<TYPE> >(this->shared_from_this()),
          variable_db->getContext("CURRENT"))));

  bool data_on_patch = false;

  TYPE* data = cell_data->getPointer();

  const SAMRAI::hier::Box& cell_data_box = cell_data->getGhostBox();

  const int box_width = region.numberCells(0);
  const int box_height = region.numberCells(1);
  const int cell_data_width = cell_data_box.numberCells(0);
  const int cell_data_height = cell_data_box.numberCells(1);
  
  if (dim == SAMRAI::tbox::Dimension(2)) {
    int buffer_offset = 0;
    int data_offset = cell_data_box.offset(region.lower());

    for (int i1 = 0; i1 < box_height; i1++) {
      for (int i0 = 0; i0 < box_width; i0++) {
        int data_index = data_offset + i0;
        buffer[buffer_offset + i0] = static_cast<double>(data[data_index]);
      }
      data_offset += cell_data_width;
      buffer_offset += box_width;
    }
    data_on_patch = true;
  }

  return data_on_patch;
}

}
}

#endif // CLEVERLEAF_PDAT_CLEVERCELLVARIABLE_C_

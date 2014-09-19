#ifndef CLEVERLEAF_PDAT_CLEVERNODEVARIABLE_C_
#define CLEVERLEAF_PDAT_CLEVERNODEVARIABLE_C_

#include "CleverNodeVariable.h"


#include <boost/make_shared.hpp>
#include "SAMRAI/hier/VariableDatabase.h"

#include "CleverNodeDataFactory.h"

namespace clever {
namespace pdat {

template<typename TYPE>
CleverNodeVariable<TYPE>::CleverNodeVariable(
    const SAMRAI::tbox::Dimension& dim,
    const std::string& name,
    int depth):
  SAMRAI::hier::Variable(name,
      boost::make_shared<CleverNodeDataFactory<TYPE> >(depth,
        SAMRAI::hier::IntVector::getZero(dim)))
{
}

template<typename TYPE>
CleverNodeVariable<TYPE>::~CleverNodeVariable(){}

template<typename TYPE>
bool CleverNodeVariable<TYPE>::fineBoundaryRepresentsVariable() const
{
  return true;
}

template<typename TYPE>
bool CleverNodeVariable<TYPE>::dataLivesOnPatchBorder() const
{
  return true;
}

template<typename TYPE>
int CleverNodeVariable<TYPE>::getDepth() const
{
  boost::shared_ptr<CleverNodeDataFactory<TYPE> > clever_node_data_factory(
      SHARED_PTR_CAST(CleverNodeDataFactory<TYPE> ,
        getPatchDataFactory()));

  return clever_node_data_factory->getDepth();
}

template<typename TYPE>
bool CleverNodeVariable<TYPE>::packDerivedDataIntoDoubleBuffer(
    double* buffer,
    const SAMRAI::hier::Patch& patch,
    const SAMRAI::hier::Box& region,
    const std::string& variable_name,
    int depth_index) const
{
  const SAMRAI::tbox::Dimension& dim(patch.getDim());

  SAMRAI::hier::VariableDatabase* variable_db = 
    SAMRAI::hier::VariableDatabase::getDatabase();

  boost::shared_ptr<clever::pdat::CleverNodeData<TYPE> > node_data(
      SHARED_PTR_CAST(clever::pdat::CleverNodeData<TYPE>,
        patch.getPatchData(boost::const_pointer_cast<clever::pdat::CleverNodeVariable<TYPE> >(this->shared_from_this()),
          variable_db->getContext("CURRENT"))));

  bool data_on_patch = false;

  TYPE* data = node_data->getPointer();

  const hier::Box& data_box = SAMRAI::pdat::NodeGeometry::toNodeBox(
      node_data->getGhostBox());

  const hier::Box& node_region = SAMRAI::pdat::NodeGeometry::toNodeBox(region);

  const int box_width = node_region.numberCells(0);
  const int box_height = node_region.numberCells(1);
  const int data_width = data_box.numberCells(0);
  const int data_height = data_box.numberCells(1);
  
  if (dim == tbox::Dimension(2)) {
    int buffer_offset = 0;
    int data_offset = data_box.offset(region.lower())
      + (data_box.size()*depth_index);

    for (int i1 = 0; i1 < box_height; i1++) {
      for (int i0 = 0; i0 < box_width; i0++) {
        int data_index = data_offset + i0;

        buffer[buffer_offset + i0] = data[data_index];
      }

      data_offset += data_width;
      buffer_offset += box_width;
    }

    data_on_patch = true;
  }

  return data_on_patch;
}

}
}

#endif // CLEVERLEAF_PDAT_CLEVERNODEVARIABLE_C_

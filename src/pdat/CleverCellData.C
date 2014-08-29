#ifndef CLEVERLEAF_PDAT_CLEVERCELLDATA_C_
#define CLEVERLEAF_PDAT_CLEVERCELLDATA_C_

#include "CleverCellData.h"

#include "SAMRAI/pdat/CellOverlap.h"

namespace clever {
namespace pdat {

template<typename TYPE>
size_t CleverCellData<TYPE>::getSizeOfData(
    const SAMRAI::hier::Box& box,
    int depth,
    const SAMRAI::hier::IntVector& ghosts)
{
  const SAMRAI::hier::Box ghost_box = SAMRAI::hier::Box::grow(box, ghosts);
  return AlignedArrayData<TYPE, 64>::getSizeOfData(ghost_box, depth);
}

template<typename TYPE>
CleverCellData<TYPE>::CleverCellData(
    const SAMRAI::hier::Box& domain,
    int depth,
    const SAMRAI::hier::IntVector& ghosts):
  SAMRAI::hier::PatchData(domain,ghosts),
  d_depth(depth)
{
  d_array_data.reset(new AlignedArrayData<TYPE, 64>(getGhostBox(), depth));
}

template<typename TYPE>
CleverCellData<TYPE>::~CleverCellData(){}

template<typename TYPE>
TYPE* CleverCellData<TYPE>::getPointer(int depth)
{
  return d_array_data->getPointer(depth);
}

template<typename TYPE>
void CleverCellData<TYPE>::fill(const TYPE& value)
{
  d_array_data->fill(value);
}

template<typename TYPE>
void CleverCellData<TYPE>::fillAll(const TYPE& value)
{
  d_array_data->fillAll(value);
}

template<typename TYPE>
void CleverCellData<TYPE>::copy(const SAMRAI::hier::PatchData& src)
{
  const CleverCellData* cell_data_source = 
    dynamic_cast<const CleverCellData<TYPE> *>(&src);

  if (cell_data_source == 0) {
    src.copy2(*this);
  } else {
    const SAMRAI::hier::Box box = 
      d_array_data->getBox() * cell_data_source->d_array_data->getBox();
    if (!box.empty()) {
      d_array_data->copy(*(cell_data_source->d_array_data), box);
    }
  }
}

template<typename TYPE>
void CleverCellData<TYPE>::copy2(SAMRAI::hier::PatchData& dst) const
{
  CleverCellData* cell_data_destination = 
   static_cast<CleverCellData<TYPE> *>(&dst); 

  const SAMRAI::hier::Box box = 
    d_array_data->getBox()*cell_data_destination->d_array_data->getBox();

  if (!box.empty()) {
    cell_data_destination->d_array_data->copy(*d_array_data, box);
  }
}

template<typename TYPE>
void CleverCellData<TYPE>::copy(
    const SAMRAI::hier::PatchData& src,
    const SAMRAI::hier::BoxOverlap& overlap)
{
  const CleverCellData* cell_data_source = 
    dynamic_cast<const CleverCellData<TYPE> *>(&src);

  const SAMRAI::pdat::CellOverlap* cell_overlap = 
    dynamic_cast<const SAMRAI::pdat::CellOverlap *>(&overlap);

  if ((cell_data_source == 0) || (cell_overlap == 0)) {
    src.copy2(*this, overlap);
  } else {
    if (cell_overlap->getTransformation().getRotation() ==
        SAMRAI::hier::Transformation::NO_ROTATE) {
      d_array_data->copy(*(cell_data_source->d_array_data),
          cell_overlap->getDestinationBoxContainer(),
          cell_overlap->getTransformation());
    } else {
      std::cerr << "CleverCellData: rotated boxes not supported in copy" 
        << std::endl;
      std::abort();
    }
  }
}

template<typename TYPE>
void CleverCellData<TYPE>::copy2(
    SAMRAI::hier::PatchData& dst, 
    const SAMRAI::hier::BoxOverlap& overlap) const
{
  CleverCellData* cell_data_destination =
    static_cast<CleverCellData<TYPE> *>(&dst);

  const SAMRAI::pdat::CellOverlap* cell_overlap =
    static_cast<const SAMRAI::pdat::CellOverlap *>(&overlap);

  if (cell_overlap->getTransformation().getRotation() ==
      SAMRAI::hier::Transformation::NO_ROTATE) {
    cell_data_destination->d_array_data->copy(*d_array_data,
        cell_overlap->getDestinationBoxContainer(),
        cell_overlap->getTransformation());
  } else {
      std::cerr << "CleverCellData: rotated boxes not supported in copy2"
        << std::endl;
      std::abort();
  }
}

template<typename TYPE>
bool CleverCellData<TYPE>::canEstimateStreamSizeFromBox() const
{
    return true;
}

template<typename TYPE>
int CleverCellData<TYPE>::getDataStreamSize(
    const SAMRAI::hier::BoxOverlap& overlap) const
{
  const SAMRAI::pdat::CellOverlap* cell_overlap = 
    static_cast<const SAMRAI::pdat::CellOverlap *>(&overlap);

  return d_array_data->getDataStreamSize(
      cell_overlap->getDestinationBoxContainer());
}

template<typename TYPE>
void CleverCellData<TYPE>::packStream(
    SAMRAI::tbox::MessageStream& stream,
        const SAMRAI::hier::BoxOverlap& overlap) const
{
  const SAMRAI::pdat::CellOverlap* cell_overlap = 
    static_cast<const SAMRAI::pdat::CellOverlap *>(&overlap);

  if(cell_overlap->getTransformation().getRotation()
      == SAMRAI::hier::Transformation::NO_ROTATE)
  {
    d_array_data->packStream(
        stream,
        cell_overlap->getDestinationBoxContainer(),
        cell_overlap->getTransformation());
  } else {
    std::cerr << "CleverCellData: rotated boxes not supported in packStream!" 
      << std::endl;
    std::abort();
  }
}

template<typename TYPE>
void CleverCellData<TYPE>::unpackStream(SAMRAI::tbox::MessageStream& stream,
    const SAMRAI::hier::BoxOverlap& overlap)
{
  const SAMRAI::pdat::CellOverlap* cell_overlap = 
    static_cast<const SAMRAI::pdat::CellOverlap *>(&overlap);

  d_array_data->unpackStream(
      stream,
      cell_overlap->getDestinationBoxContainer(),
      cell_overlap->getSourceOffset());
}

template<typename TYPE>
int CleverCellData<TYPE>::getDepth() const
{
  return d_depth;
}

}
}

#endif // CLEVERLEAF_PDAT_CLEVERCELLDATA_C_

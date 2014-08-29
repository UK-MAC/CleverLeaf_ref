#ifndef CLEVERLEAF_PDAT_ALIGNEDARRAYDATAOPERATIONUTILITIES_H_
#define CLEVERLEAF_PDAT_ALIGNEDARRAYDATAOPERATIONUTILITIES_H_

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"

namespace clever {
namespace pdat {

template<typename TYPE, int ALIGNMENT> class AlignedArrayData;

/*!
 * @brief Class ArrayDataOperationUtilities<TYPE, OP> provides
 * generic looping operations for all array-based patch data types.
 * The operations are templated on data type
 * and the operation that will be performed on individual array
 * elements in the innermost loop.
 *
 * @see ArrayData
 */

template<typename TYPE, int ALIGNMENT>
class AlignedArrayDataOperationUtilities
{
public:
   /*!
    * Perform operation on a subset of data components of source and
    * destination array data objects and put results in destination array data
    * object.
    *
    * @param dst    Reference to destination array data object.
    * @param src    Const reference to source array data object.
    * @param opbox  Const reference to Box indicating index space region of
    *               operation.
    * @param src_shift  Const reference to IntVector indicating shift required
    *                   to put source index space region into destination index
    *                   space region.
    * @param dst_start_depth  Integer specifying starting depth component of
    *                         operation in destination array.
    * @param src_start_depth  Integer specifying starting depth component of
    *                         operation in source array.
    * @param num_depth  Integer number of depth components on which to perform
    *                   operation.
    * @param op  Const reference to object that performs operations on
    *            individual data array elements.
    *
    * @pre (dst.getDim() == src.getDim()) &&
    *      (dst.getDim() == opbox.getDim()) &&
    *      (dst.getDim() == src_shift.getDim())
    * @pre num_depth >= 0
    * @pre (0 <= dst_start_depth) &&
    *      (dst_start_depth + num_depth <= dst.getDepth())
    * @pre (0 <= src_start_depth) &&
    *      (src_start_depth + num_depth <= src.getDepth())
    */
   static void
   doArrayDataOperationOnBox(
      AlignedArrayData<TYPE, ALIGNMENT>& dst,
      const AlignedArrayData<TYPE, ALIGNMENT>& src,
      const SAMRAI::hier::Box& opbox,
      const SAMRAI::hier::IntVector& src_shift,
      int dst_start_depth,
      int src_start_depth,
      int num_depth);

   /*!
    * Perform operation on all data components of array data object and
    * corresponding buffer data, putting results in either the array data
    * object or buffer.
    *
    * @param arraydata   Const reference to array data object.
    * @param buffer      Const pointer to first element in buffer.
    * @param opbox       Const reference to Box indicating operation region in
    *                    index space of array data object.
    * @param src_is_buffer  Boolean value indicating whether buffer is source
    *                       data for operation; if true results will be placed
    *                       in array data object, otherwise results will go in
    *                       buffer.
    * @param op  Const reference to object that performs operations on
    *            individual data array elements.
    *
    * @pre arraydata.getDim() == opbox.getDim()
    * @pre buffer != 0
    * @pre opbox.isSpatiallyEqual(opbox * arraydata.getBox())
    */
   template<bool src_is_buffer>
   static void
   doArrayDataBufferOperationOnBox(
      const AlignedArrayData<TYPE, ALIGNMENT>& arraydata,
      const TYPE * buffer,
      const SAMRAI::hier::Box& opbox);
private:
   // the following are not implemented:
   AlignedArrayDataOperationUtilities();
   ~AlignedArrayDataOperationUtilities();
   AlignedArrayDataOperationUtilities(
       const AlignedArrayDataOperationUtilities&);
   AlignedArrayDataOperationUtilities& operator = (
      const AlignedArrayDataOperationUtilities&);
};

}
}

#include "AlignedArrayDataOperationUtilities.C"

#endif // CLEVERLEAF_PDAT_ALIGNEDARRAYDATAOPERATIONUTILITIES_H_

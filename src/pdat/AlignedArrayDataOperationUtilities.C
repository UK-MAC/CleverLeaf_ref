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
#ifndef CLEVERLEAF_PDAT_ALIGNEDARRAYDATAOPERATIONUTILITIES_C_
#define CLEVERLEAF_PDAT_ALIGNEDARRAYDATAOPERATIONUTILITIES_C_

#include "AlignedArrayDataOperationUtilities.h"
#include "AlignedArrayData.h"

#include "SAMRAI/tbox/Dimension.h"

#define MAX_DIMENSION 3

namespace clever {
namespace pdat {

template<typename TYPE, int ALIGNMENT>
void AlignedArrayDataOperationUtilities<TYPE, ALIGNMENT>::doArrayDataOperationOnBox(
   AlignedArrayData<TYPE, ALIGNMENT>& dst,
   const AlignedArrayData<TYPE, ALIGNMENT>& src,
   const SAMRAI::hier::Box& opbox,
   const SAMRAI::hier::IntVector& src_shift,
   int dst_start_depth,
   int src_start_depth,
   int num_depth)
{
  const SAMRAI::tbox::Dimension& dim(dst.getDim());

  TYPE * const dst_ptr = dst.getPointer();
  const TYPE * const src_ptr = src.getPointer();

  const SAMRAI::hier::Box& dst_box(dst.getBox());
  const SAMRAI::hier::Box& src_box(src.getBox());

  int box_w[MAX_DIMENSION];
  int dst_w[MAX_DIMENSION];
  int src_w[MAX_DIMENSION];
  int dim_counter[MAX_DIMENSION];
  for (int i = 0; i < dim.getValue(); ++i) {
    box_w[i] = opbox.numberCells(i);
    dst_w[i] = dst_box.numberCells(i);
    src_w[i] = src_box.numberCells(i);
    dim_counter[i] = 0;
  }

  const int dst_offset = dst.getOffset();
  const int src_offset = src.getOffset();

  /*
   * Data on the opbox can be decomposed into a set of
   * contiguous array sections representing data in a straight line
   * in the 0 coordinate direction.
   *
   * num_d0_blocks is the number of such array sections.
   * dst_begin, src_begin are the array indices for the first
   * data items in each array section to be copied.
   */

  const int num_d0_blocks = opbox.size() / box_w[0];

  int dst_begin = dst_box.offset(opbox.lower())
    + dst_start_depth * dst_offset;
  int src_begin = src_box.offset(opbox.lower() - src_shift)
    + src_start_depth * src_offset;

  /*
   * Loop over the depth sections of the data arrays.
   */

  for (int d = 0; d < num_depth; ++d) {

    int dst_counter = dst_begin;
    int src_counter = src_begin;

    int dst_b[MAX_DIMENSION];
    int src_b[MAX_DIMENSION];
    for (int nd = 0; nd < dim.getValue(); ++nd) {
      dst_b[nd] = dst_counter;
      src_b[nd] = src_counter;
    }

    /*
     * Loop over each contiguous block of data.
     */

    for (int nb = 0; nb < num_d0_blocks; ++nb) {

      std::copy(&src_ptr[src_counter], &src_ptr[src_counter+box_w[0]], &dst_ptr[dst_counter]);

      int dim_jump = 0;

      /*
       * After each contiguous block is copied, calculate the
       * beginning array index for the next block.
       */

      for (int j = 1; j < dim.getValue(); ++j) {
        if (dim_counter[j] < box_w[j] - 1) {
          ++dim_counter[j];
          dim_jump = j;
          break;
        } else {
          dim_counter[j] = 0;
        }
      }

      if (dim_jump > 0) {

        int dst_step = 1;
        int src_step = 1;
        for (int k = 0; k < dim_jump; ++k) {
          dst_step *= dst_w[k];
          src_step *= src_w[k];
        }
        dst_counter = dst_b[dim_jump - 1] + dst_step;
        src_counter = src_b[dim_jump - 1] + src_step;

        for (int m = 0; m < dim_jump; ++m) {
          dst_b[m] = dst_counter;
          src_b[m] = src_counter;
        }

      }  // if dim_jump > 0

    }  // nb loop over contiguous data blocks

    /*
     * After copy is complete on a full box for one depth index,
     * advance by the offset values.
     */

    dst_begin += dst_offset;
    src_begin += src_offset;

  }  // d loop over depth indices
}

/*
 *************************************************************************
 *
 * Function that performs specified operation involving source and
 * destination data pointers and puts result in destination array
 * data object using explicit dimension-generic looping constructs.
 *
 *************************************************************************
 */

template<typename TYPE, int ALIGNMENT>
template<bool src_is_buffer>
void AlignedArrayDataOperationUtilities<TYPE, ALIGNMENT>::doArrayDataBufferOperationOnBox(
   const AlignedArrayData<TYPE, ALIGNMENT>& arraydata,
   const TYPE* buffer,
   const SAMRAI::hier::Box& opbox)
{
  const SAMRAI::tbox::Dimension& dim(arraydata.getDim());

  TYPE * const dst_ptr =
    (src_is_buffer ? const_cast<TYPE *>(arraydata.getPointer())
     : const_cast<TYPE *>(buffer));
  const TYPE * const src_ptr =
    (src_is_buffer ? buffer : arraydata.getPointer());

  const SAMRAI::hier::Box& array_d_box(arraydata.getBox());
  const int array_d_depth = arraydata.getDepth();

  int box_w[MAX_DIMENSION];
  int dat_w[MAX_DIMENSION];
  int dim_counter[MAX_DIMENSION];
  for (int i = 0; i < dim.getValue(); ++i) {
    box_w[i] = opbox.numberCells(i);
    dat_w[i] = array_d_box.numberCells(i);
    dim_counter[i] = 0;
  }

  const int dat_offset = arraydata.getOffset();
  const int buf_offset = box_w[0];

  /*
   * Data on the opbox can be decomposed into a set of
   * contiguous array sections representing data in a straight line
   * in the 0 coordinate direction.
   *
   * num_d0_blocks is the number of such array sections.
   * dat_begin, buf_begin are the array indices for the first
   * data items in each array section to be copied.
   */

  const int num_d0_blocks = opbox.size() / box_w[0];

  int dat_begin = array_d_box.offset(opbox.lower());
  int buf_begin = 0;

  /*
   * Loop over the depth sections of the data arrays.
   */

  for (int d = 0; d < array_d_depth; ++d) {

    int dat_counter = dat_begin;
    int buf_counter = buf_begin;

    int& dst_counter = (src_is_buffer ? dat_counter : buf_counter);
    int& src_counter = (src_is_buffer ? buf_counter : dat_counter);

    int dat_b[MAX_DIMENSION];
    for (int nd = 0; nd < dim.getValue(); ++nd) {
      dat_b[nd] = dat_counter;
    }

    /*
     * Loop over each contiguous block of data.
     */

    std::copy(&src_ptr[src_counter],
        &src_ptr[src_counter+box_w[0]],
        &dst_ptr[dst_counter]);

    int dim_jump = 0;

    /*
     * After each contiguous block is packed, calculate the
     * beginning array index for the next block.
     */

    for (int j = 1; j < dim.getValue(); ++j) {
      if (dim_counter[j] < box_w[j] - 1) {
        ++dim_counter[j];
        dim_jump = j;
        break;
      } else {
        dim_counter[j] = 0;
      }
    }

    if (dim_jump > 0) {

      int dat_step = 1;
      for (int k = 0; k < dim_jump; ++k) {
        dat_step *= dat_w[k];
      }
      dat_counter = dat_b[dim_jump - 1] + dat_step;

      for (int m = 0; m < dim_jump; ++m) {
        dat_b[m] = dat_counter;
      }

    }  // if dim_jump > 0

    buf_counter += buf_offset;

    /*
     * After packing is complete on a full box for one depth index,
     * advance by the offset value.
     */

    dat_begin += dat_offset;
    buf_begin = buf_counter;

  }  // d loop over depth indices
}

}
}

#endif // CLEVERLEAF_PDAT_ALIGNEDARRAYDATAOPERATIONUTILITIES_C_

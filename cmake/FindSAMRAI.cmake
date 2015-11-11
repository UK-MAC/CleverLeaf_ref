#.rst:
# FindSAMRAI
# ---------
#
# Find SAMRAI library
#
# This module finds an installed version of the SAMRAI library (see
# https://computation-rnd.llnl.gov/SAMRAI/). The list of libraries searched for
# includes all libraries built as part of the normal SAMRAI install.
#
# This module sets the following variables:
#
# ::
#
#   SAMRAI_FOUND - set to true if the SAMRAI library is found
#   SAMRAI_INCLUDE_DIRS - directory containing the SAMRAI header files
#   SAMRAI_ALGS_LIBRARIES - path of the SAMRAI_algs library
#   SAMRAI_APPU_LIBRARIES - path of the SAMRAI_appu library
#   SAMRAI_GEOM_LIBRARIES - path of the SAMRAI_geom library
#   SAMRAI_HIER_LIBRARIES - path of the SAMRAI_hier library
#   SAMRAI_MATH_LIBRARIES - path of the SAMRAI_math library
#   SAMRAI_MESH_LIBRARIES - path of the SAMRAI_mesh library
#   SAMRAI_PDAT_LIBRARIES - path of the SAMRAI_pdat library
#   SAMRAI_SOLV_LIBRARIES - path of the SAMRAI_solv library
#   SAMRAI_TBOX_LIBRARIES - path of the SAMRAI_tbox library
#   SAMRAI_XFER_LIBRARIES - path of the SAMRAI_xfer library

include (FindPackageHandleStandardArgs)

find_path (SAMRAI_PREFIX
  NAMES include/SAMRAI/SAMRAI_config.h
)

find_path (SAMRAI_INCLUDE_DIRS
  NAMES SAMRAI/SAMRAI_config.h 
  HINTS ${SAMRAI_PREFIX}/include
)

find_library (SAMRAI_ALGS_LIBRARIES
  NAMES libSAMRAI_algs.a
  HINTS ${SAMRAI_PREFIX}/lib
)

find_library (SAMRAI_APPU_LIBRARIES
  NAMES libSAMRAI_appu.a
  HINTS ${SAMRAI_PREFIX}/lib
)

find_library (SAMRAI_GEOM_LIBRARIES
  NAMES libSAMRAI_geom.a
  HINTS ${SAMRAI_PREFIX}/lib
)

find_library (SAMRAI_HIER_LIBRARIES
  NAMES libSAMRAI_hier.a
  HINTS ${SAMRAI_PREFIX}/lib
)

find_library (SAMRAI_MATH_LIBRARIES
  NAMES libSAMRAI_math.a
  HINTS ${SAMRAI_PREFIX}/lib
)

find_library (SAMRAI_MESH_LIBRARIES
  NAMES libSAMRAI_mesh.a
  HINTS ${SAMRAI_PREFIX}/lib
)

find_library (SAMRAI_PDAT_LIBRARIES
  NAMES libSAMRAI_pdat.a
  HINTS ${SAMRAI_PREFIX}/lib
)

find_library (SAMRAI_SOLV_LIBRARIES
  NAMES libSAMRAI_solv.a
  HINTS ${SAMRAI_PREFIX}/lib
)

find_library (SAMRAI_TBOX_LIBRARIES
  NAMES libSAMRAI_tbox.a
  HINTS ${SAMRAI_PREFIX}/lib
)

find_library (SAMRAI_XFER_LIBRARIES
  NAMES libSAMRAI_xfer.a
  HINTS ${SAMRAI_PREFIX}/lib
)

find_package_handle_standard_args (SAMRAI
  DEFAULT_MSG
  SAMRAI_INCLUDE_DIRS
  SAMRAI_ALGS_LIBRARIES
  SAMRAI_APPU_LIBRARIES
  SAMRAI_GEOM_LIBRARIES
  SAMRAI_HIER_LIBRARIES
  SAMRAI_MATH_LIBRARIES
  SAMRAI_MESH_LIBRARIES
  SAMRAI_PDAT_LIBRARIES
  SAMRAI_SOLV_LIBRARIES
  SAMRAI_TBOX_LIBRARIES
  SAMRAI_XFER_LIBRARIES
)

mark_as_advanced (
  SAMRAI_INCLUDE_DIRS
  SAMRAI_ALGS_LIBRARIES
  SAMRAI_APPU_LIBRARIES
  SAMRAI_GEOM_LIBRARIES
  SAMRAI_HIER_LIBRARIES
  SAMRAI_MATH_LIBRARIES
  SAMRAI_MESH_LIBRARIES
  SAMRAI_PDAT_LIBRARIES
  SAMRAI_SOLV_LIBRARIES
  SAMRAI_TBOX_LIBRARIES
  SAMRAI_XFER_LIBRARIES
)

set (SAMRAI_FOUND TRUE)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(geom)
add_subdirectory(hydro)
add_subdirectory(pdat)

add_executable(cleverleaf
  main.C)

target_link_libraries(cleverleaf
  hydro
  geom
  pdat
  ${SAMRAI_ALGS_LIBRARIES}
  ${SAMRAI_XFER_LIBRARIES}
  ${SAMRAI_APPU_LIBRARIES}
  ${SAMRAI_MATH_LIBRARIES}
  ${SAMRAI_MESH_LIBRARIES}
  ${SAMRAI_GEOM_LIBRARIES}
  ${SAMRAI_PDAT_LIBRARIES}
  ${SAMRAI_SOLV_LIBRARIES}
  ${SAMRAI_HIER_LIBRARIES}
  ${SAMRAI_TBOX_LIBRARIES}
  ${HDF5_LIBRARIES})

install(TARGETS cleverleaf
  RUNTIME DESTINATION bin)

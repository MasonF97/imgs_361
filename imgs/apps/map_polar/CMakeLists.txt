rit_add_executable(map_polar
  SOURCES
    map_polar.cpp
)

target_link_libraries(map_polar
  rit::ipcv_geometric_transformation 
  Boost::filesystem
  Boost::program_options
  opencv_core
  opencv_highgui
  opencv_imgcodecs
)

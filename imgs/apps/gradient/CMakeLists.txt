rit_add_executable(gradient 
  SOURCES
    gradient.cpp
)

target_link_libraries(gradient
  rit::ipcv_spatial_filtering
  Boost::filesystem
  Boost::program_options
  opencv_core
  opencv_highgui
  opencv_imgcodecs
)

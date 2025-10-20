rit_add_library(ipcv_spatial_filtering
  SOURCES
    Filter2D.cpp
    Sobel.cpp
  HEADERS
    Filter2D.h
    Sobel.h
)

target_link_libraries(ipcv_spatial_filtering 
  PUBLIC 
    opencv_core
  PRIVATE
)

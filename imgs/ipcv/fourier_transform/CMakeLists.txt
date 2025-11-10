rit_add_library(ipcv_fourier_transform
  SOURCES
    fourier_transform_display.cpp
  HEADERS
    fourier_transform_display.h
)

target_link_libraries(ipcv_fourier_transform 
  PUBLIC 
    opencv_core
  PRIVATE
    rit::ipcv_utils

)

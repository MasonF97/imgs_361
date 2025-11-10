rit_add_executable(fft_display
  SOURCES
    fft_display.cpp
)

target_link_libraries(fft_display
  rit::ipcv_fourier_transform 
  Boost::filesystem
  Boost::program_options
  opencv_core
  opencv_highgui
  opencv_imgcodecs
)

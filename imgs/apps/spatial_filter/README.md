# Uploaded Files
These are all of the uploaded files

## 1
File name in submission:
spatial_filter.cpp 

Destination path:
./imgs/apps/spatial_filter/spatial_filter.cpp 

## 2
File name in submission:
Filter2D.cpp

Destination path:
./imgs/ipcv/spatial_filtering/Filter2D.cpp

## 3
File name in submission:
gradient.cpp

Destination path:
./imgs/apps/gradient/gradient.cpp

## 4
File name in submission:
Sobel.cpp

Destination path:
./imgs/ipcv/spatial_filtering/Sobel.cpp

## 5
File name in submission:
Sobel.h

Destination path:
./imgs/ipcv/spatial_filtering/Sobel.h


# Files to Update
These are all of the files that need to be updated. Sorry, I could only upload 6 files.

## 1
File: "./CMakeLists.txt"

Add this line after line 32(Not sure if this is needed, it may just be me, but I didn't want to risk it):
```
include_directories(${EIGEN3_INCLUDE_DIR})
```

## 2
File: "./imgs/ipcv/spatial_filtering/Filter2D.h"
Add "WRAP" and "ISOLATED" to the border modes class and set the default border mode to "WRAP"

## 3
File: "./imgs/ipcv/spatial_filtering/CMakeLists.txt"
Replace with this:
```
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
```

## 4
Create this folder: "./imgs/apps/gradient"

## 5
Create File: "./imgs/apps/gradient/CMakeLists.txt"
With the following:
```
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
```

## 6
File: "./imgs/apps/CMakeLists.txt"
Add the following:
```
add_subdirectory(gradient)
```


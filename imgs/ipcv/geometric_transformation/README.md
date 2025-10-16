# Files
I could only upload 6 files, so some require manual editing, sorry.

## 1
File name in submission:
CMakeLists.txt

Destination path: 
./CMakeLists.txt

## 2
File name in submission:
MapQ2Q.cpp

Destination path: 
./imgs/ipcv/geometric_transformation/MapQ2Q.cpp

## 3
File name in submission:
MapPolar.cpp

Destination path: 
./imgs/ipcv/geometric_transformation/MapPolar.cpp

## 4
File name in submission:
MapPolar.h

Destination path: 
./imgs/ipcv/geometric_transformation/MapPolar.h

## 5
Update this file: 
./imgs/ipcv/geometric_transformation/CMakeLists.txt
```
rit_add_library(ipcv_geometric_transformation
  SOURCES
    MapGCP.cpp
    MapPolar.cpp
    MapQ2Q.cpp
    MapRST.cpp
    Remap.cpp
  HEADERS
    MapGCP.h
    MapPolar.h
    MapQ2Q.h
    MapRST.h
    Remap.h
    GeometricTransformation.h
)

target_link_libraries(ipcv_geometric_transformation 
  PUBLIC 
    Eigen3::Eigen
    opencv_core
  PRIVATE
)
```

## 6
Add the following to this file: 
./imgs/ipcv/geometric_transformation/GeometricTransformation.h
```
#include "imgs/ipcv/geometric_transformation/MapPolar.h"
```

## 7
Add the following to this file:: 
./imgs/apps/CMakeLists.txt
```
add_subdirectory(map_polar)
```

## 8
Create this folder:
./imgs/apps/map_polar

## 9
Update this file: 
./imgs/apps/map_polar/CMakeLists.txt
```
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

```

## 10
File name in submission:
imgs.apps.map_polar.map_polar.cpp

Destination path: 
./imgs/apps/map_polar/map_polar.cpp


# To Run
Navigate to the build directory and run:
```
bin/map_polar -i <path-to-image>
```
or to do log polar:
```
bin/map_polar -i <path-to-image> -L
```
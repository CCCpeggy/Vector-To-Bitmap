cmake_minimum_required(VERSION 2.8)

project(VectorToBitmapFLTK)
set(SRC_DIR ${PROJECT_SOURCE_DIR}/src)
set(INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(LIB_DIR ${PROJECT_SOURCE_DIR}/lib)
set(CVSystem_DIR ${SRC_DIR}/CVSystem)

include_directories(${INCLUDE_DIR})
include_directories(${INCLUDE_DIR}/glad4.6/include)
include_directories(${INCLUDE_DIR}/glm-0.9.8.5/glm)

add_Definitions("-D_XKEYCHECK_H")

add_executable(VectorToBitmapFLTK
    ${SRC_DIR}/CallBacks.h
    ${SRC_DIR}/CallBacks.cpp
    ${SRC_DIR}/ControlPoint.h
    ${SRC_DIR}/ControlPoint.cpp
    ${SRC_DIR}/main.cpp
    ${SRC_DIR}/Object.h
    ${SRC_DIR}/Track.h
    ${SRC_DIR}/Track.cpp
    ${SRC_DIR}/PreView.h
    ${SRC_DIR}/PreView.cpp
    ${SRC_DIR}/TextureView.h
    ${SRC_DIR}/TextureView.cpp
    ${SRC_DIR}/TrainWindow.h
    ${SRC_DIR}/TrainWindow.cpp
    ${SRC_DIR}/Shader/ShaderObject.cpp
    ${SRC_DIR}/Shader/BasicTextureShader.cpp
    ${SRC_DIR}/Shader/ScreentoneShader.cpp
    ${SRC_DIR}/Common.cpp
    ${INCLUDE_DIR}/Shader/ShaderObject.h
    ${INCLUDE_DIR}/Shader/BasicTextureShader.h
    ${INCLUDE_DIR}/Shader/ScreentoneShader.h
    ${INCLUDE_DIR}/Common.h
    ${INCLUDE_DIR}/MyObj.h
    ${INCLUDE_DIR}/glad4.6/src/glad.c
    ${CVSystem_DIR}/ScreentoneSegmentation.cpp
    ${CVSystem_DIR}/SystemParams.cpp
    ${CVSystem_DIR}/Triangulator1.cpp
    ${CVSystem_DIR}/PixelsInTriangle.cpp
    ${CVSystem_DIR}/PointInTrianglesTest.cpp
    ${CVSystem_DIR}/CurveInterpolation.cpp
    ${CVSystem_DIR}/CurveFitting.cpp 
    ${CVSystem_DIR}/UtilityFunctions.cpp
    ${CVSystem_DIR}/CSSSmoothing.cpp
    ${CVSystem_DIR}/CurveRDP.cpp
    ${CVSystem_DIR}/LinesSorter.cpp
    ${CVSystem_DIR}/CGALTools.cpp)

add_library(Utilities 
    ${SRC_DIR}/Utilities/3DUtils.h
    ${SRC_DIR}/Utilities/3DUtils.cpp
    ${SRC_DIR}/Utilities/ArcBallCam.h
    ${SRC_DIR}/Utilities/ArcBallCam.cpp
    ${SRC_DIR}/Utilities/Pnt3f.h
    ${SRC_DIR}/Utilities/Pnt3f.cpp)

target_link_libraries(VectorToBitmapFLTK debug 
    ${LIB_DIR}/Debug/fltk_formsd.lib 
    ${LIB_DIR}/Debug/fltk_gld.lib
    ${LIB_DIR}/Debug/fltk_imagesd.lib
    ${LIB_DIR}/Debug/fltk_jpegd.lib
    ${LIB_DIR}/Debug/fltk_pngd.lib
    ${LIB_DIR}/Debug/fltk_zd.lib
    ${LIB_DIR}/Debug/fltkd.lib)

target_link_libraries(VectorToBitmapFLTK optimized 
    ${LIB_DIR}/Release/fltk_forms.lib
    ${LIB_DIR}/Release/fltk_gl.lib
    ${LIB_DIR}/Release/fltk_images.lib
    ${LIB_DIR}/Release/fltk_jpeg.lib
    ${LIB_DIR}/Release/fltk_png.lib
    ${LIB_DIR}/Release/fltk_z.lib
    ${LIB_DIR}/Release/fltk.lib)

target_link_libraries(VectorToBitmapFLTK 
    ${LIB_DIR}/OpenGL32.lib
    ${LIB_DIR}/glu32.lib)

target_link_libraries(VectorToBitmapFLTK Utilities)

# Potrace
set(Potrace_INCLUDE_DIR D:/VG_Library/Potrace/Potrace)
set(Potrace_LIB_DIR D:/VG_Library/Potrace/x64)

include_directories(${Potrace_INCLUDE_DIR})
target_link_libraries(VectorToBitmapFLTK debug 
    ${Potrace_LIB_DIR}/Debug/Potrace_d.lib)

target_link_libraries(VectorToBitmapFLTK optimized 
    ${Potrace_LIB_DIR}/Release/Potrace.lib)

# CGAL
set(CGAL_INCLUDE_DIRS
    D:/VG_Library/CGAL-4.14.3/include
    D:/VG_Library/CGAL-4.14.3/VC16/include
    D:/VG_Library/CGAL-4.14.3/auxiliary/gmp/include)
set(CGAL_GMP_LIB_DIRS D:/VG_Library/CGAL-4.14.3/auxiliary/gmp/lib)
set(CGAL_LIB_DIR  D:/VG_Library/CGAL-4.14.3/VC16/lib)

include_directories(${CGAL_INCLUDE_DIRS})
target_link_libraries(VectorToBitmapFLTK ${CGAL_GMP_LIB_DIRS}/libgmp-10.lib)
target_link_libraries(VectorToBitmapFLTK ${CGAL_GMP_LIB_DIRS}/libmpfr-4.lib)
link_directories(${CGAL_GMP_LIB_DIRS})
link_directories(${CGAL_LIB_DIR})
link_directories(${CGAL_LIB_DIR}/Debug)
link_directories(${CGAL_LIB_DIR}/Release)

# OpenCV
set(OpenCV_DIR D:/VG_Library/opencv440/build/x64/vc15/lib)

find_package( OpenCV REQUIRED )
message( STATUS "OpenCV library status:" )
message( STATUS "    config: ${OpenCV_DIR}" )
message( STATUS "    version: ${OpenCV_VERSION}" )
message( STATUS "    libraries: ${OpenCV_LIBS}" )
message( STATUS "    include path: ${OpenCV_INCLUDE_DIRS}" )

include_directories(${OpenCV_INCLUDE_DIRS})
target_link_libraries(VectorToBitmapFLTK ${OpenCV_LIBS})

# Boost
SET(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "D:/VG_Library/boost_1_74_0")
SET(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "D:/VG_Library/boost_1_74_0/libs")
SET(BOOST_ROOT "D:/G_Library/boost_1_74_0")
set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 
find_package(Boost 1.48.0)

include_directories(${Boost_INCLUDE_DIRS}) 
target_link_libraries(VectorToBitmapFLTK D:/VG_Library/boost_1_74_0/stage/lib/libboost_chrono-vc142-mt-gd-x64-1_74.lib)

# Eigen
SET(Eigen_INCLUDE_DIRS "D:/VG_Library/Eigen-3.3.7")

include_directories(${Eigen_INCLUDE_DIRS})

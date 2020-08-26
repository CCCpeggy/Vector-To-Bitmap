cmake_minimum_required(VERSION 2.8)

project(VectorToBitmapFLTK)
set(SRC_DIR ${PROJECT_SOURCE_DIR}/src)
set(INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(LIB_DIR ${PROJECT_SOURCE_DIR}/lib)
set(OpenCV_DIR D:/VG_Library/opencv440/build/x64/vc15/lib)

find_package( OpenCV REQUIRED )
message( STATUS "OpenCV library status:" )
message( STATUS "    config: ${OpenCV_DIR}" )
message( STATUS "    version: ${OpenCV_VERSION}" )
message( STATUS "    libraries: ${OpenCV_LIBS}" )
message( STATUS "    include path: ${OpenCV_INCLUDE_DIRS}" )

include_directories(${INCLUDE_DIR})
include_directories(${INCLUDE_DIR}/glad4.6/include)
include_directories(${INCLUDE_DIR}/glm-0.9.8.5/glm)
include_directories(${OpenCV_INCLUDE_DIRS})

add_Definitions("-D_XKEYCHECK_H")

# find required opencv
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
    ${INCLUDE_DIR}/glad4.6/src/glad.c)

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
    ${LIB_DIR}/Debug/fltkd.lib
    ${OpenCV_LIBS})

target_link_libraries(VectorToBitmapFLTK optimized 
    ${LIB_DIR}/Release/fltk_forms.lib
    ${LIB_DIR}/Release/fltk_gl.lib
    ${LIB_DIR}/Release/fltk_images.lib
    ${LIB_DIR}/Release/fltk_jpeg.lib
    ${LIB_DIR}/Release/fltk_png.lib
    ${LIB_DIR}/Release/fltk_z.lib
    ${LIB_DIR}/Release/fltk.lib
    ${OpenCV_LIBS})

target_link_libraries(VectorToBitmapFLTK 
    ${LIB_DIR}/OpenGL32.lib
    ${LIB_DIR}/glu32.lib)

target_link_libraries(VectorToBitmapFLTK Utilities)
    
#ifndef COMMON_H
#define COMMON_H

// define
#define M_EPS 1e-8

// stl
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <fstream> 
#include <list>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <time.h>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <cassert>

// opengl
#include <windows.h> // OpenGL needs windows.h
// #include "GL/gl.h"
#include <glad/glad.h>
#include <glm/glm.hpp>
#include "GL/glu.h"

// opencv
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/types_c.h>

// FLTK
#include <Fl/fl.h>
#include <Fl/Fl_Double_Window.h>
#include <Fl/Fl_Button.h>
#include <Fl/Fl_Group.H>
#include <Fl/Fl_Value_Slider.H>
#include <Fl/Fl_Browser.H>
#include <Fl/Fl_Gl_Window.h>
#include <FL/Fl_Box.h>
#include <FL/Fl_Input.h>
#include <Fl/Fl_File_Chooser.H>
#include <Fl/math.h>

// Boost
#include <boost/chrono.hpp>

// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Spatial_sort_traits_adapter_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Segment_2.h>

#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>

#include <CGAL/Triangulation_conformer_2.h>

// Potrace
extern "C"
{
#include "potracelib.h"
#include "platform.h"
#include "bitmap_io.h"
#include "potrace_main.h"
}

// CVSystem
#include "CVSystem/TriangleType.h"
#include "CVSystem/MyPoint.h"
#include "CVSystem/MyIndexedBezierCurves.h"
#include "CVSystem/MyIndexedTriangle.h"
#include "CVSystem/MyQuad.h"
#include "CVSystem/MyLine.h"
#include "CVSystem/MyTriangle.h"
#include "CVSystem/MyIndexedLine.h"
#include "CVSystem/TriangleType.h"
#include "CVSystem/ScreentoneSegmentation.h"
#include "CVSystem/SystemParams.h"
#include "CVSystem/CSSSmoothing.h"
#include "CVSystem/Triangulator1.h"
#include "CVSystem/Delaunay_mesh_face_base_info_2.h"
#include "CVSystem/CGALTools.h"
#include "CVSystem/PixelsInTriangle.h"
#include "CVSystem/PointInTrianglesTest.h"
#include "CVSystem/CurveInterpolation.h"
#include "CVSystem/CurveFitting.h"
#include "CVSystem/UtilityFunctions.h"
#include "CVSystem/CurveRDP.h"
#include "CVSystem/LineCloud.h"
#include "CVSystem/LinesSorter.h"
#include "CVSystem/Triangulator1.h"
#include "CVSystem/NanoFLANN.h"

// Utilities
#include "../src/Utilities/3DUtils.H"
#include "../src/Utilities/ArcBallCam.H"

// Eigen
#include <Eigen/Dense>

// Path
class Path {
public:
	static std::string ShaderPath;
};


typedef struct _TextureData
{
	_TextureData() : width(0), height(0), data(0), channel(0), idx(0) {}
	int width;
	int height;
	int channel;
	unsigned char* data;
	GLuint idx;
} TextureData;

class Common
{
public:
	static void DumpInfo(void);
	static void ShaderLog(GLuint shader);
	static void PrintGLError();
	static TextureData LoadPng(const char* path);
	static GLuint LoadTexture(const char* path, int& width, int& height);
	static void LoadTexture(const char* path, TextureData& texData);
	static const char** LoadShaderSource(const char* file);
	static void FreeShaderSource(const char** srcp);
	static bool CheckShaderCompiled(GLuint shader);
	static bool CheckProgramLinked(GLuint program);
	static bool CheckFrameBufferStatus();
	static bool CheckGLError();
	static void SavePng(const char* path, GLuint tex, int width, int height);
	static bool FileExists(const char* filePath);

};

#endif
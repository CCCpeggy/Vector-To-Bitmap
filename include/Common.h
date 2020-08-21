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
#include <math.h>

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

// CVSystem
#include "../Include/CVSystem/TriangleType.h"
#include "../Include/CVSystem/MyPoint.h"
#include "../Include/CVSystem/MyIndexedBezierCurves.h"
#include "../Include/CVSystem/MyIndexedTriangle.h"
#include "../Include/CVSystem/MyQuad.h"
#include "../Include/CVSystem/TriangleType.h"


// Utilities
#include "../src/Utilities/3DUtils.H"
#include "../src/Utilities/ArcBallCam.H"
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

	//// String
	// split
	static std::vector<std::string>& split(const std::string& s, char delim, std::vector<std::string>& elems);
	// split
	static std::vector<std::string> split(const std::string& s, char delim);
};

#endif
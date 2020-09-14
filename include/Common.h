#ifndef COMMON_H
#define COMMON_H

// define
#define M_EPS 1e-8

// stl
#include "stllibrary.h"

// opengl
#include "opengl.h"

// opencv
#include "opencv.h"

// FLTK
#include "fltk.h"

// Boost
#include <boost/chrono.hpp>

// CGAL
#include "cgal.h"

// Potrace
#include "potrace.h"

// CVSystem
#include "cvsystem.h"

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
#include "../Include/Common.h"

std::string Path::ShaderPath = "../shader/";

// STB
#define STB_IMAGE_IMPLEMENTATION
#include "../include/STB/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../include/STB/stb_image_write.h"

void Common::DumpInfo(void)
{
	printf("Vendor: %s\n", glGetString(GL_VENDOR));
	printf("Renderer: %s\n", glGetString(GL_RENDERER));
	printf("Version: %s\n", glGetString(GL_VERSION));
	printf("GLSL: %s\n", glGetString(GL_SHADING_LANGUAGE_VERSION));
}

void Common::ShaderLog(GLuint shader)
{
	GLint isCompiled = 0;
	glGetShaderiv(shader, GL_COMPILE_STATUS, &isCompiled);
	if (isCompiled == GL_FALSE)
	{
		GLint maxLength = 0;
		glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &maxLength);

		// The maxLength includes the NULL character
		GLchar* errorLog = new GLchar[maxLength];
		glGetShaderInfoLog(shader, maxLength, &maxLength, &errorLog[0]);

		printf("%s\n", errorLog);
		delete[] errorLog;
	}
}

void Common::PrintGLError()
{
	GLenum code = glGetError();
	switch (code)
	{
	case GL_NO_ERROR:
		std::cout << "GL_NO_ERROR" << std::endl;
		break;
	case GL_INVALID_ENUM:
		std::cout << "GL_INVALID_ENUM" << std::endl;
		break;
	case GL_INVALID_VALUE:
		std::cout << "GL_INVALID_VALUE" << std::endl;
		break;
	case GL_INVALID_OPERATION:
		std::cout << "GL_INVALID_OPERATION" << std::endl;
		break;
	case GL_INVALID_FRAMEBUFFER_OPERATION:
		std::cout << "GL_INVALID_FRAMEBUFFER_OPERATION" << std::endl;
		break;
	case GL_OUT_OF_MEMORY:
		std::cout << "GL_OUT_OF_MEMORY" << std::endl;
		break;
	case GL_STACK_UNDERFLOW:
		std::cout << "GL_STACK_UNDERFLOW" << std::endl;
		break;
	case GL_STACK_OVERFLOW:
		std::cout << "GL_STACK_OVERFLOW" << std::endl;
		break;
	default:
		std::cout << "GL_ERROR" << std::endl;
	}
}

bool Common::CheckShaderCompiled(GLuint shader)
{
	GLint isCompiled = 0;
	glGetShaderiv(shader, GL_COMPILE_STATUS, &isCompiled);
	if (isCompiled == GL_FALSE)
	{
		GLint maxLength = 0;
		glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &maxLength);

		// The maxLength includes the NULL character
		GLchar* errorLog = new GLchar[maxLength];
		glGetShaderInfoLog(shader, maxLength, &maxLength, &errorLog[0]);

		printf("%s\n", errorLog);
		delete[] errorLog;

		return false;
	}

	return true;
}

bool Common::CheckProgramLinked(GLuint program)
{
	GLint isLinked = 0;
	glGetProgramiv(program, GL_LINK_STATUS, &isLinked);
	if (isLinked == GL_FALSE) {

		GLint maxLength = 0;
		glGetProgramiv(program, GL_INFO_LOG_LENGTH, &maxLength);

		GLchar* errorLog = new GLchar[maxLength];
		glGetProgramInfoLog(program, sizeof(errorLog), NULL, errorLog);
		fprintf(stderr, "Program link error: %s\n", errorLog);

		return false;
	}
	return true;
}

bool Common::CheckFrameBufferStatus()
{
	GLenum fboStatus = glCheckFramebufferStatus(GL_FRAMEBUFFER);

	if (fboStatus != GL_FRAMEBUFFER_COMPLETE) {
		printf("FBO error: %d\n", fboStatus);
		return false;
	}

	return true;
}

bool Common::CheckGLError()
{
	GLenum errCode = glGetError();
	if (errCode != GL_NO_ERROR)
	{
		const GLubyte* errString = gluErrorString(errCode);
		printf("%s\n", errString);

		return false;
	}
	return true;
}

GLuint Common::LoadTexture(const char* path, int& width, int& height)
{
	TextureData data = LoadPng(path);

	width = data.width;
	height = data.height;

	GLuint texture;
	glGenTextures(1, &texture);
	glBindTexture(GL_TEXTURE_2D, texture);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, data.data);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glBindTexture(GL_TEXTURE_2D, 0);
	stbi_image_free(data.data);
	return texture;
}

void Common::LoadTexture(const char* path, TextureData& texData)
{
	texData = LoadPng(path);

	glGenTextures(1, &texData.idx);
	glBindTexture(GL_TEXTURE_2D, texData.idx);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, texData.width, texData.height, 0, GL_RGBA, GL_UNSIGNED_BYTE, texData.data);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glBindTexture(GL_TEXTURE_2D, 0);
	stbi_image_free(texData.data);
}

void Common::SavePng(const char* path, GLuint tex, int width, int height)
{
	int size = width * height * 4;
	uint8_t* data = new uint8_t[size];
	glBindTexture(GL_TEXTURE_2D, tex);
	glGetTexImage(GL_TEXTURE_2D, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);
	stbi_write_png(path, width, height, 4, data, width * 4);
	glBindTexture(GL_TEXTURE_2D, 0);
	delete[] data;
}


//Read shader file
const char** Common::LoadShaderSource(const char* file)
{
	FILE* fp = fopen(file, "rb");
	fseek(fp, 0, SEEK_END);
	long sz = ftell(fp);
	fseek(fp, 0, SEEK_SET);
	char* src = new char[sz + 1];
	fread(src, sizeof(char), sz, fp);
	src[sz] = '\0';
	const char** srcp = new const char* [1];
	srcp[0] = src;
	return srcp;
}

//Release 2-dimension array
void Common::FreeShaderSource(const char** srcp)
{
	delete srcp[0];
	delete srcp;
}

TextureData Common::LoadPng(const char* path)
{
	TextureData texture;
	texture.data = stbi_load(path, &texture.width, &texture.height, &texture.channel, 4);
	return texture;
}

bool Common::FileExists(const char* filePath)
{
	struct stat stFileInfo;
	int intStat;
	intStat = stat(filePath, &stFileInfo);
	if (intStat == 0) {
		return true;
	}
	else {
		return false;
	}
}

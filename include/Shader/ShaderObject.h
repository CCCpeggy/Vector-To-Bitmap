#ifndef __SHADEROBJECT__
#define __SHADEROBJECT__

#include "../Common.h"

class ShaderObject
{
public:
	ShaderObject();
	~ShaderObject();

	virtual bool Init();
	virtual bool AddShader(GLenum shaderType, std::string fileName);
	virtual bool Finalize();
	virtual void Enable();
	virtual void Disable();
	virtual GLuint GetUniformLocation(const GLchar* name);

private:
	std::list<GLuint> shaderList;
	GLuint program;
};

#endif

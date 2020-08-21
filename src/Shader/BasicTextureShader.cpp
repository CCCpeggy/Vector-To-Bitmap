#include "../../Include/Shader/BasicTextureShader.h"
#include "../../Include/Common.h"

BasicTextureShader::BasicTextureShader()
{
}


BasicTextureShader::~BasicTextureShader()
{
}

bool BasicTextureShader::Init()
{
	if (!ShaderObject::Init())
	{
		return false;
	}
	if (!AddShader(GL_VERTEX_SHADER, Path::ShaderPath + "BasicShader.vs.glsl"))
	{
		return false;
	}
	
	if (!AddShader(GL_FRAGMENT_SHADER, Path::ShaderPath + "BasicShader.fs.glsl"))
	{
		return false;
	}
	
	if (!Finalize())
	{
		return false;
	}

	return true;
}

void BasicTextureShader::SetTexture(GLuint tex)
{
	glBindTexture(GL_TEXTURE_2D, tex);
}

void BasicTextureShader::Enable()
{
	ShaderObject::Enable();
}

void BasicTextureShader::Disable()
{
	ShaderObject::Disable();
	glBindTexture(GL_TEXTURE_2D, 0);
}

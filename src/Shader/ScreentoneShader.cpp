#include "../../Include/Shader/ScreentoneShader.h"
#include "../../Include/Common.h"

ScreentoneShader::ScreentoneShader(){
}


ScreentoneShader::~ScreentoneShader()
{
}

bool ScreentoneShader::Init()
{
	if (!ShaderObject::Init())
	{
		return false;
	}

	if (!AddShader(GL_VERTEX_SHADER, Path::ShaderPath + "Screentone.vs.glsl"))
	{
		return false;
	}

	if (!AddShader(GL_FRAGMENT_SHADER, Path::ShaderPath + "Screentone.fs.glsl"))
	{
		return false;
	}
	
	if (!Finalize())
	{
		return false;
	}

	triangleTypeLocation = GetUniformLocation("triangleType");
	if (triangleTypeLocation == -1)
	{
		puts("Get uniform loaction error: triangleType");
		return false;
	}

	triangleTypeLocation = GetUniformLocation("triangleType");
	if (triangleTypeLocation == -1)
	{
		puts("Get uniform loaction error: triangleType");
		return false;
	}

	inImgWidthLocation = GetUniformLocation("inImgWidth");
	if (inImgWidthLocation == -1)
	{
		puts("Get uniform loaction error: inImgWidth");
		return false;
	}

	inImgHeightLocation = GetUniformLocation("inImgHeight");
	if (inImgHeightLocation == -1)
	{
		puts("Get uniform loaction error: inImgHeight");
		return false;
	}

	outImgWidthLocation = GetUniformLocation("outImgWidth");
	if (outImgWidthLocation == -1)
	{
		puts("Get uniform loaction error: outImgWidth");
		return false;
	}

	outImgHeightLocation = GetUniformLocation("outImgHeight");
	if (outImgHeightLocation == -1)
	{
		puts("Get uniform loaction error: outImgHeight");
		return false;
	}

	texWidthLocation = GetUniformLocation("texWidth");
	if (texWidthLocation == -1)
	{
		puts("Get uniform loaction error: texWidth");
		return false;
	}

	texHeightLocation = GetUniformLocation("texHeight");
	if (texHeightLocation == -1)
	{
		puts("Get uniform loaction error: texHeight");
		return false;
	}

	winWidthLocation = GetUniformLocation("windowWidth");
	if (texWidthLocation == -1)
	{
		puts("Get uniform loaction error: windowWidth");
		return false;
	}

	winHeightLocation = GetUniformLocation("windowHeight");
	if (texHeightLocation == -1)
	{
		puts("Get uniform loaction error: windowHeight");
		return false;
	}

	return true;
}

void ScreentoneShader::SetTriangleType(CVSystem::TriangleType type)
{
	glUniform1i(triangleTypeLocation, type);
}

void ScreentoneShader::SetTexture(GLuint tex)
{
	glBindTexture(GL_TEXTURE_2D, tex);
}

void ScreentoneShader::SetInImgSize(int width, int height)
{
	glUniform1i(inImgWidthLocation, width);
	glUniform1i(inImgHeightLocation, height);
}

void ScreentoneShader::SetOutImgSize(int width, int height)
{
	glUniform1i(outImgWidthLocation, width);
	glUniform1i(outImgHeightLocation, height);
}

void ScreentoneShader::SetOutImgWidth(int width)
{
	glUniform1i(outImgWidthLocation, width);
}

void ScreentoneShader::SetOutImgHeight(int height)
{
	glUniform1i(outImgHeightLocation, height);
}

void ScreentoneShader::SetTexSize(int width, int height)
{
	glUniform1i(texWidthLocation, width);
	glUniform1i(texHeightLocation, height);
}

void ScreentoneShader::SetWinSize(int width, int height)
{
	glUniform1i(winWidthLocation, width);
	glUniform1i(winHeightLocation, height);
}

void ScreentoneShader::Enable()
{
	ShaderObject::Enable();
}

void ScreentoneShader::Disable()
{
	ShaderObject::Disable();
	glBindTexture(GL_TEXTURE_2D, 0);
}

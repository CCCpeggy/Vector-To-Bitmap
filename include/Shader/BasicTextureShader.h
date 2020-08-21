#ifndef __BASICTEXSHADER__
#define __BASICTEXSHADER__
#include "../../Include/Shader/ShaderObject.h"

class BasicTextureShader: public ShaderObject
{
public:
	BasicTextureShader();
	~BasicTextureShader();

	bool Init();
	void SetTexture(GLuint tex);
	void Enable();
	void Disable();

private:
};

#endif
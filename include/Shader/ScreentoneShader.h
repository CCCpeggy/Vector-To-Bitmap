#ifndef __SCREENTONESHADER__
#define __SCREENTONESHADER__
#include "../../Include/Shader/ShaderObject.h"
#include "../../Include/CVSystem/TriangleType.h"

class ScreentoneShader : public ShaderObject
{
public:
	ScreentoneShader();
	~ScreentoneShader();

	bool Init();
	void SetTriangleType(CVSystem::TriangleType type);
	void SetTexture(GLuint tex);
	void SetInImgSize(int width, int height);
	void SetOutImgSize(int width, int height);
	void SetOutImgWidth(int width);
	void SetOutImgHeight(int height);
	void SetTexSize(int width, int height);
	void SetWinSize(int width, int height);
	void Enable();
	void Disable();

private:
	GLuint triangleTypeLocation;
	GLuint inImgWidthLocation;
	GLuint inImgHeightLocation;
	GLuint outImgWidthLocation;
	GLuint outImgHeightLocation;
	GLuint texWidthLocation;
	GLuint texHeightLocation;
	GLuint winWidthLocation;
	GLuint winHeightLocation;
};

#endif
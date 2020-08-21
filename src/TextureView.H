#include "../Include/Common.h"
#include "../Include/Shader/BasicTextureShader.h"
#include "../Include/Shader/ScreentoneShader.h"
#pragma once
class TrainWindow;

#pragma warning(push)
#pragma warning(disable:4312)
#pragma warning(disable:4311)
#pragma warning(pop)


class TextureView : public Fl_Gl_Window
{
	public:
		// note that we keep the "standard widget" constructor arguments
		TextureView(int x, int y, int w, int h, const char* l = 0);

		// overrides of important window things
		virtual int handle(int);
		virtual void draw();
		void openTexture(std::string);

	public:
		TrainWindow*		tw;				// The parent of this display window
		BasicTextureShader	basicTeTextureShader;
		TextureData			loadTex;
		bool				firstDraw;
};
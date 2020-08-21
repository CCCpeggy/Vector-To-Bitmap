#include "../Include/Common.h"
#include "../Include/Shader/BasicTextureShader.h"
#include "../Include/Shader/ScreentoneShader.h"
#include "../Include/MyObj.h"
#pragma once
class TrainWindow;
class CTrack;

#pragma warning(push)
#pragma warning(disable:4312)
#pragma warning(disable:4311)
#pragma warning(pop)


class PreView : public Fl_Gl_Window
{
	public:
		// note that we keep the "standard widget" constructor arguments
		PreView(int x, int y, int w, int h, const char* l = 0);

		// overrides of important window things
		virtual int handle(int);
		virtual void draw();
		virtual void render();
		void setImgOutputSize(int width, int height);
		void openMyObjFile(std::string filePath);

		// pick a point (for when the mouse goes down)
		void doPick();

	public:
		int					selectedCube;  // simple - just remember which cube is selected

		TrainWindow*		tw;				// The parent of this display window
		CTrack*				m_pTrack;		// The track of the entire scene
		TextureData			fbxTex;

		BasicTextureShader	basicTeTextureShader;
		ScreentoneShader	screenToneShader;
		GLuint				fbo;
		bool				firstDraw;
		MyObj				myObj;
		int					windowWidth;
		int					windowHeight;
};
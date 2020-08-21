#include "TextureView.H"
#include "TrainWindow.H"
#include "CallBacks.H"


TextureView::
TextureView(int x, int y, int w, int h, const char* l)
	: Fl_Gl_Window(x, y, w, h, l),
	firstDraw(true)
{
	mode( FL_RGB|FL_ALPHA|FL_DOUBLE | FL_STENCIL );
}

int TextureView::handle(int event)
{
	// remember what button was used
	static int last_push;

	switch (event) {
		// Mouse button being pushed event
	case FL_PUSH:
		Fl_File_Chooser* fc = new Fl_File_Chooser(".", "Image Files (*.{png})", 0, "Input File");
		fc->callback(openTextureFile, tw);
		fc->show();
		break;
	}
	return Fl_Gl_Window::handle(event);
}

void TextureView::draw()
{
	//initialized glad
	if (firstDraw) {
		if (!gladLoadGL()) {
			throw std::runtime_error("Could not initialize GLAD!");
			return;
		}
		firstDraw = false;

		Common::LoadTexture("../asset/in.png", loadTex);

		basicTeTextureShader.Init();
	}

	glClearColor(1, 1, 1, 0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glEnable(GL_TEXTURE_2D);
	basicTeTextureShader.Enable();
	basicTeTextureShader.SetTexture(loadTex.idx);
	glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
	basicTeTextureShader.Disable();
}

void TextureView::openTexture(std::string filePath)
{
	Common::LoadTexture(filePath.c_str(), loadTex);
}


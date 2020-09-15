#include "PreView.H"
#include "TrainWindow.H"
#include "TextureView.H"


PreView::
PreView(int x, int y, int w, int h, const char* l)
	: Fl_Gl_Window(x, y, w, h, l),
	firstDraw(true),
	windowWidth(w),
	windowHeight(h)
{
	mode( FL_RGB|FL_ALPHA|FL_DOUBLE | FL_STENCIL );
}

int PreView::handle(int event)
{
	// remember what button was used
	static int last_push;

	switch(event) {
		// Mouse button being pushed event
		case FL_PUSH:
			last_push = Fl::event_button();
			// if the left button be pushed is left mouse button
			if (last_push == FL_LEFT_MOUSE  ) {
				doPick();
				damage(1);
				return 1;
			};
			break;

	   // Mouse button release event
		case FL_RELEASE: // button release
			damage(1);
			last_push = 0;
			return 1;

		// Mouse button drag event
		case FL_DRAG:

			// Compute the new control point position
			if ((last_push == FL_LEFT_MOUSE) && (selectedCube >= 0)) {
				ControlPoint* cp = &m_pTrack->points[selectedCube];

				double r1x, r1y, r1z, r2x, r2y, r2z;
				getMouseLine(r1x, r1y, r1z, r2x, r2y, r2z);

				double rx, ry, rz;
				mousePoleGo(r1x, r1y, r1z, r2x, r2y, r2z, 
								static_cast<double>(cp->pos.x), 
								static_cast<double>(cp->pos.y),
								static_cast<double>(cp->pos.z),
								rx, ry, rz,
								(Fl::event_state() & FL_CTRL) != 0);

				cp->pos.x = (float) rx;
				cp->pos.y = (float) ry;
				cp->pos.z = (float) rz;
				damage(1);
			}
			break;

		// in order to get keyboard events, we need to accept focus
		case FL_FOCUS:
			return 1;

		// every time the mouse enters this window, aggressively take focus
		case FL_ENTER:	
			focus(this);
			break;

		case FL_KEYBOARD:
		 		int k = Fl::event_key();
				int ks = Fl::event_state();
				if (k == 'p') {
					// Print out the selected control point information
					if (selectedCube >= 0) 
						printf("Selected(%d) (%g %g %g) (%g %g %g)\n",
								 selectedCube,
								 m_pTrack->points[selectedCube].pos.x,
								 m_pTrack->points[selectedCube].pos.y,
								 m_pTrack->points[selectedCube].pos.z,
								 m_pTrack->points[selectedCube].orient.x,
								 m_pTrack->points[selectedCube].orient.y,
								 m_pTrack->points[selectedCube].orient.z);
					else
						printf("Nothing Selected\n");

					return 1;
				};
				break;
	}

	return Fl_Gl_Window::handle(event);
}

void PreView::draw()
{
	//initialized glad
	if (firstDraw) {
		if (!gladLoadGL()) {
			throw std::runtime_error("Could not initialize GLAD!");
			return;
		}
		firstDraw = false;
		myObj.InitBuffer();
		tw->textureView->openTexture("../asset/in.png");

		basicTeTextureShader.Init();

		screenToneShader.Init();
		screenToneShader.Enable();
		screenToneShader.SetTexSize(tw->textureView->loadTex.width, tw->textureView->loadTex.height);
		screenToneShader.SetWinSize(windowWidth, windowHeight);
		screenToneShader.Disable();
		// frame buffer init
		glGenFramebuffers(1, &fbo);

		openMyObjFile("D:\\workplace\\MangaVectorData\\Mesh\\starFK_1200_B.myobj");
		tw->setUIImgSize(myObj.img_width, myObj.img_height);

		setImgOutputSize(myObj.img_width, myObj.img_height);
	}

	glBindFramebuffer(GL_FRAMEBUFFER, fbo);
	render();
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	glClearColor(0, 0, .3f, 0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glEnable(GL_TEXTURE_2D);
	basicTeTextureShader.Enable();
	basicTeTextureShader.SetTexture(fbxTex.idx);
	glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
	basicTeTextureShader.Disable();

	//myObj.DrawLines();
}

void PreView::render()
{
	glClearColor(1.0, 1.0, 1.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	screenToneShader.Enable();
	screenToneShader.SetTexture(tw->textureView->loadTex.idx);
	screenToneShader.SetTriangleType(CVSystem::TriangleType::TRIANGLE_WHITE);
	myObj.DrawWhite2();
	screenToneShader.SetTriangleType(CVSystem::TriangleType::TRIANGLE_SCREENTONE);
	myObj.DrawSC2();
	screenToneShader.SetTriangleType(CVSystem::TriangleType::TRIANGLE_BLACK);
	myObj.DrawBlack();
	screenToneShader.SetTriangleType(CVSystem::TriangleType::TRIANGLE_WHITE);
	myObj.DrawWhite();
	screenToneShader.SetTriangleType(CVSystem::TriangleType::TRIANGLE_SCREENTONE);
	myObj.DrawSC();
	screenToneShader.Disable();
}

void PreView::setImgOutputSize(int width=0, int height=0)
{
	if (width > 0)
		fbxTex.width = width;
	if (height > 0)
		fbxTex.height = height;
	if (fbxTex.idx) {
		glDeleteTextures(1, &fbxTex.idx);
	}

	glGenTextures(1, &fbxTex.idx);
	glBindTexture(GL_TEXTURE_2D, fbxTex.idx);

	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, fbxTex.width, fbxTex.height, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	glBindFramebuffer(GL_FRAMEBUFFER, fbo);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, fbxTex.idx, 0);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	glBindTexture(GL_TEXTURE_2D, 0);

	screenToneShader.Enable();
	screenToneShader.SetOutImgSize(fbxTex.width, fbxTex.height);
	screenToneShader.Disable();
}

void PreView::openMyObjFile(std::string filePath)
{
	myObj.LoadFile(filePath);

	screenToneShader.Enable();
	screenToneShader.SetInImgSize(myObj.img_width, myObj.img_height);
	screenToneShader.Disable();
}

void PreView::openBitmapFile(std::string filePath)
{
	myObj.ReadFromBitmap(filePath, std::vector<std::string>{ "mask1.png" });

	screenToneShader.Enable();
	screenToneShader.SetInImgSize(myObj.img_width, myObj.img_height);
	screenToneShader.Disable();
}

void PreView::
doPick()
//========================================================================
{
	// since we'll need to do some GL stuff so we make this window as 
	// active window
	make_current();		

	// where is the mouse?
	int mx = Fl::event_x(); 
	int my = Fl::event_y();

	// get the viewport - most reliable way to turn mouse coords into GL coords
	int viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);

	// Set up the pick matrix on the stack - remember, FlTk is
	// upside down!
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity ();
	gluPickMatrix((double)mx, (double)(viewport[3]-my), 
						5, 5, viewport);

	// now draw the objects - but really only see what we hit
	GLuint buf[100];
	glSelectBuffer(100,buf);
	glRenderMode(GL_SELECT);
	glInitNames();
	glPushName(0);

	// draw the cubes, loading the names as we go
	for(size_t i=0; i<m_pTrack->points.size(); ++i) {
		glLoadName((GLuint) (i+1));
		m_pTrack->points[i].draw();
	}

	// go back to drawing mode, and see how picking did
	int hits = glRenderMode(GL_RENDER);
	if (hits) {
		// warning; this just grabs the first object hit - if there
		// are multiple objects, you really want to pick the closest
		// one - see the OpenGL manual 
		// remember: we load names that are one more than the index
		selectedCube = buf[3]-1;
	} else // nothing hit, nothing selected
		selectedCube = -1;

	printf("Selected Cube %d\n",selectedCube);
}
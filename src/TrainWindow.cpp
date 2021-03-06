#include "TrainWindow.H"
#include "PreView.H"
#include "TextureView.H"
#include "CallBacks.H"

TrainWindow::
TrainWindow(const int x, const int y) 
	: Fl_Double_Window(x,y,800,600,"Train and Roller Coaster")
{
	// make all of the widgets
	begin();	// add to this widget
	{
		int pty=5;			// where the last widgets were drawn

		preView = new PreView(5,5,590,590);
		preView->tw = this;
		preView->m_pTrack = &m_Track;
		this->resizable(preView);

		// to make resizing work better, put all the widgets in a group
		widgets = new Fl_Group(600,5,190,590);
		widgets->begin();

		//runButton = new Fl_Button(605,pty,60,20,"Run");
		//togglify(runButton);

		//Fl_Button* fb = new Fl_Button(700,pty,25,20,"@>>");
		//fb->callback((Fl_Callback*)forwCB,this);
		//Fl_Button* rb = new Fl_Button(670,pty,25,20,"@<<");
		//rb->callback((Fl_Callback*)backCB,this);
		//
		//arcLength = new Fl_Button(730,pty,65,20,"ArcLength");
		//togglify(arcLength,1);
  //
		//pty+=25;
		//speed = new Fl_Value_Slider(655,pty,140,20,"speed");
		//speed->range(0,10);
		//speed->value(2);
		//speed->align(FL_ALIGN_LEFT);
		//speed->type(FL_HORIZONTAL);

		//pty += 30;

		//// camera buttons - in a radio button group
		//Fl_Group* camGroup = new Fl_Group(600,pty,195,20);
		//camGroup->begin();
		//worldCam = new Fl_Button(605, pty, 60, 20, "World");
  //      worldCam->type(FL_RADIO_BUTTON);		// radio button
  //      worldCam->value(1);			// turned on
  //      worldCam->selection_color((Fl_Color)3); // yellow when pressed
		//worldCam->callback((Fl_Callback*)damageCB,this);
		//trainCam = new Fl_Button(670, pty, 60, 20, "Train");
  //      trainCam->type(FL_RADIO_BUTTON);
  //      trainCam->value(0);
  //      trainCam->selection_color((Fl_Color)3);
		//trainCam->callback((Fl_Callback*)damageCB,this);
		//topCam = new Fl_Button(735, pty, 60, 20, "Top");
  //      topCam->type(FL_RADIO_BUTTON);
  //      topCam->value(0);
  //      topCam->selection_color((Fl_Color)3);
		//topCam->callback((Fl_Callback*)damageCB,this);
		//camGroup->end();

		//pty += 30;

		//// browser to select spline types
		//// TODO: make sure these choices are the same as what the code supports
		//splineBrowser = new Fl_Browser(605,pty,120,75,"Spline Type");
		//splineBrowser->type(2);		// select
		//splineBrowser->callback((Fl_Callback*)damageCB,this);
		//splineBrowser->add("Linear");
		//splineBrowser->add("Cardinal Cubic");
		//splineBrowser->add("Cubic B-Spline");
		//splineBrowser->select(2);

		//pty += 110;

		//// add and delete points
		//Fl_Button* ap = new Fl_Button(605,pty,80,20,"Add Point");
		//ap->callback((Fl_Callback*)addPointCB,this);
		//Fl_Button* dp = new Fl_Button(690,pty,80,20,"Delete Point");
		//dp->callback((Fl_Callback*)deletePointCB,this);

		//pty += 25;
		//// reset the points
		//resetButton = new Fl_Button(735,pty,60,20,"Reset");
		//resetButton->callback((Fl_Callback*)resetCB,this);
		//Fl_Button* loadb = new Fl_Button(605,pty,60,20,"Load");
		//loadb->callback((Fl_Callback*) loadCB, this);
		//Fl_Button* saveb = new Fl_Button(670,pty,60,20,"Save");
		//saveb->callback((Fl_Callback*) saveCB, this);

		//pty += 25;
		//// roll the points
		//Fl_Button* rx = new Fl_Button(605,pty,30,20,"R+X");
		//rx->callback((Fl_Callback*)rpxCB,this);
		//Fl_Button* rxp = new Fl_Button(635,pty,30,20,"R-X");
		//rxp->callback((Fl_Callback*)rmxCB,this);
		//Fl_Button* rz = new Fl_Button(670,pty,30,20,"R+Z");
		//rz->callback((Fl_Callback*)rpzCB,this);
		//Fl_Button* rzp = new Fl_Button(700,pty,30,20,"R-Z");
		//rzp->callback((Fl_Callback*)rmzCB,this);
		// pty += 25;

		widthInput = new Fl_Input(650, pty, 40, 20, "witdh: ");
		widthInput->value("0");
		widthInput->callback((Fl_Callback*)setWidth, this);
		heightInput = new Fl_Input(750, pty, 40, 20, "height: ");
		heightInput->value("0");
		heightInput->callback((Fl_Callback*)setHeight, this);
		pty += 25;

		pty += 5;
		textureView = new TextureView(605, pty, 50, 50);
		textureView->tw = this;
		pty += 60;

		Fl_Button* openMyObjBtn = new Fl_Button(605, pty, 180, 20, "Open MyObj");
		openMyObjBtn->callback((Fl_Callback*)openMyObjFileDialog, this);
		pty += 30;

		Fl_Button* openBitmapjBtn = new Fl_Button(605, pty, 180, 20, "Open Bitmap");
		openBitmapjBtn->callback((Fl_Callback*)openBitmapFileDialog, this);
		pty += 30;

		Fl_Button* saveImgBtn = new Fl_Button(605, pty, 180, 20, "Save");
		saveImgBtn->callback((Fl_Callback*)saveImg, this);

		// we need to make a little phantom widget to have things resize correctly
		Fl_Box* resizebox = new Fl_Box(600,595,200,5);
		widgets->resizable(resizebox);

		widgets->end();
	}
	end();	// done adding to this widget
}

void TrainWindow::
togglify(Fl_Button* b, int val)
//========================================================================
{
	b->type(FL_TOGGLE_BUTTON);		// toggle
	b->value(val);		// turned off
	b->selection_color((Fl_Color)3); // yellow when pressed	
	b->callback((Fl_Callback*)damageCB,this);
}

void TrainWindow::setUIImgSize(int _width, int _height)
{
	widthInput->value(std::to_string(_width).c_str());
	heightInput->value(std::to_string(_height).c_str());
}

//************************************************************************
//
// *
//========================================================================
void TrainWindow::
damageMe()
//========================================================================
{
	if (preView->selectedCube >= ((int)m_Track.points.size()))
		preView->selectedCube = 0;
	preView->damage(1);
}

void TrainWindow::
advanceTrain(float dir)
{

}
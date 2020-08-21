#pragma once

#pragma warning(push)
#pragma warning(disable:4312)
#pragma warning(disable:4311)
#include "../include/Common.h"
#pragma warning(pop)

// we need to know what is in the world to show
#include "Track.H"

// other things we just deal with as pointers, to avoid circular references
class PreView;
class TextureView;

class TrainWindow : public Fl_Double_Window {
	public:
		TrainWindow(const int x=50, const int y=50);

	public:
		// call this method when things change
		void damageMe();

		// this moves the train forward on the track - its up to you to do this
		// correctly. it gets called from the idle callback loop
		// it should handle forward and backwards
		void advanceTrain(float dir = 1);

		// simple helper function to set up a button
		void togglify(Fl_Button*, int state=0);

		void setUIImgSize(int width, int height);

	public:
		// keep track of the stuff in the world
		CTrack				m_Track;

		// the widgets that make up the Window
		PreView*			preView;
		TextureView*		textureView;

		Fl_Group*			widgets;	// all widgets, grouped for resizing ease

		// utility buttons
		Fl_Button*			resetButton;

		// which viewpoint are we drawing from
		Fl_Button*			worldCam;
		Fl_Button*			trainCam;
		Fl_Button*			topCam;

		// the type of the spline (use its value to determine)
		Fl_Browser*			splineBrowser;

		// are we animating the train?
		Fl_Button*			runButton;
		// if we're animating it, how fast should it go?
		Fl_Value_Slider*	speed;
		Fl_Button*			arcLength;		// do we use arc length for speed?

		// UI
		Fl_Input*			heightInput;
		Fl_Input*			widthInput;

};

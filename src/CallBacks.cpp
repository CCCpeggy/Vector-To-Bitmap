#pragma once

#include "TrainWindow.H"
#include "PreView.H"
#include "TextureView.H"
#include "CallBacks.H"

#pragma warning(push)
#pragma warning(disable:4312)
#pragma warning(disable:4311)
#pragma warning(pop)

//***************************************************************************
//
// * Reset the control points back to their base setup
//===========================================================================
void resetCB(Fl_Widget*, TrainWindow* tw)
//===========================================================================
{
	tw->m_Track.resetPoints();
	tw->preView->selectedCube = -1;
	tw->m_Track.trainU = 0;
	tw->damageMe();
}

//***************************************************************************
//
// * any time something changes, you need to force a redraw
//===========================================================================
void damageCB(Fl_Widget*, TrainWindow* tw)
{
	tw->damageMe();
}

//***************************************************************************
//
// * Callback that adds a new point to the spline
// idea: add the point AFTER the selected point
//===========================================================================
void addPointCB(Fl_Widget*, TrainWindow* tw)
//===========================================================================
{
	// get the number of points
	size_t npts = tw->m_Track.points.size();
	// the number for the new point
	size_t newidx = (tw->preView->selectedCube>=0) ? tw->preView->selectedCube : 0;

	// pick a reasonable location
	size_t previdx = (newidx + npts -1) % npts;
	Pnt3f npos = (tw->m_Track.points[previdx].pos + tw->m_Track.points[newidx].pos) * .5f;

	tw->m_Track.points.insert(tw->m_Track.points.begin() + newidx,npos);

	// make it so that the train doesn't move - unless its affected by this control point
	// it should stay between the same points
	if (ceil(tw->m_Track.trainU) > ((float)newidx)) {
		tw->m_Track.trainU += 1;
		if (tw->m_Track.trainU >= npts) tw->m_Track.trainU -= npts;
	}

	tw->damageMe();
}

//***************************************************************************
//
// * Callback that deletes a point from the spline
//===========================================================================
void deletePointCB(Fl_Widget*, TrainWindow* tw)
//===========================================================================
{
	if (tw->m_Track.points.size() > 4) {
		if (tw->preView->selectedCube >= 0) {
			tw->m_Track.points.erase(tw->m_Track.points.begin() + tw->preView->selectedCube);
		} else
			tw->m_Track.points.pop_back();
	}
	tw->damageMe();
}
//***************************************************************************
//
// * Advancing the train
//===========================================================================
void forwCB(Fl_Widget*, TrainWindow* tw)
{
	tw->advanceTrain(2);
	tw->damageMe();
}
//***************************************************************************
//
// * Reverse the movement of the train
//===========================================================================
void backCB(Fl_Widget*, TrainWindow* tw)
//===========================================================================
{
	tw->advanceTrain(-2);
	tw->damageMe();
}



static unsigned long lastRedraw = 0;
//***************************************************************************
//
// * Callback for idling - if things are sitting, this gets called
// if the run button is pushed, then we need to make the train go.
// This is taken from the old "RunButton" demo.
// another nice problem to have - most likely, we'll be too fast
// don't draw more than 30 times per second
//===========================================================================
void runButtonCB(TrainWindow* tw)
//===========================================================================
{
	if (tw->runButton->value()) {	// only advance time if appropriate
		if (clock() - lastRedraw > CLOCKS_PER_SEC/30) {
			lastRedraw = clock();
			tw->advanceTrain();
			tw->damageMe();
		}
	}
}

//***************************************************************************
//
// * Load the control points from the files
//===========================================================================
void loadCB(Fl_Widget*, TrainWindow* tw)
//===========================================================================
{
	const char* fname = 
		fl_file_chooser("Pick a Track File","*.txt","TrackFiles/track.txt");
	if (fname) {
		tw->m_Track.readPoints(fname);
		tw->damageMe();
	}
}
//***************************************************************************
//
// * Save the control points
//===========================================================================
void saveCB(Fl_Widget*, TrainWindow* tw)
//===========================================================================
{
	const char* fname = 
		fl_input("File name for save (should be *.txt)","TrackFiles/");
	if (fname)
		tw->m_Track.writePoints(fname);
}

//***************************************************************************
//
// * Rotate the selected control point about x axis
//===========================================================================
void rollx(TrainWindow* tw, float dir)
{
	int s = tw->preView->selectedCube;
	if (s >= 0) {
		Pnt3f old = tw->m_Track.points[s].orient;
		float si = sin(((float)M_PI_4) * dir);
		float co = cos(((float)M_PI_4) * dir);
		tw->m_Track.points[s].orient.y = co * old.y - si * old.z;
		tw->m_Track.points[s].orient.z = si * old.y + co * old.z;
	}
	tw->damageMe();
} 

//***************************************************************************
//
// * Rotate the selected control point about x axis by one more degree
//===========================================================================
void rpxCB(Fl_Widget*, TrainWindow* tw)
//===========================================================================
{
	rollx(tw,1);
}
//***************************************************************************
//
// * Rotate the selected control point  about x axis by less one degree
//===========================================================================
void rmxCB(Fl_Widget*, TrainWindow* tw)
//===========================================================================
{
	rollx(tw,-1);
}

//***************************************************************************
//
// * Rotate the selected control point  about z axis
//===========================================================================
void rollz(TrainWindow* tw, float dir)
//===========================================================================
{
	int s = tw->preView->selectedCube;
	if (s >= 0) {

		Pnt3f old = tw->m_Track.points[s].orient;

		float si = sin(((float)M_PI_4) * dir);
		float co = cos(((float)M_PI_4) * dir);

		tw->m_Track.points[s].orient.y = co * old.y - si * old.x;
		tw->m_Track.points[s].orient.x = si * old.y + co * old.x;
	}

	tw->damageMe();
}

//***************************************************************************
//
// * Rotate the selected control point  about the z axis one more degree
//===========================================================================
void rpzCB(Fl_Widget*, TrainWindow* tw)
//===========================================================================
{
	rollz(tw,1);
}

//***************************************************************************
//
// *  Rotate the selected control point  about the z axis one less degree
//===========================================================================
void rmzCB(Fl_Widget*, TrainWindow* tw)
//===========================================================================
{
	rollz(tw, -1);
}

void setWidth(Fl_Input* input, TrainWindow* tw)
{
	int value = atoi(input->value());
	if (value > 0) {
		tw->preView->setImgOutputSize(value, 0);
		tw->preView->redraw();
	}
}

void setHeight(Fl_Input* input, TrainWindow* tw)
{
	int value = atoi(input->value());
	if (value > 0) {
		tw->preView->setImgOutputSize(0, value);
		tw->preView->redraw();
	}
}

void saveImg(Fl_Widget*, TrainWindow* tw) {
	TextureData* fbxTex = &tw->preView->fbxTex;
	Common::SavePng("out.png", fbxTex->idx, fbxTex->width, fbxTex->height);
}

void openMyObjFileDialog(Fl_Widget* widget, TrainWindow* tw)
{
	Fl_File_Chooser* fc = new Fl_File_Chooser(".", "Myobj Files (*.{myobj})", Fl_File_Chooser::SINGLE, "Input File");
	fc->callback(openMyObjFile, tw);
	fc->show();
}

void openBitmapFileDialog(Fl_Widget*, TrainWindow* tw)
{
	Fl_File_Chooser* fc = new Fl_File_Chooser(".", "Image Files (*.{bmp,gif,jpg,png})", Fl_File_Chooser::SINGLE, "Input File");
	fc->callback(openBitmapFile, tw);
	fc->show();
}

void openMyObjFile(Fl_File_Chooser* w, void* vtw)
{
	TrainWindow* tw = (TrainWindow*)vtw;
	tw->preView->openMyObjFile(w->value());
	tw->preView->redraw();
}

void openBitmapFile(Fl_File_Chooser* w, void* vtw)
{
	TrainWindow* tw = (TrainWindow*)vtw;
	tw->preView->openBitmapFile(w->value());
	tw->preView->redraw();
}

void openTextureFile(Fl_File_Chooser* w, void* vtw)
{
	TrainWindow* tw = (TrainWindow*)vtw;
	tw->textureView->openTexture(w->value());
	tw->textureView->redraw();
	tw->preView->redraw();
}


#include "cinder/app/AppNative.h"
#include "cinder/Camera.h"
#include "cinder/MayaCamUI.h"
#include "cinder/params/Params.h"
#include "cinder/gl/gl.h"
#include "cinder/Vector.h"
#include "StamFluidSystem.h"

using namespace ci;
using namespace ci::app;
using namespace std;

class WaterSimulationApp : public AppNative {

public:
	void prepareSettings(Settings *settings){
		//settings->enableHighDensityDisplay();
		settings->setWindowSize(800,600);
		settings->setResizable(true);
		settings->setFrameRate(60.f);
	};
	void setup();
	void resize();
	void update();
	void draw();

	void mouseDown( MouseEvent event );	
	void mouseDrag(MouseEvent event);
	void keyDown(KeyEvent event);
	

	MayaCamUI				mMayaCam;

	params::InterfaceGlRef	mParams;
	StamFluidSystem			stamfluid;
	////FlipFluidSytem		flipfluid;

	Vec3f					boxCenter;
	Vec3i					boxDimension;
	
	Vec3f					gravity;
	Vec3f					randomForce;
	Vec3f					randomForceX;
	
	bool					mAnimate;
	bool					mDrawGravity;
	bool					mAddNewForce;
	bool					mDrawBox;
	bool					mDrawGrid;
	bool					mDrawVelocity;
};

void WaterSimulationApp::setup(){

	boxDimension=Vec3f(30.0f,30.0f,30.0f);
	boxCenter=boxDimension/2;

	randomForce=Vec3f::zero();
	randomForceX=Vec3f::zero();

	//// setup our default camera, looking down the z-axis
	CameraPersp	cam;
	cam.setEyePoint(Vec3f(200, 50, 50));
	cam.setCenterOfInterestPoint(boxCenter);
	cam.setPerspective(60.0, getWindowAspectRatio(), 1.0, 1000.0);
	mMayaCam.setCurrentCam( cam);

	Matrix44f mvM=mMayaCam.getCamera().getModelViewMatrix();
	gravity = Vec3f(mvM.m10,mvM.m11,mvM.m12);//col4:camera center, row: right,up,back
	gravity*= 9.8f;
	// Setup the parameters
	mDrawGravity=false;
	mAddNewForce=false;
	mAnimate	=false;
	mDrawBox	=true;

	mParams = params::InterfaceGl::create( getWindow(), "Water Simulation Parameters", toPixels( Vec2i( 300,100  ) ) );
	mParams->addParam("Draw Box", &mDrawBox);
	mParams->addParam("Draw Grid",&mDrawGrid);
	mParams->addSeparator();
	mParams->addParam("External Force Position",&randomForceX);
	mParams->addParam("External Force Dir",&randomForce);
	mParams->addParam("Draw Gravity",&mDrawGravity);
	mParams->addParam("Draw Velocity", &mDrawVelocity);
	mParams->addSeparator();
	mParams->addParam("Animate",&mAnimate);
	mParams->addText("status","label=` `");
	mParams->addParam("Time Elapsed",&stamfluid.elapsed,"",true);
//	Fluid System setup
	stamfluid.reset(boxDimension);
}


void WaterSimulationApp::mouseDown( MouseEvent event ){
	mMayaCam.mouseDown(event.getPos());
}
void WaterSimulationApp::mouseDrag(MouseEvent event){
	mMayaCam.mouseDrag(event.getPos(),event.isLeftDown(),event.isMiddleDown(),event.isRightDown());
	Matrix44f mvM=mMayaCam.getCamera().getModelViewMatrix();
	gravity = Vec3f(mvM.m10,mvM.m11,mvM.m12);//col4:camera center, row: right,up,back
	gravity*= 9.8f;
}
void WaterSimulationApp::keyDown( KeyEvent event ){
	switch(event.getChar()){
	case 'r':stamfluid.reset(boxDimension);break;
	case 'f':setFullScreen( ! isFullScreen() ); break;
	}
	switch(event.getCode()){
	case KeyEvent::KEY_SPACE:mAnimate=!mAnimate;break;
	case KeyEvent::KEY_ESCAPE:quit();break;
	}
	
}


void WaterSimulationApp::update(){
	////change status
	if(mAnimate){
		stamfluid.step(1.f/60.0f);
		mParams->setOptions("status","label=`In animation:`");
	}else{
		mParams->setOptions("status","label=`Idle:`");
	}
}

void WaterSimulationApp::resize(){
	CameraPersp cam = mMayaCam.getCamera();
	cam.setAspectRatio( getWindowAspectRatio() );
	mMayaCam.setCurrentCam( cam );

}

void WaterSimulationApp::draw(){

	gl::clear( Color( 0.11f, 0.13f, 0.1f ) );
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	mParams->draw();

	gl::setMatrices(mMayaCam.getCamera());
	if(mDrawBox){
		glColor4f(1, 0.75, 0.5, 1);//orange
		gl::drawStrokedCube(boxCenter,boxDimension);
	}
	if(mDrawGravity){
		glColor4f(1.0f,0.3f,0.0f,0.5f);//red
		if(mDrawGravity)glColor4f(0.0f,1.0f,0.5f,1);//green
		gl::drawVector(boxCenter,(boxCenter-gravity),3.0f,1.0f);
	}
	
	//draw random force 
	if(randomForce.lengthSquared()!=0){
		glColor4f(1.0f,0.0f,0.0f,1);
		gl::drawVector(randomForceX,randomForce,1.0f,0.5f);
	}
}



CINDER_APP_NATIVE( WaterSimulationApp, RendererGl )
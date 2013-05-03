/*
Stam Semi Lagrangian Fluid Grid System
[].Stam
[].Stagger Grid
[].Shader
*/
#pragma once 
#include "Array3.h"
#include "Source.h"
#include "cinder/app/AppBasic.h"
#include <vector>
using namespace ci;
using namespace ci::app;
using namespace std;

class StamFluidSystem{

public :
	StamFluidSystem(){}
	~StamFluidSystem(){
		delete u0;delete v0; delete density0;delete pressure;
		delete w0;delete u1; delete density1;
		delete v1;delete w1; delete divergence;
	}
	ci::Vec3i	gridDim;
	ci::Vec3i	cellDim;
	double		elapsed;

	Array3f *u0,*v0,*w0,
		*u1,*v1,*w1;
	Array3f *density0,*density1;
	Array3f *divergence, *pressure;
	vector<Source> sources;

	Vec3f getVelocity(Vec3f pos);
	float interpolate(Vec3f pos,  const Array3f *s);

	void step(float dt);
	void reset(ci::Vec3i dim){
		console()<<"RESET: stam fluid sytem.-----"<<endl;
		gridDim=dim;
		cellDim=Vec3i(1,1,1);
		gridDim+=2;//add boundary cell
		u0=new Array3f(gridDim.x,gridDim.y,gridDim.z);
		u1=new Array3f(gridDim.x,gridDim.y,gridDim.z);
		v0=new Array3f(gridDim.x,gridDim.y,gridDim.z);
		v1=new Array3f(gridDim.x,gridDim.y,gridDim.z);
		w0=new Array3f(gridDim.x,gridDim.y,gridDim.z);
		w1=new Array3f(gridDim.x,gridDim.y,gridDim.z);

		density0=new Array3f(gridDim.x,gridDim.y,gridDim.z);
		density1=new Array3f(gridDim.x,gridDim.y,gridDim.z);

		divergence=new Array3f(gridDim.x,gridDim.y,gridDim.z);
		pressure=new Array3f(gridDim.x,gridDim.y,gridDim.z);
		gridDim-=2;//remove boundary cell.
		if(!sources.empty())sources.clear();
		elapsed=0;
	}

};



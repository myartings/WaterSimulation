/*
Stam Semi Lagrangian Fluid Grid System
[].Stam
[].Stagger Grid
[].Shader
*/
#pragma once 
#include "Array3.h"
#include "cinder/app/AppBasic.h"
#include "FluidSystem.h"
using namespace ci;
using namespace ci::app;
using namespace std;

class StamFluidSystem:public FluidSystem{

public :
	StamFluidSystem(){}
	~StamFluidSystem(){}

	Array3fRef u1,v1,w1;

	

private:
	
	void diffuse(int bnd, Array3fRef &x, const Array3fRef&x0,float diff,float dt);
	void transport(int bnd, Array3fRef &x,Array3fRef &x0, Array3fRef &u,Array3fRef &v,Array3fRef &w, float dt);
	void project(Array3fRef &u,Array3fRef &v,Array3fRef &w,Array3fRef &div,Array3fRef &p);
	void linear_solve(int bnd,Array3fRef &x, const Array3fRef &x0, float a, float div);
	void set_boundary(int bnd,Array3fRef &x);
	void swap_velocity();

	void  getVelocity(Vec3f pos,Vec3f &vel);
	void  traceParticle(Vec3f x0,float h, Vec3f &x1);
	
	//virtual 
	void reset_derived(Vec3i dim){
		console()<<"stam fluid system reset"<<endl;
		u1.reset(new Array3f(gridDim.x,gridDim.y,gridDim.z));
		v1.reset(new Array3f(gridDim.x,gridDim.y,gridDim.z));
		w1.reset(new Array3f(gridDim.x,gridDim.y,gridDim.z));
	}
	void interpolate_index(Vec3f &pos,Vec3i &index, Vec3f &ifloat);
	void source_step(float dt);
	void velocity_step(float dt);
	void scalar_step(float dt);
	void drawVelocity();
};



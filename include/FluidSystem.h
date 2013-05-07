#pragma once
#include "cinder/Vector.h"
#include "Array3.h"
#include "Source.h"
#include <vector>
using namespace ci;
using namespace ci::app;
using namespace std;

class FluidSystem{

public:
	FluidSystem(){}
	
	Vec3i	gridDim;
	Vec3i	cellDim;
	virtual ~FluidSystem(){deleteMemory();}
	double	elapsed;

	bool					mDrawVelocity;
	float					viscosity;
	float					diffusion;
	float					buoyancy;
	int						iterations;

	unique_ptr<Array3f> u0,v0,w0;
	unique_ptr<Array3f> temperature0,temperature1;
	unique_ptr<Array3f> divergence, pressure;
	vector<Source> sources;

	void step(float dt){
		source_step(dt);
		velocity_step(dt);
		scalar_step(dt);
		elapsed+=dt;
	}
	void reset(Vec3i dim){

		viscosity=1e-6f;
		buoyancy=0.1f;
		diffusion=1e-6f;
		iterations=4;
		elapsed=0;

		mDrawVelocity =true;

		gridDim=dim;
		cellDim=Vec3i(1,1,1);
		gridDim+=2;//add boundary cell		
		
		reset_derived(dim);
		console()<<"fluid system reset"<<endl;
		
		u0.reset(new Array3f(gridDim.x,gridDim.y,gridDim.z));
		v0.reset(new Array3f(gridDim.x,gridDim.y,gridDim.z));
		w0.reset(new Array3f(gridDim.x,gridDim.y,gridDim.z));
		
		temperature0.reset(new Array3f(gridDim.x,gridDim.y,gridDim.z));
		temperature1.reset(new Array3f(gridDim.x,gridDim.y,gridDim.z));
		
		divergence.reset(new Array3f(gridDim.x,gridDim.y,gridDim.z));
		pressure.reset(new Array3f(gridDim.x,gridDim.y,gridDim.z));
	
		gridDim-=2;//remove boundary cell.
		if(!sources.empty())sources.clear();
		
	}
	void deleteMemory(){
		sources.clear();
	}
	void draw(){
		if(mDrawVelocity)drawVelocity();	
	}
	
	void addSource(Vec3f x,float dt,unique_ptr<Array3f> &s,float force){
	
		Vec3i index;
		Vec3f ifloat;
		//TODO merge the calculation for u0,v0,w0
	
		interpolate_index(x,index,ifloat);
		if (index.x>=0 && index.x<s->nx){ 
			if(index.y>=0 && index.y<s->ny){
				if(index.z>=0 && index.z<s->nz)//i,j,k
					(*s)(index)+=(1-ifloat.x)*(1-ifloat.y)*(1-ifloat.z)*dt*force;
				if(index.z+1>=0 && index.z+1<s->nz)//i,j,k+1
					(*s)(index.x,index.y,index.z+1)+= (1-ifloat.x)*(1-ifloat.y)*ifloat.z* dt *force;

			}else{if(index.y+1>=0 && index.y+1<s->ny){
				if(index.z>=0 && index.z<s->nz)//i,j+1,k
					(*s)(index.x,index.y+1,index.z)+= (1-ifloat.x)*ifloat.y*(1-ifloat.z)* dt *force;
				if(index.z+1>=0 && index.z+1<s->nz)//i,j+1,k+1
					(*s)(index.x,index.y+1,index.z+1)+= (1-ifloat.x)*ifloat.y*ifloat.z* dt *force;	
				}
			}

		}else{if (index.x+1>=0 && index.x+1<s->nx){ 
				if(index.y>=0 && index.y<s->ny){
					if(index.z>=0 && index.z<s->nz)//i+1,j,k
						(*s)(index.x+1,index.y,index.z) += ifloat.x*(1-ifloat.y)*(1-ifloat.z)*dt*force;	
					if(index.z+1>=0 && index.z+1<s->nz)//i+1,j,k+1
						(*s)(index.x+1,index.y,index.z+1)+= ifloat.x*(1-ifloat.y)*ifloat.z* dt *force;
				}else{if(index.y+1>=0 && index.y+1<s->ny){
					if(index.z>=0 && index.z<s->nz)//i+1,j+1,k
						(*s)(index.x+1,index.y+1,index.z) += ifloat.x*ifloat.y*(1-ifloat.z)*dt*force;
					if(index.z+1>=0 && index.z+1<s->nz)//i+1,j+1,k+1
						(*s)(index.x+1,index.y+1,index.z+1) += ifloat.x*ifloat.y*ifloat.z*dt *force;
					}
				}
			}
		}
    		
	}

private:
	virtual void reset_derived(Vec3i dim)=0;
	virtual void interpolate_index(Vec3f &pos,Vec3i &index, Vec3f &ifloat)=0;
	
	virtual void source_step(float dt)=0;
	virtual void velocity_step(float dt)=0;
	virtual void scalar_step(float dt)=0;

	virtual void drawVelocity()=0;
};
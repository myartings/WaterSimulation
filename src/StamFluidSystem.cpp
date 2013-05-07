#include "StamFluidSystem.h"


void StamFluidSystem::interpolate_index(Vec3f &pos,Vec3i &index, Vec3f &ifloat)
{
	index.x =(int)floor(pos.x/cellDim.x - 0.5 );
	index.y =(int)floor(pos.y/cellDim.y - 0.5 );
	index.z =(int)floor(pos.z/cellDim.z - 0.5 );
	//it minus 0.5 because the i,j was position in the center of the cell.
	//if directly - i,j,k, produce value > =1 
	ifloat.x = pos.x/cellDim.x - 0.5f - index.x;
	ifloat.y = pos.y/cellDim.y - 0.5f - index.y;
	ifloat.z = pos.z/cellDim.z - 0.5f - index.z;
}
/*return velocity at given point using interpolation on u0,v0,w0. 
@param pos.
@return: interpolate velocity
*/
void StamFluidSystem::getVelocity( Vec3f pos,Vec3f &vel )
{
	Vec3i index;
	Vec3f ifloat;
	interpolate_index(pos,index,ifloat);
	/*vel.x=(*u0)(Vec3f(pos));
	vel.y=(*v0)(Vec3f(pos));
	vel.z=(*w0)(Vec3f(pos));
	*/
	vel.x=u0->trilerp(index.x,index.y,index.z,ifloat.x,ifloat.y, ifloat.z);
	vel.y=v0->trilerp(index.x,index.y,index.z,ifloat.x,ifloat.y, ifloat.z);
	vel.z=w0->trilerp(index.x,index.y,index.z,ifloat.x,ifloat.y, ifloat.z);
	
}
void StamFluidSystem::traceParticle( Vec3f x0,float h, Vec3f &x1 )
{
	Vec3f vel = Vec3f::zero();
	x1.set( x0 );
	getVelocity(x1,vel);
	vel*=h;      
	x1+=vel;
}



void StamFluidSystem::source_step(float dt){
	for(Source s : sources){
		//TODO:add contribution
	}
}
void StamFluidSystem::velocity_step( float dt )
{
	diffuse(1, u1, u0, viscosity, dt);
	diffuse(2, v1, v0, viscosity, dt);
	diffuse(3, w1, w0, viscosity, dt);
	swap_velocity();
	/*		//project(u,v,u0,v0);
	project(u, v, divergence, pressure);
	float [][]temp=U0;
	U0=U1;
	U1=temp;
	transport(1, u, u0, u0, v0, dt);
	transport(2, v, v0, u0, v0, dt);
	//        project(u, v, divergence, pressure);

	// project(u, v, u0,v0);
	*/
}
void StamFluidSystem::scalar_step( float dt )
{

}

void StamFluidSystem::drawVelocity()
{
	Vec3f velo=Vec3f::zero(),pos;
	//cell center
	glColor3f(0,1,0);
	for(int i=0;i<=gridDim.x+1;i++){
		for(int j=0;j<=gridDim.y+1;j++){
			for(int k=0;k<=gridDim.z+1;k++){
				pos=Vec3f(i,j,k);
				pos*=cellDim;
				getVelocity(pos,velo);
				if(velo.lengthSquared()!=0){
					pos-=0.5f;
					velo+=pos;
					gl::drawVector(pos,velo);
				}
			}
		}
	}

	/*//boundary limit
	glColor3f(1,1,1);
	gl::drawLine(Vec3f(-0.5f,-0.5f,-0.5f),Vec3f(0.5f+gridDim.x,-0.5f,-0.5f));
	gl::drawLine(Vec3f(-0.5f,-0.5f,-0.5f),Vec3f(-0.5f,0.5f+gridDim.y,-0.5f));
	gl::drawLine(Vec3f(-0.5f,-0.5f,-0.5f),Vec3f(-0.5f,-0.5f,0.5f+gridDim.z));
	*/
	//stagger cell border
}

void StamFluidSystem::diffuse( int bnd, unique_ptr<Array3f> &x, const unique_ptr<Array3f>&x0,float diff,float dt )
{
	float a1 = dt * diff *gridDim.x*gridDim.y,
		a2 = dt * diff *gridDim.x*gridDim.z,
		a3 = dt * diff *gridDim.y*gridDim.z,
		av = (a1+a2+a3)/3.0f,
		div=(1.0f+6.0f*av);
	linear_solve(bnd, x,x0,av,div,dt);
}

void StamFluidSystem::linear_solve(int bnd, unique_ptr<Array3f> &x, const unique_ptr<Array3f> &x0, float a, float div,float dt )
{
	int i,j,k,m;		
	for (m = 0; m < iterations; m++) {
		for (i = 1; i <= gridDim.x; i++) {
			for (j = 1; j <= gridDim.y; j++) {
				for(k=1; k<=gridDim.z; k++){

					(*x)(i,j,k) = ((*x0)(i,j,k) + a
						* (	 (*x)(i-1, j,   k)
						+ (*x)(i+1, j,   k)
						+ (*x)(i,   j-1, k)
						+ (*x)(i,   j+1, k)
						+ (*x)(i,   j,   k-1)
						+ (*x)(i,   j,   k+1)
						))
						/div;
				}
			}
		}
		set_boundary(bnd, x);
	}
}

void StamFluidSystem::set_boundary( int bnd,unique_ptr<Array3f> &x )
{

	//x->n* = N*+2, 
	int i,j,k;
	for (j = 1; j <x->ny-1; j++) {
		for(k=1;k <x->nz-1; k++){
			(*x)(0, j,k) = (bnd == 1) ? -(*x)(1, j,k) : (*x)(1, j,k);
			(*x)(x->nx-1, j,k) = (bnd == 1) ? -(*x)(x->nx-2,j,k): (*x)(x->nx-2,j,k);
		}
	}
	for (i = 1; i <x->nx-1; i++) {
		for(k=1;k <x->nz-1; k++){
			(*x)(i,0,k) = (bnd == 2) ? -(*x)(i,1,k) : (*x)(i,1,k);
			(*x)(i,x->ny-1,k) = (bnd == 2) ? -(*x)(i,x->ny-2,k): (*x)(i,x->ny-2,k);
		}
	}	
	for (j = 1; j <x->ny-1; j++) {
		for(i=1;i <x->nx-1; i++){
			(*x)(i,j,0) = (bnd == 3) ? -(*x)(i,j,1) : (*x)(i,j,1);
			(*x)(i,j,x->nz-1) = (bnd == 3) ? -(*x)(i,j,x->nz-2): (*x)(i,j,x->nz-2);
		}
	}	


	(*x)(0,      0,      0      ) = 0.33f * ((*x)(1,0,0) + (*x)(0,1,0)+ (*x)(0,0,1));
	(*x)(x->nx-1,0,      0      ) = 0.33f * ((*x)(x->nx-1,1,0) + (*x)(x->nx-2,0,0)+(*x)(x->nx-1,0,1));
	(*x)(0,      x->ny-1,0      ) = 0.33f * ((*x)(1,x->ny-1,0) + (*x)(0,x->ny-2,0)+(*x)(0,x->ny-1,1));
	(*x)(0,      0,      x->nz-1) = 0.33f * ((*x)(0,1,x->nz-1) + (*x)(0,0,x->nz-2)+(*x)(1,0,x->nz-1));


	(*x)(x->nx-1,x->ny-1, 0      ) = 0.33f * ((*x)(x->nx-2,x->ny-1,0) + (*x)(x->nx-1,x->ny-2,0)+(*x)(x->nx-1,x->ny-1,1));
	(*x)(x->nx-1,0,       x->nz-1) = 0.33f * ((*x)(x->nx-2,0,x->nz-1) + (*x)(x->nx-1,0,x->nz-2)+(*x)(x->nx-1,1,x->nz-1));
	(*x)(0,      x->ny-1, x->nz-1) = 0.33f * ((*x)(0,x->ny-1,x->nz-2) + (*x)(0,x->ny-2,x->nz-1)+(*x)(1,x->ny-1,x->nz-1));
	(*x)(x->nx-1,x->ny-1, x->nz-1) = 0.33f * ((*x)(x->nx-2,x->ny-1,x->nz-1) + (*x)(x->nx-1,x->ny-2,x->nz-1)+(*x)(x->nx-1,x->ny-1,x->nz-2));

}

void StamFluidSystem::swap_velocity()
{	u0.swap(u1);
v0.swap(v1);
w0.swap(w1);
}






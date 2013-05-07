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
	project(u1, v1, w1, divergence, pressure);
	swap_velocity();
	transport(1, u1, u0, u0, v0, w0, dt);
	transport(2, v1, v0, u0, v0, w0, dt);
	transport(3, w1, w0, u0, v0, w0, dt);
	//project(u1, v1, w1, divergence, pressure);
	
}
void StamFluidSystem::scalar_step( float dt )
{

}

void StamFluidSystem::drawVelocity()
{
	Vec3f velo=Vec3f::zero(),pos;
	//cell center
	glColor3f(0,1,0);
	for(int i=11;i<=20;i++){
		for(int j=11;j<=20;j++){
			for(int k=11;k<=20;k++){
				pos=Vec3f(i,j,k);
				pos*=cellDim;
				getVelocity(pos,velo);
				if(velo.lengthSquared()!=0){
					pos-=0.5f;
					velo+=pos;
					gl::drawVector(pos,velo,0.1f,0.02f);
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

void StamFluidSystem::diffuse( int bnd, Array3fRef &x, const Array3fRef&x0,float diff,float dt )
{
	float a1 = dt * diff *gridDim.x*gridDim.y,
		a2 = dt * diff *gridDim.x*gridDim.z,
		a3 = dt * diff *gridDim.y*gridDim.z,
		av = (a1+a2+a3)/3.0f,
		div=(1.0f+6.0f*av);
	linear_solve(bnd, x,x0,av,div);
}
void StamFluidSystem::linear_solve(int bnd, Array3fRef &x, const Array3fRef &x0, float a, float div)
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
void StamFluidSystem::set_boundary( int bnd,Array3fRef &x )
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

void StamFluidSystem::project( Array3fRef &u,Array3fRef &v,Array3fRef &w,Array3fRef &div,Array3fRef &p )
{

	int i, j, k, m;
	for (i = 1; i <=gridDim.x; i++) {
		for (j = 1; j <=gridDim.y; j++){
			for(k=1;k<=gridDim.z;k++){
				(*div)(i,j,k) = -0.5f*
					(((*u)(i + 1, j,k)
					- (*u)(i - 1, j,k))/gridDim.x 
					+((*v)(i, j + 1,k)
					- (*v)(i, j - 1,k))/gridDim.y
					+((*w)(i, j,  k+1)
					- (*w)(i, j,  k-1))/gridDim.z);
				(*p)(i,j,k) = 0;
			}
		}
	}

	set_boundary(0, div);
	set_boundary(0, p);
	linear_solve(0,p,div,1.0f,6.0f);

	for (i = 1; i <= gridDim.x; i++) {
		for (j = 1; j <= gridDim.y; j++) {
			for(k=1;k<=gridDim.z;k++){
				(*u)(i, j,k) -= 0.5 * ((*p)(i + 1,j,k )-(*p)(i - 1,j,k))*gridDim.x;
				(*v)(i, j,k) -= 0.5 * ((*p)(i, j + 1,k)-(*p)(i, j - 1,k))*gridDim.y;
				(*w)(i, j,k) -= 0.5 * ((*p)(i, j,  k+1)-(*p)(i, j, k-1))*gridDim.z;
			}
		}
	}
	set_boundary(1, u);
	set_boundary(2, v);
	set_boundary(3, w);

}
void StamFluidSystem::transport( int bnd, Array3fRef &x,Array3fRef &x0, Array3fRef &u,Array3fRef &v,Array3fRef &w, float dt )
{ 
	int i,j,k;
	Vec3i index;
	Vec3f p0,p1=Vec3f::zero(),ifloat;
	for(i=1;i<=gridDim.x;i++){
		for(j=1;j<=gridDim.y;j++){
			for(k=1;k<=gridDim.z;k++){
				p0=Vec3f(i,j,k);
				p0+=0.5f;
				p0*=cellDim;
				traceParticle(p0,-dt,p1);
				this->interpolate_index(p1,index,ifloat);
				(*x)(i,j,k)=x0->trilerp(index.x,index.y,index.z,ifloat.x,ifloat.y,ifloat.z);
			}
		}
	}
	set_boundary(bnd,x);
}






#include "StamFluidSystem.h"


float StamFluidSystem::interpolate( Vec3f pos, const Array3f *s )
{
	float result,fx,fy,fz;
	int i,j,k;

	i =(int)floor(pos.x/cellDim.x - 0.5 );
	j =(int)floor(pos.y/cellDim.y - 0.5 );
	k =(int)floor(pos.z/cellDim.z - 0.5 );
	//it minus 0.5 because the i,j was position in the center of the cell.
	//if directly - i,j,k, produce value > = 1
	fx = pos.x/cellDim.x - 0.5f - i;
	fy = pos.y/cellDim.y - 0.5f - j;
	fz = pos.z/cellDim.z - 0.5f - k;
	result = s->trilerp(i,j,k,fx,fy,fz);
	return result;
}

/*return velocity at given point using interpolation on u0,v0,w0. 
@param pos.
@return: interpolate velocity
*/
cinder::Vec3f StamFluidSystem::getVelocity( Vec3f pos )
{
	Vec3f vel=Vec3f::zero();
	vel.x=interpolate(pos,u0);
	vel.y=interpolate(pos,v0);
	vel.z=interpolate(pos,w0);
	return vel;
}

void StamFluidSystem::step( float dt )
{
	elapsed+=dt;
}


#ifndef SOURCE_H
#define SOURCE_H
#include "cinder/Vector.h"
using namespace ci;

struct Source{
	/*position of heat source*/
	Vec3f location;
	/*heating and cooling rate*/
	float amount;
	Source(void){}
	Source(Vec3f location_, float amount_){
		location=location_;
		amount=amount_;
	}
};

#endif
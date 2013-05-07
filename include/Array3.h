/* Reference: Cinder Wave Simulation Array3 Class
Modification: For 3D and Staggered Grid 
*/
#ifndef ARRAY3
#define ARRAY3

#include <cstdio>
#include <cmath>
#include <cstring>
#include "cinder/Vector.h"

template<class T>
struct Array3{
	int nx,ny,nz;
	int size;
	T *data;

	Array3()
		:nx(0), ny(0),nz(0),size(0), data(0)
	{}

	Array3(int nx_, int ny_, int nz_)
		:nx(0), ny(0),nz(0),size(0), data(0)
	{ init(nx_, ny_, nz_ ); }

	void init(int nx_, int ny_, int nz_)
	{
		delete_memory();
		nx=nx_;
		ny=ny_;
		nz=nz_;
		size=nx*ny*nz;
		data=new T[size];
		zero();
	}

	~Array3()
	{delete_memory(); }

	void delete_memory()
	{
		delete[] data;
		data=0;
		nx=ny=nz=size=0;
	}

	const T &operator()(int i, int j, int k) const
	{ return data[i+j*nx*nz+k*nx]; }
	//x:row, y:column, z extra
	T &operator()(int i, int j,int k)
	{ return data[i+j*nx*nz+k*nx]; }

	const T &operator()(ci::Vec3i index) const
	{ return data[index.x+index.y*nx*nz+index.z*nx]; }
	//x:row, y:column, z extra
	T &operator()(ci::Vec3i index)
	{ return data[index.x+index.y*nx*nz+index.z*nx]; }


	//trilinear interpolation
	T trilerp(int i, int j, int k, T fx, T fy, T fz) const
	{ 
		float result =0;
		//TODO separate this into if statement
		/*if (i>=0 && i<this->nx){ 
			if(j>=0 && j<this->ny){
				if(k>=0 && k<this->nz)//i,j,k
					result +=(1-fx)*(1-fy)*(1-fz)*(*this)(i,j,k);
				if(k+1>=0 && k+1<this->nz)//i,j,k+1
					result+= (1-fx)*(1-fy)*fz*(*this)(i,j,k+1);

			}else{if(j+1>=0 && j+1<this->ny){
				if(k>=0 && k<this->nz)//i,j+1,k
					result+= (1-fx)*fy*(1-fz)*(*this)(i,j+1,k);
				if(k+1>=0 && k+1<this->nz)//i,j+1,k+1
					result+= (1-fx)*fy*fz* (*this)(i,j+1,k+1);	
			}
			}

		}else{if (i+1>=0 && i+1<this->nx){ 
			if(j>=0 && j<this->ny){
				if(k>=0 && k<this->nz)//i+1,j,k
					result += fx*(1-fy)*(1-fz)*(*this)(i+1,j,k);	
				if(k+1>=0 && k+1<this->nz)//i+1,j,k+1
					result+= fx*(1-fy)*fz*(*this)(i+1,j,k+1);
			}else{if(j+1>=0 && j+1<this->ny){
				if(k>=0 && k<this->nz)//i+1,j+1,k
					result+= fx*fy*(1-fz)*(*this)(i+1,j+1,k);
				if(k+1>=0 && k+1<this->nz)//i+1,j+1,k+1
					result+= fx*fy*fz*(*this)(i+1,j+1,k+1);
			}
			}
			}
		}*/


		/* if(i>=0&&j>=0&&k>=0&&i<this->nx&&j<this->ny&&k<this->nz){
		result=
		(1-fx)*((1-fy)*((1-fz)*(*this)(i,j,k)
		+fz*((k<this->nz-1)?((*this)(i,j,k+1)):0))
		+fy*((1-fz)*((j<this->ny-1)?((*this)(i,j+1,k)):0)
		+fz*((j<this->ny-1&&k<this->nz-1)?((*this)(i,j+1,k+1)):0)))+
		fx*((1-fy)*((1-fz)*((i<this->nx-1)?((*this)(i+1,j,k)):0)
		+fz*((i<this->nx-1&&k<this->nz-1)?((*this)(i+1,j,k+1)):0))
		+fy*((1-fz)*((i<this->nx-1&&j<this->ny-1)?((*this)(i+1,j+1,k)):0)
		+fz*((i<this->nx-1&&j<this->ny-1&&k<this->nz-1)?((*this)(i+1,j+1,k+1)):0)));
		*/			 
		result = (1-fx)*((1-fy)*((1-fz)*(*this)(i,j,k)+fz*(*this)(i,j,k+1))
		+fy*((1-fz)*(*this)(i,j+1,k)+fz*(*this)(i,j+1,k+1)))+
		fx*((1-fy)*((1-fz)*(*this)(i+1,j,k)+fz*(*this)(i+1,j,k+1))
		+fy*((1-fz)*(*this)(i+1,j+1,k)+fz*(*this)(i+1,j+1,k+1)));
		
		return  result;

	}

	void copy_to(Array3 &a) const
	{ std::memcpy(a.data, data, size*sizeof(T)); }

	//TODO: what does infnorm do?
	T infnorm() const
	{ 
		T r=0;
		for(int i=0; i<size; ++i)
			if(!(std::fabs(data[i])<=r))
				r=std::fabs(data[i]);
		return r;
	}

	//set default data to zero
	void zero()
	{ std::memset(data, 0, size*sizeof(T)); }

	//where could you dot product a whole 3D matrix?
	double dot(const Array3 &a) const
	{
		double r=0;
		for(int i=0; i<size; ++i)
			r+=data[i]*a.data[i];
		return r;
	}

	void increment(double scale, const Array3 &a)
	{ for(int i=0; i<size; ++i) data[i]+=scale*a.data[i]; }

	void scale_and_increment(double scale, const Array3 &a)
	{ for(int i=0; i<size; ++i) data[i]=scale*data[i]+a.data[i]; }

	T average(){
		T sum=0;
		for(int i=0;i<size;++i) sum+=data[i];
		return sum/size;
	}
	/*void write_matlab(FILE *fp, const char *variable_name)
	{
	fprintf(fp, "%s=[", variable_name);
	for(int i=0; i<nx; ++i){
	for(int j=0; j<ny; ++j)
	fprintf(fp, "%lg ", (double)(*this)(i,j));
	fprintf(fp, ";\n");
	}
	fprintf(fp,"];\n");
	}*/
};

typedef Array3<float> Array3f;
typedef Array3<double> Array3d;
typedef Array3<char> Array3c;

typedef std::unique_ptr<struct Array3<float>> Array3fRef;
#endif

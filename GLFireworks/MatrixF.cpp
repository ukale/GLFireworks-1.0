/*
 * Copyright 2010 Kinetic Bytes
 */

#include <OpenGLES/ES1/gl.h>
#include "MatrixF.h"

////////////////////////////////////////////////////
#if defined(WIN32) && defined(_DEBUG)
	#define _CRTDBG_MAP_ALLOC
	#include <stdlib.h>
	#include <crtdbg.h>
	#define DEBUG_NEW new( _NORMAL_BLOCK, __FILE__, __LINE__ )
	#define new DEBUG_NEW
#endif
////////////////////////////////////////////////////

MatrixF MatrixF::orthoMatrix;
MatrixF MatrixF::perspectiveMatrix;
MatrixF MatrixF::mvMatrix;       //Model View
MatrixF MatrixF::pMatrix;        //Projection
MatrixF MatrixF::mvpMatrix;      //Model View Projection
std::stack<MatrixF> MatrixF::mvMatrixStack;
std::stack<MatrixF> MatrixF::pMatrixStack;

////////////////////////////////////////////////////
//MatrixF implementation

static float identityM[16] = {
    1,0,0,0,
    0,1,0,0,
    0,0,1,0,
    0,0,0,1
};

//set a 4x4 Identity Matrix into the passed in float array
void MatrixF::setIdentityM(float* f) 
{
    memcpy(f, identityM, 16*sizeof(float));
    //memset(f, 0, 16*sizeof(float));
    //f[0] = f[5] = f[10] = f[15] = 1.0f;
}

//Multiply two 4x4 matrices in column major order
void MatrixF::multiplyMM(float* rm, float* m1, float* m2)
{
    float d=m1[0],e=m1[1],g=m1[2],f=m1[3],h=m1[4],i=m1[5],j=m1[6],k=m1[7],l=m1[8],o=m1[9],m=m1[10],n=m1[11],p=m1[12],r=m1[13],s=m1[14],a=m1[15];
    float A=m2[0],B=m2[1],t=m2[2],u=m2[3],v=m2[4],w=m2[5],x=m2[6],y=m2[7],z=m2[8],C=m2[9],D=m2[10],E=m2[11],q=m2[12],F=m2[13],G=m2[14],b=m2[15];
    
    rm[0]=A*d+B*h+t*l+u*p;  rm[1]=A*e+B*i+t*o+u*r;  rm[2]=A*g+B*j+t*m+u*s;  rm[3]=A*f+B*k+t*n+u*a;
    rm[4]=v*d+w*h+x*l+y*p;  rm[5]=v*e+w*i+x*o+y*r;  rm[6]=v*g+w*j+x*m+y*s;  rm[7]=v*f+w*k+x*n+y*a;
    rm[8]=z*d+C*h+D*l+E*p;  rm[9]=z*e+C*i+D*o+E*r;  rm[10]=z*g+C*j+D*m+E*s; rm[11]=z*f+C*k+D*n+E*a;
    rm[12]=q*d+F*h+G*l+b*p; rm[13]=q*e+F*i+G*o+b*r; rm[14]=q*g+F*j+G*m+b*s; rm[15]=q*f+F*k+G*n+b*a;

/*    
	//returns result in rm
	for (int i=0; i<4; i++)
	{
		float f0 = m1[i], f1 = m1[4+i], f2 = m1[8+i], f3 = m1[12+i];
		rm[i]    = f0*m2[0]  + f1*m2[1]  + f2*m2[2]  + f3*m2[3];
		rm[4+i]  = f0*m2[4]  + f1*m2[5]  + f2*m2[6]  + f3*m2[7];
		rm[8+i]  = f0*m2[8]  + f1*m2[9]  + f2*m2[10] + f3*m2[11];
		rm[12+i] = f0*m2[12] + f1*m2[13] + f2*m2[14] + f3*m2[15];
	}
*/
}

//Multiply a 4x4 Matrix in column major order with a 4x1 column vector
void MatrixF::multiplyMV(float* rv, float* m, float* v)
{
	//returns result in rv
	float f0 = v[0], f1 = v[1], f2 = v[2], f3 = v[3];
	rv[0]  = m[0]*f0 + m[4]*f1 + m[8]*f2  + m[12]*f3;
	rv[1]  = m[1]*f0 + m[5]*f1 + m[9]*f2  + m[13]*f3;
	rv[2]  = m[2]*f0 + m[6]*f1 + m[10]*f2 + m[14]*f3;
	rv[3]  = m[3]*f0 + m[7]*f1 + m[11]*f2 + m[15]*f3;
}

//Multiply a 4x4 Matrix in column major order with a 3x1 column vector
void MatrixF::multiplyVec3(float* rv, float* m, float* v)
{
    float d=v[0], e=v[1], b=v[2];
    rv[0]=m[0]*d+m[4]*e+m[8]*b+m[12];
    rv[1]=m[1]*d+m[5]*e+m[9]*b+m[13];
    rv[2]=m[2]*d+m[6]*e+m[10]*b+m[14];
}

//Multiply a 4x4 Matrix in column major order with a 2x1 column vector with translation
//useful in transforming a 2d positional vector where translation is required (ex 2d position)
void MatrixF::multiplyVec2_Pos(float* rv, float* m, float* v)
{
    float d=v[0], e=v[1];
    rv[0]=m[0]*d+m[4]*e+m[12];
    rv[1]=m[1]*d+m[5]*e+m[13];
}

//Multiply a 4x4 Matrix in column major order with a 2x1 column vector without translation
//useful in transforming a 2d directional vector where no translation is needed (ex velocity)
void MatrixF::multiplyVec2_Vec(float* rv, float* m, float* v)
{
    float d=v[0], e=v[1];
    rv[0]=m[0]*d+m[4]*e;
    rv[1]=m[1]*d+m[5]*e;
}

/*
 * Rotates matrix m by angle a (in degrees) around the axis (x, y, z)
 * @param rm returns the result
 * @param m source matrix
 * @param a angle to rotate in degrees
 * @param x scale factor x
 * @param y scale factor y
 * @param z scale factor z
 */
void MatrixF::rotateM(float* rm, float* m, float a, float x, float y, float z)
{
    float r[16];
    setRotateM(r, a, x, y, z);
    multiplyMM(rm, m, r);
}

/**
 * Rotates matrix m in place by angle a (in degrees)
 * around the axis (x, y, z)
 * @param m source matrix
 * @param a angle to rotate in degrees
 * @param x scale factor x
 * @param y scale factor y
 * @param z scale factor z
 */
void MatrixF::rotateM(float* m, float a, float x, float y, float z) 
{
    float temp[32];
    setRotateM(temp, a, x, y, z);
    multiplyMM(temp+16, m, temp);
	memcpy(m, temp+16, 16*sizeof(float));
}

/**
 * Rotates matrix m by angle a (in degrees) around the axis (x, y, z)
 * @param rm returns the result
 * @param a angle to rotate in degrees in degrees
 * @param x scale factor x
 * @param y scale factor y
 * @param z scale factor z
 */
void MatrixF::setRotateM(float* rm, float a, float x, float y, float z)
{
    rm[3] = 0;
    rm[7] = 0;
    rm[11]= 0;
    rm[12]= 0;
    rm[13]= 0;
    rm[14]= 0;
    rm[15]= 1;
	a = DegToRad(a);
    float s = sinf(a);
    float c = cosf(a);
    if (1.0f == x && 0.0f == y && 0.0f == z)
	{
        rm[5] = c;   rm[10]= c;
        rm[6] = s;   rm[9] = -s;
        rm[1] = 0;   rm[2] = 0;
        rm[4] = 0;   rm[8] = 0;
        rm[0] = 1;
    }
	else if (0.0f == x && 1.0f == y && 0.0f == z)
	{
        rm[0] = c;   rm[10]= c;
        rm[8] = s;   rm[2] = -s;
        rm[1] = 0;   rm[4] = 0;
        rm[6] = 0;   rm[9] = 0;
        rm[5] = 1;
    }
	else if (0.0f == x && 0.0f == y && 1.0f == z)
	{
        rm[0] = c;   rm[5] = c;
        rm[1] = s;   rm[4] = -s;
        rm[2] = 0;   rm[6] = 0;
        rm[8] = 0;   rm[9] = 0;
        rm[10]= 1;
    }
	else
	{
		float len = sqrtf(x*x + y*y + z*z);
        if (1.0f != len) 
		{
            float recipLen = 1.0f / len;
            x *= recipLen;
            y *= recipLen;
            z *= recipLen;
        }
        float nc = 1.0f - c;
        float xy = x * y;
        float yz = y * z;
        float zx = z * x;
        float xs = x * s;
        float ys = y * s;
        float zs = z * s;
        rm[0] =  x*x*nc +  c;
        rm[4] =  xy*nc - zs;
        rm[8] =  zx*nc + ys;
        rm[1] =  xy*nc + zs;
        rm[5] =  y*y*nc +  c;
        rm[9] =  yz*nc - xs;
        rm[2] =  zx*nc - ys;
        rm[6] =  yz*nc + xs;
        rm[10] = z*z*nc +  c;
    }
}

void MatrixF::scaleM(float* m, float x, float y, float z)
{
    m[0]*=x; m[1]*=x; m[2]*=x;  m[3]*=x;
    m[4]*=y; m[5]*=y; m[6]*=y;  m[7]*=y;
    m[8]*=z; m[9]*=z; m[10]*=z; m[11]*=z;
}

void MatrixF::scaleMM(float* sm, float* m, float x, float y, float z)
{
    sm[0]=m[0]*x; sm[1]=m[1]*x; sm[2]=m[2]*x;   sm[3]=m[3]*x;
    sm[4]=m[4]*y; sm[5]=m[5]*y; sm[6]=m[6]*y;   sm[7]=m[7]*y;
    sm[8]=m[8]*z; sm[9]=m[9]*z; sm[10]=m[10]*z; sm[11]=m[11]*z;
    sm[12]=m[12]; sm[13]=m[13]; sm[14]=m[14]; sm[15]=m[15];
}

/**
* Translates matrix m by x, y, and z, putting the result in tm
* @param tm returns the result
* @param m source matrix
* @param x translation factor x
* @param y translation factor y
* @param z translation factor z
*/
void MatrixF::translateM(float* tm, float* m, float x, float y, float z) 
{
    float g=m[0],f=m[1],h=m[2],i=m[3],j=m[4],k=m[5],l=m[6],o=m[7],q=m[8],n=m[9],p=m[10],r=m[11];
    memcpy(tm, m, 12*sizeof(float));
    //tm[0]=g;tm[1]=f;tm[2]=h;tm[3]=i;tm[4]=j;tm[5]=k;tm[6]=l;tm[7]=o;tm[8]=q;tm[9]=n;tm[10]=p;tm[11]=r;
    tm[12]=g*x+j*y+q*z+m[12];
    tm[13]=f*x+k*y+n*z+m[13];
    tm[14]=h*x+l*y+p*z+m[14];
    tm[15]=i*x+o*y+r*z+m[15];

    /*
    memcpy(tm, m, 12*sizeof(float));
	//for (int i=0; i<12; i++)
	//	tm[i] = m[i];

    var g=a[0],f=a[1],h=a[2],i=a[3],j=a[4],k=a[5],l=a[6],o=a[7],m=a[8],n=a[9],p=a[10],r=a[11];c[0]=g;c[1]=f;c[2]=h;c[3]=i;c[4]=j;c[5]=k;c[6]=l;c[7]=o;c[8]=m;c[9]=n;c[10]=p;c[11]=r;c[12]=g*d+j*e+m*b+a[12];c[13]=f*d+k*e+n*b+a[13];c[14]=h*d+l*e+p*b+a[14];c[15]=i*d+o*e+r*b+a[15];return c};

	for (int i=0 ; i<4 ; i++) 
		tm[12 + i] = m[i] * x + m[4 + i] * y + m[8 + i] * z + m[12 + i];
    */
}

/**
* Translates matrix m by x, y, and z in place.
* @param m matrix
* @param x translation factor x
* @param y translation factor y
* @param z translation factor z
*/
void MatrixF::translateM(float* m, float x, float y, float z) 
{
    m[12] += m[0]*x + m[4]*y + m[8]*z;
    m[13] += m[1]*x + m[5]*y + m[9]*z;
    m[14] += m[2]*x + m[6]*y + m[10]*z;
    m[15] += m[3]*x + m[7]*y + m[11]*z;
    /*
	for (int i=0; i<4; i++) 
	{
		m[12 + i] += m[i] * x + m[4 + i] * y + m[8 + i] * z;
	}
    */
}

/*
mat4.rotate=function(a,b,c,d){var e=c[0],g=c[1];c=c[2];var f=Math.sqrt(e*e+g*g+c*c);if(!f)return null;if(f!=1){f=1/f;e*=f;g*=f;c*=f}var h=Math.sin(b),i=Math.cos(b),j=1-i;b=a[0];f=a[1];var k=a[2],l=a[3],o=a[4],m=a[5],n=a[6],p=a[7],r=a[8],s=a[9],A=a[10],B=a[11],t=e*e*j+i,u=g*e*j+c*h,v=c*e*j-g*h,w=e*g*j-c*h,x=g*g*j+i,y=c*g*j+e*h,z=e*c*j+g*h;e=g*c*j-e*h;g=c*c*j+i;if(d){if(a!=d){d[12]=a[12];d[13]=a[13];d[14]=a[14];d[15]=a[15]}}else d=a;d[0]=b*t+o*u+r*v;d[1]=f*t+m*u+s*v;d[2]=k*t+n*u+A*v;d[3]=l*t+p*u+B*
    v;d[4]=b*w+o*x+r*y;d[5]=f*w+m*x+s*y;d[6]=k*w+n*x+A*y;d[7]=l*w+p*x+B*y;d[8]=b*z+o*e+r*g;d[9]=f*z+m*e+s*g;d[10]=k*z+n*e+A*g;d[11]=l*z+p*e+B*g;return d};
*/

void MatrixF::rotateMX(float* m, float a)
{
    a = DegToRad(a);
    float d=sinf(a), b=cosf(a);
    float e=m[4], g=m[5], f=m[6], h=m[7], i=m[8], j=m[9], k=m[10], l=m[11];
    m[4]=e*b+i*d;   m[5]=g*b+j*d;   m[6]=f*b+k*d;   m[7]=h*b+l*d;
    m[8]=e*-d+i*b;  m[9]=g*-d+j*b;  m[10]=f*-d+k*b; m[11]=h*-d+l*b;
}

void MatrixF::rotateMY(float* m, float a)
{
    a = DegToRad(a);
    float d=sinf(a), b=cosf(a);
    float e=m[0], g=m[1], f=m[2], h=m[3], i=m[8], j=m[9], k=m[10], l=m[11];
    m[0]=e*b+i*-d;  m[1]=g*b+j*-d;  m[2]=f*b+k*-d;  m[3]=h*b+l*-d;
    m[8]=e*d+i*b;   m[9]=g*d+j*b;   m[10]=f*d+k*b;  m[11]=h*d+l*b;
}

void MatrixF::rotateMZ(float* m, float a)
{
    a = DegToRad(a);
    float d=sinf(a), b=cosf(a);
    float s = -d;
    float e=m[0], g=m[1], f=m[2], h=m[3], i=m[4], j=m[5], k=m[6], l=m[7];
    m[0]=e*b+i*d;   m[1]=g*b+j*d;   m[2]=f*b+k*d;   m[3]=h*b+l*d;
    m[4]=e*s+i*b;  m[5]=g*s+j*b;  m[6]=f*s+k*b;  m[7]=h*s+l*b;
}

void MatrixF::rotateMZ_Radians(float* m, float a)
{
    float d=sinf(a), b=cosf(a);
    float s = -d;
    float e=m[0], g=m[1], f=m[2], h=m[3], i=m[4], j=m[5], k=m[6], l=m[7];
    m[0]=e*b+i*d;   m[1]=g*b+j*d;   m[2]=f*b+k*d;   m[3]=h*b+l*d;
    m[4]=e*s+i*b;  m[5]=g*s+j*b;  m[6]=f*s+k*b;  m[7]=h*s+l*b;
}

float MatrixF::determinantM(float* m)
{
    float b=m[0],c=m[1],d=m[2],e=m[3],g=m[4],f=m[5],h=m[6],i=m[7],j=m[8],k=m[9],l=m[10],o=m[11],q=m[12],n=m[13],p=m[14],a=m[15];
    return q*k*h*e-j*n*h*e-q*f*l*e+g*n*l*e+j*f*p*e-g*k*p*e-q*k*d*i+j*n*d*i+q*c*l*i-b*n*l*i-j*c*p*i+b*k*p*i+q*f*d*o-g*n*d*o-q*c*h*o+b*n*h*o+g*c*p*o-b*f*p*o-j*f*d*a+g*k*d*a+j*c*h*a-b*k*h*a-g*c*l*a+b*f*l*a;
}

/*
 mat4.inverse=function(a,b){b||(b=a);var c=a[0],d=a[1],e=a[2],g=a[3],f=a[4],h=a[5],i=a[6],j=a[7],k=a[8],l=a[9],o=a[10],m=a[11],n=a[12],p=a[13],r=a[14],s=a[15],A=c*h-d*f,B=c*i-e*f,t=c*j-g*f,u=d*i-e*h,v=d*j-g*h,w=e*j-g*i,x=k*p-l*n,y=k*r-o*n,z=k*s-m*n,C=l*r-o*p,D=l*s-m*p,E=o*s-m*r,q=1/(A*E-B*D+t*C+u*z-v*y+w*x);b[0]=(h*E-i*D+j*C)*q;b[1]=(-d*E+e*D-g*C)*q;b[2]=(p*w-r*v+s*u)*q;b[3]=(-l*w+o*v-m*u)*q;b[4]=(-f*E+i*z-j*y)*q;b[5]=(c*E-e*z+g*y)*q;b[6]=(-n*w+r*t-s*B)*q;b[7]=(k*w-o*t+m*B)*q;b[8]=(f*D-h*z+j*x)*q;
 b[9]=(-c*D+d*z-g*x)*q;b[10]=(n*v-p*t+s*A)*q;b[11]=(-k*v+l*t-m*A)*q;b[12]=(-f*C+h*y-i*x)*q;b[13]=(c*C-d*y+e*x)*q;b[14]=(-n*u+p*B-r*A)*q;b[15]=(k*u-l*B+o*A)*q;return b};mat4.toRotationMat=function(a,b){b||(b=mat4.create());b[0]=a[0];b[1]=a[1];b[2]=a[2];b[3]=a[3];b[4]=a[4];b[5]=a[5];b[6]=a[6];b[7]=a[7];b[8]=a[8];b[9]=a[9];b[10]=a[10];b[11]=a[11];b[12]=0;b[13]=0;b[14]=0;b[15]=1;return b};
 */


/**
 * Transposes a 4 x 4 matrix.
 *
 * @param mTrans the array that holds the output transposed matrix
 * @param m the input array
 */
void MatrixF::transposeM(float* mTrans, float* m)
{
    mTrans[0]=m[0];     mTrans[1]=m[4];     mTrans[2]=m[8];     mTrans[3]=m[12];
    mTrans[4]=m[1];     mTrans[5]=m[5];     mTrans[6]=m[9];     mTrans[7]=m[13];
    mTrans[8]=m[2];     mTrans[9]=m[6];     mTrans[10]=m[10];   mTrans[11]=m[14];
    mTrans[12]=m[3];    mTrans[13]=m[7];    mTrans[14]=m[11];   mTrans[15]=m[15];
    /*
	for (int i = 0; i < 4; i++) 
	{
		int mBase = i * 4;
		mTrans[i] = m[mBase];
		mTrans[i + 4] = m[mBase + 1];
		mTrans[i + 8] = m[mBase + 2];
		mTrans[i + 12] = m[mBase + 3];
	}
    */
}

/**
 * Inverts a 4 x 4 matrix.
 *
 * @param mInv the array that holds the output inverted matrix
 * @param m the input array
 * @return true if the matrix could be inverted, false if it could not.
 */
bool MatrixF::invertM(float* mInv, float* m)
{
    // Invert a 4 x 4 matrix using Cramer's Rule

    // array of transpose source matrix
    float src[16];

    // transpose matrix
    transposeM(src, m);

    // temp array for pairs
    float tmp[12];

    // calculate pairs for first 8 elements (cofactors)
    tmp[0] = src[10] * src[15];
    tmp[1] = src[11] * src[14];
    tmp[2] = src[9] * src[15];
    tmp[3] = src[11] * src[13];
    tmp[4] = src[9] * src[14];
    tmp[5] = src[10] * src[13];
    tmp[6] = src[8] * src[15];
    tmp[7] = src[11] * src[12];
    tmp[8] = src[8] * src[14];
    tmp[9] = src[10] * src[12];
    tmp[10] = src[8] * src[13];
    tmp[11] = src[9] * src[12];

    // Holds the destination matrix while we're building it up.
    float dst[16];

    // calculate first 8 elements (cofactors)
    dst[0] = tmp[0] * src[5] + tmp[3] * src[6] + tmp[4] * src[7];
    dst[0] -= tmp[1] * src[5] + tmp[2] * src[6] + tmp[5] * src[7];
    dst[1] = tmp[1] * src[4] + tmp[6] * src[6] + tmp[9] * src[7];
    dst[1] -= tmp[0] * src[4] + tmp[7] * src[6] + tmp[8] * src[7];
    dst[2] = tmp[2] * src[4] + tmp[7] * src[5] + tmp[10] * src[7];
    dst[2] -= tmp[3] * src[4] + tmp[6] * src[5] + tmp[11] * src[7];
    dst[3] = tmp[5] * src[4] + tmp[8] * src[5] + tmp[11] * src[6];
    dst[3] -= tmp[4] * src[4] + tmp[9] * src[5] + tmp[10] * src[6];
    dst[4] = tmp[1] * src[1] + tmp[2] * src[2] + tmp[5] * src[3];
    dst[4] -= tmp[0] * src[1] + tmp[3] * src[2] + tmp[4] * src[3];
    dst[5] = tmp[0] * src[0] + tmp[7] * src[2] + tmp[8] * src[3];
    dst[5] -= tmp[1] * src[0] + tmp[6] * src[2] + tmp[9] * src[3];
    dst[6] = tmp[3] * src[0] + tmp[6] * src[1] + tmp[11] * src[3];
    dst[6] -= tmp[2] * src[0] + tmp[7] * src[1] + tmp[10] * src[3];
    dst[7] = tmp[4] * src[0] + tmp[9] * src[1] + tmp[10] * src[2];
    dst[7] -= tmp[5] * src[0] + tmp[8] * src[1] + tmp[11] * src[2];

    // calculate pairs for second 8 elements (cofactors)
    tmp[0] = src[2] * src[7];
    tmp[1] = src[3] * src[6];
    tmp[2] = src[1] * src[7];
    tmp[3] = src[3] * src[5];
    tmp[4] = src[1] * src[6];
    tmp[5] = src[2] * src[5];
    tmp[6] = src[0] * src[7];
    tmp[7] = src[3] * src[4];
    tmp[8] = src[0] * src[6];
    tmp[9] = src[2] * src[4];
    tmp[10] = src[0] * src[5];
    tmp[11] = src[1] * src[4];

    // calculate second 8 elements (cofactors)
    dst[8] = tmp[0] * src[13] + tmp[3] * src[14] + tmp[4] * src[15];
    dst[8] -= tmp[1] * src[13] + tmp[2] * src[14] + tmp[5] * src[15];
    dst[9] = tmp[1] * src[12] + tmp[6] * src[14] + tmp[9] * src[15];
    dst[9] -= tmp[0] * src[12] + tmp[7] * src[14] + tmp[8] * src[15];
    dst[10] = tmp[2] * src[12] + tmp[7] * src[13] + tmp[10] * src[15];
    dst[10] -= tmp[3] * src[12] + tmp[6] * src[13] + tmp[11] * src[15];
    dst[11] = tmp[5] * src[12] + tmp[8] * src[13] + tmp[11] * src[14];
    dst[11] -= tmp[4] * src[12] + tmp[9] * src[13] + tmp[10] * src[14];
    dst[12] = tmp[2] * src[10] + tmp[5] * src[11] + tmp[1] * src[9];
    dst[12] -= tmp[4] * src[11] + tmp[0] * src[9] + tmp[3] * src[10];
    dst[13] = tmp[8] * src[11] + tmp[0] * src[8] + tmp[7] * src[10];
    dst[13] -= tmp[6] * src[10] + tmp[9] * src[11] + tmp[1] * src[8];
    dst[14] = tmp[6] * src[9] + tmp[11] * src[11] + tmp[3] * src[8];
    dst[14] -= tmp[10] * src[11] + tmp[2] * src[8] + tmp[7] * src[9];
    dst[15] = tmp[10] * src[10] + tmp[4] * src[8] + tmp[9] * src[9];
    dst[15] -= tmp[8] * src[9] + tmp[11] * src[10] + tmp[5] * src[8];

    // calculate determinant
    float det = src[0] * dst[0] + src[1] * dst[1] + src[2] * dst[2] + src[3] * dst[3];
    if (det == 0.0f) return false;

    //calculate matrix inverse
    det = 1.0f / det;
    for (int j = 0; j < 16; j++)
        mInv[j] = dst[j] * det;

    return true;
}

#define MAT_SINGU_THRESHOLD 0.000001f
bool MatrixF::invertMGauss(float* mInv, float* m)
{
    float mat[16], inv[16];
    //we assume the matrix is invertible
    memcpy(inv, identityM, 16*sizeof(float));
    memcpy(mat, m, 16*sizeof(float));
    
    int i, k;
    float temp1, temp2;
    for (i=0; i<4; i++)
    {
        if (mat[5*i] == 0.0f)   //<= MAT_SINGU_THRESHOLD)
        {
            //pivot
            for (int a=i+1; a<4; a++)
            {
                if (mat[i*4+a] == 0.0f) continue;
                
                int g = i, h = a;
                temp1 = mat[g]; mat[g] = mat[h]; mat[h] = temp1;
                temp1 = inv[g]; inv[g] = inv[h]; inv[h] = temp1;
                g += 4; h +=4 ;
                temp1 = mat[g]; mat[g] = mat[h]; mat[h] = temp1;
                temp1 = inv[g]; inv[g] = inv[h]; inv[h] = temp1;
                g += 4; h +=4 ;
                temp1 = mat[g]; mat[g] = mat[h]; mat[h] = temp1;
                temp1 = inv[g]; inv[g] = inv[h]; inv[h] = temp1;
                g += 4; h +=4 ;
                temp1 = mat[g]; mat[g] = mat[h]; mat[h] = temp1;
                temp1 = inv[g]; inv[g] = inv[h]; inv[h] = temp1;
                
                break;
            }
        }
        //eliminate
        if (mat[5*i] == 0.0f) return false; //singular - no inverse
        temp1 = 1.0f/mat[5*i];
        mat[i] *= temp1;
        mat[i+4] *= temp1;
        mat[i+8] *= temp1;
        mat[i+12] *= temp1;
        inv[i] *= temp1;
        inv[i+4] *= temp1;
        inv[i+8] *= temp1;
        inv[i+12] *= temp1;
        
        for (k=0; k<4; k++)
        {
            if (k==i) continue;
            temp2 = mat[i*4+k];
            mat[k] -= (mat[i] * temp2);
            mat[k+4] -= (mat[i+4] * temp2);
            mat[k+8] -= (mat[i+8] * temp2);
            mat[k+12] -= (mat[i+12] * temp2);
            inv[k] -= (inv[i] * temp2);
            inv[k+4] -= (inv[i+4] * temp2);
            inv[k+8] -= (inv[i+8] * temp2);
            inv[k+12] -= (inv[i+12] * temp2);
        }
    }
    memcpy(mInv, inv, 16*sizeof(float));
    return true;
}

void MatrixF::multiplyVec3Normalize(float* rv, float* m, float* v)
{
    float d=v[0], e=v[1], b=v[2];
    float x = m[0]*d+m[4]*e+m[8]*b+m[12];
    float y = m[1]*d+m[5]*e+m[9]*b+m[13];
    float z = m[2]*d+m[6]*e+m[10]*b+m[14];
    float len = sqrtf(x*x + y*y + z*z);
    if (len != 0.0f) len = 1.0f/len;
    rv[0] = x*len;
    rv[1] = y*len;
    rv[2] = z*len;
}

void MatrixF::frustum(float left, float right, float bottom, float top, float zNear, float zFar, MatrixF& dest)
{
    float h = 1.0f / (right - left);
    float i = 1.0f / (top - bottom);
    float j = 1.0f / (zFar - zNear);
    
    float* f = dest.f;
    memset(f, 0, 16*sizeof(float));
    f[0] = zNear*2.0f*h;
    f[5] = zNear*2.0f*i;
    f[8] = (right + left)*h;
    f[9] = (top + bottom)*i;
    f[10] = -(zFar + zNear)*j;
    f[11] = -1;
    f[14] = -(zFar*zNear*2.0f)*j;
}

void MatrixF::perspective(float fovy, float ratio, float zNear, float zFar, MatrixF& dest)
{
    float size = zNear * tanf(fovy * PIBy180);
    float size_by_ratio = size/ratio;
    MatrixF::frustum(-size, size, -size_by_ratio, size_by_ratio, zNear, zFar, dest);
}

void MatrixF::ortho(float left, float right, float bottom, float top, float zNear, float zFar, MatrixF& dest)
{
    float h = 1.0f / (right - left);
    float i = 1.0f / (top - bottom);
    float j = 1.0f / (zFar - zNear);

    float* f = dest.f;
    memset(f, 0, 16*sizeof(float));
    f[0] = 2.0f*h;
    f[5] = 2.0f*i;
    f[10] = -2.0f*j;
    f[12] = -(right + left)*h;
    f[13] = -(top + bottom)*i;
    f[14] = -(zFar + zNear)*j;
    f[15] = 1.0f;
};

/*
 * Copyright 2014 Kinetic Bytes
 */

#ifndef __MATRIXF_H_
#define __MATRIXF_H_

#include <Math.h>
#include <memory.h>
#include <stack>

#define PI 3.1415926535898f
#define PI_2 6.28318530718f
#define PI_by_2 1.5707963267949f
#define PIBy180 0.01745329252f
#define _180ByPI 57.29577951f

#define DegToRad(x) ((x)*PIBy180)
#define RadToDeg(x) ((x)*_180ByPI)

//////////////////////////////////
//4x4 Matrix

class MatrixF
{
public:
	static void setIdentityM(float* f);
	static void multiplyMM(float* rm, float* m1, float* m2);
	static void multiplyMV(float* rv, float* m, float* v);
    static void multiplyVec3(float* rv, float* m, float* v);
    static void multiplyVec2_Pos(float* rv, float* m, float* v);
    static void multiplyVec2_Vec(float* rv, float* m, float* v);
	static void rotateM(float* rm, float* m, float a, float x, float y, float z);
	static void rotateM(float* m, float a, float x, float y, float z);
	static void setRotateM(float* rm, float a, float x, float y, float z);
    static void rotateMX(float* m, float a);
    static void rotateMY(float* m, float a);
    static void rotateMZ(float* m, float a);
    static void rotateMZ_Radians(float* m, float a);
    static void scaleM(float* m, float x, float y, float z);
    static void scaleMM(float* sm, float* m, float x, float y, float z);
	static void translateM(float* m, float x, float y, float z);
	static void translateM(float* tm, float* m, float x, float y, float z);
    static float determinantM(float* m);
	static bool invertM(float* mInv, float* m);
    static bool invertMGauss(float* mInv, float* m);
	static void transposeM(float* mTrans, float* m);
    static void multiplyVec3Normalize(float* rv, float* m, float* v);
    
    static void frustum(float left, float right, float bottom, float top, float zNear, float zFar, MatrixF& dest);
    static void ortho(float left, float right, float bottom, float top, float zNear, float zFar, MatrixF& dest);
    static void perspective(float fovy, float ratio, float zNear, float zFar, MatrixF& dest);
    
    //matrix stack manipulation
    static void pushModelView()
    {
        mvMatrixStack.push(mvMatrix);
    }
    static void popModelView()
    {
        mvMatrix = mvMatrixStack.top();
        mvMatrixStack.pop();
    }
    static void pushProjection()
    {
        pMatrixStack.push(pMatrix);
    }
    static void popProjection()
    {
        pMatrix = pMatrixStack.top();
        pMatrixStack.pop();
    }
    
	MatrixF() { memset(f, 0, 16*sizeof(float)); }
    MatrixF(const MatrixF& src) { memcpy(f, src.f, 16*sizeof(float)); } //copy constructor
    MatrixF& operator=(const MatrixF& src) { memcpy(f, src.f, 16*sizeof(float)); return *this; }  //assignment operator overload
    void copyOf(const MatrixF& src) { memcpy(f, src.f, 16*sizeof(float)); }

	//data (holds a 4x4 matrix as an array of floats
public:
	float f[16];
    static MatrixF orthoMatrix;
    static MatrixF perspectiveMatrix;
    static MatrixF mvMatrix;       //Model View
    static MatrixF pMatrix;        //Projection
    static MatrixF mvpMatrix;      //Model View Projection
    
private:
    static std::stack<MatrixF> mvMatrixStack;
    static std::stack<MatrixF> pMatrixStack;
};

////////////////////////////
//GLUT Helpers
void gluLookAt(float eyeX, float eyeY, float eyeZ,
        float centerX, float centerY, float centerZ, 
		float upX, float upY, float upZ);

#endif	//__MATRIXF_H_
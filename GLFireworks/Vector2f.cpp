/*
 * Copyright 2014 Kinetic Bytes
 */

#include "Vector2f.h"

////////////////////////////////////////////////////
#if defined(WIN32) && defined(_DEBUG)
	#define _CRTDBG_MAP_ALLOC
	#include <stdlib.h>
	#include <crtdbg.h>
	#define DEBUG_NEW new( _NORMAL_BLOCK, __FILE__, __LINE__ )
	#define new DEBUG_NEW
#endif
////////////////////////////////////////////////////

Vector2f& Vector2f::Normalize()					
{
	float f = x*x + y*y;
	if (f > 0.0f)
	{
		f = 1.0f/sqrtf(f);
		x*=f; y*=f;		//convert to unitvector
	}
	return *this;
}

Vector2f Vector2f::UnitNormal()
{
    float f = sqrtf(x*x + y*y);
    if (f > 0.0f)
    {
        f = 1.0f/f;
        return Vector2f(-y*f, x*f);
    }
    else
        return Vector2f();
}

void Vector2f::UnitNormal(Vector2f& outVec)
{
	float f = sqrtf(x*x + y*y);
	if (f > 0.0f)
	{
        f = 1.0f/f;
        outVec.x = -y*f; outVec.y = x*f;
	}
}

void Vector2f::UnitNormalInPlace()
{
    float f = sqrtf(x*x + y*y);
    if (f > 0.0f)
    {
        f = 1.0f/f;
        float t = -y*f; y = x*f; x = t;
    }
}

Vector2f Vector2f::UnitVector()
{
	float f = sqrtf(x*x + y*y);
	return (f > 0.0f) ? Vector2f(x/f, y/f) : Vector2f();
}

Vector2f& Vector2f::UnitVector(Vector2f& outVec)				
{
	float f = sqrtf(x*x + y*y);
	if (f > 0.0f) {outVec.x = x/f; outVec.y = y/f;}
	else outVec.Clear();
	return outVec;
}

Vector2f Vector2f::Rotate2D(float ang)
{
	float sa = sinf(ang), ca = cosf(ang);
	float xrot = x*ca - y*sa;
	float yrot = x*sa + y*ca;
	return Vector2f(xrot, yrot);
}

Vector2f Vector2f::Rotate2D(float ang, Vector2f& out) const
{
    float sa = sinf(ang), ca = cosf(ang);
    out.x = x*ca - y*sa;
    out.y = x*sa + y*ca;
    return out;
}

Vector2f& Vector2f::Rotate2DInPlace(float ang)
{
	float sa = sinf(ang), ca = cosf(ang);
	float xrot = x*ca - y*sa;
	float yrot = x*sa + y*ca;
	x = xrot; y = yrot;
	return *this;
}

Vector2f& Vector2f::Rotate2DInPlace(float sina, float cosa)
{
    float xrot = x*cosa - y*sina;
    float yrot = x*sina + y*cosa;
    x = xrot; y = yrot;
    return *this;
}

void Vector2f::FastShrink(float factor)
{
	float f = InvSqrt(x*x + y*y);
	x -= (x*factor*f);
	y -= (y*factor*f);
}

float Vector2f::InvSqrt(float x) 
{ 
    float xhalf = 0.5f*x;
    int i = *(int*)&x; // get bits for floating value 
    i = 0x5f3759df - (i>>1); // gives initial guess y0 
    x = *(float*)&i;   // convert bits back to float 
    x = x*(1.5f-xhalf*x*x); // Newton step, repeating increases accuracy 
    return x; 
}

Vector2f& Vector2f::ApplyHyperbola(float threshold)
{
    float dist = sqrtf(threshold*threshold + Length2())-threshold;
    Normalize().mulInPlace(dist);
    return *this;
}

/********************************
 intersection point of two line segments given their end points P1, P2 and P3, P4
 ---------------------------------------------------------------------------------------------------------
 Any point Pa on segment P1P2 is Pa = P1 + ta(P1P2)	where 0 <= ta <= 1
 Any point Pb on segment P3P4 is Pb = P3 + tb(P3P4)	where 0 <= tb <= 1
 At intersection, Pa = Pb, so
 
 P1 + ta(P1P2) = P3 + tb(P3P4)
 or:  P1x + ta(P1P2x) = P3x + tb(P3P4x)
 P1y + ta(P1P2y) = P3y + tb(P3P4y)
 ta = (P3x + tb(P3P4x) - P1x) / P1P2x	.....from first equ.
 
 so,
 tb = [(P1y-P3y)P1P2x + (P3x-P1x)P1P2y]/[P3P4y.P1P2x - P3P4x.P1P2y]
 ta = [(P1y-P3y)P3P4x + (P3x-P1x)P3P4y]/[P3P4y.P1P2x - P3P4x.P1P2y]
 
 or,
 tb = (P1P2 X P1P3) / (P3P4 X P1P2);
 ta = (P3P4 X P1P3) / (P3P4 X P1P2);
 
 if the denom [P3P4 X P1P2] is 0, the segments are parallel
 if 0 <= ta <= 1 then the intersection point lies on segment P1P2
 if 0 <= tb <= 1 then the intersection point lies on segment P3P4
 ********************************/
int Vector2f::LinesIntersect(const Vector2f& P1, const Vector2f& P2, const Vector2f& P3, const Vector2f& P4, Vector2f& ptIntersection)
{
	//return 0 if intersecting point lies on both segments
	//return 1 if intersecting point lies on segment 1
	//return 2 if intersecting point lies on segment 2
	//return 3 if intersecting point lies outside both segments
	//return -1 if lines are parallel
	Vector2f P1P2 = P2.sub(P1);
	Vector2f P3P4 = P4.sub(P3);
	Vector2f P1P3 = P3.sub(P1);
	float denom = P3P4.Cross(P1P2);
	if(denom == 0.0f) return -1;	//parallel
	float ta = P3P4.Cross(P1P3)/denom;
	float tb = P1P2.Cross(P1P3)/denom;
	bool ba = (ta >= 0.0f && ta <= 1.0f);
	bool bb = (tb >= 0.0f && tb <= 1.0f);
	//calculate point of intersection
	ptIntersection = P1.add(P1P2.mul(ta));
	return (ba && bb) ? 0 : (ba ? 1 : (bb ? 2 : 3));
}

//use this if point of intersection is not needed
int Vector2f::LinesIntersect(const Vector2f& P1, const Vector2f& P2, const Vector2f& P3, const Vector2f& P4)
{
    //return 0 if intersecting point lies on both segments
    //return 1 if intersecting point lies on segment 1
    //return 2 if intersecting point lies on segment 2
    //return 3 if intersecting point lies outside both segments
    //return -1 if lines are parallel
    Vector2f P1P2 = P2.sub(P1);
    Vector2f P3P4 = P4.sub(P3);
    Vector2f P1P3 = P3.sub(P1);
    float denom = P3P4.Cross(P1P2);
    if (denom == 0.0f) return -1;	//parallel
    float one_by_denom = 1.0f/denom;
    float ta = P3P4.Cross(P1P3)*one_by_denom;
    float tb = P1P2.Cross(P1P3)*one_by_denom;
    bool ba = (ta >= 0.0f && ta <= 1.0f);
    bool bb = (tb >= 0.0f && tb <= 1.0f);
    return (ba && bb) ? 0 : (ba ? 1 : (bb ? 2 : 3));
}

bool Vector2f::PointInTriangle(const Vector2f& P, const Vector2f& A, const Vector2f& B, const Vector2f& C)
{
    float ax, ay, bx, by, cx, cy, apx, apy, bpx, bpy, cpx, cpy;
    float cCROSSap, bCROSScp, aCROSSbp;
    
    ax = C.x - B.x;  ay = C.y - B.y;
    bx = A.x - C.x;  by = A.y - C.y;
    cx = B.x - A.x;  cy = B.y - A.y;
    apx= P.x - A.x;  apy= P.y - A.y;
    bpx= P.x - B.x;  bpy= P.y - B.y;
    cpx= P.x - C.x;  cpy= P.y - C.y;
    
    aCROSSbp = ax*bpy - ay*bpx;
    cCROSSap = cx*apy - cy*apx;
    bCROSScp = bx*cpy - by*cpx;
    
    return ((aCROSSbp <= 0.0f) && (bCROSScp <= 0.0f) && (cCROSSap <= 0.0f));    //using anticlockwise points
}

float Vector2f::Distance2(Vector2f& to)
{
    float xx = x-to.x, yy = y-to.y;
    return xx*xx + yy*yy;
}

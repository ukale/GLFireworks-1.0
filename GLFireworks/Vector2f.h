/*
 * Copyright 2014 Kinetic Bytes
 */

#ifndef __VECTOR2F_H_
#define __VECTOR2F_H_

#include <Math.h>

#define Round_F5(x) (((int)((x)*100000.0f))*0.00001f)
#define Round_F4(x) (((int)((x)*10000.0f))*0.0001f)
#define Round_F3(x) (((int)((x)*1000.0f))*0.001f)
#define Round_F2(x) (((int)((x)*100.0f))*0.01f)
#define Round_F1(x) (((int)((x)*10.0f))*0.1f)

class Vector2f 
{
public:
	Vector2f() { x = y = 0.0f; }
	//Vector2f(Vector2f& src) { x = src.x; y = src.y; }
	Vector2f(float xx, float yy) { x = xx; y = yy; }
	
	void set(float xx, float yy) { x = xx; y = yy; }
	void set(const Vector2f& src) { x = src.x; y = src.y; }
	void Clear() { x = y = 0.0f; }

	Vector2f Clone() { return Vector2f(x, y); }
	void CopyOf(const Vector2f& src) { x = src.x; y = src.y; }
	Vector2f& operator=(const Vector2f& src) { x=src.x; y=src.y; return *this; }
	
	float Dot(const Vector2f& v) const { return x * v.x + y * v.y; }
    
    //assuming i X j = k using right hand thumb rule
    float Cross(const Vector2f& v) const { return x * v.y - y * v.x; }
	//float Cross(const Vector2f& v) const { return y * v.x - x * v.y; }    //this appears to be incorrect

    Vector2f CrossVec3Z(float z) const
    {
        //this is a simplified version of A(vec3) X B(vec3) where A = {0,0,z} and B = {x,y,0}
        return Vector2f(-z*y, z*x);
    }
    Vector2f& CrossVec3Z(float z, Vector2f& out) const
    {
        //this is a simplified version of A(vec3) X B(vec3) where A = {0,0,z} and B = {x,y,0}
        out.x = -z*y;
        out.y = z*x;
        return out;
    }
    
    void round1()
    {
        x = ((int)(x*10.0f))*0.1f;
        y = ((int)(y*10.0f))*0.1f;
    }
    void round2()
    {
        x = ((int)(x*100.0f))*0.01f;
        y = ((int)(y*100.0f))*0.01f;
    }
    void round3()
    {
        x = ((int)(x*1000.0f))*0.001f;
        y = ((int)(y*1000.0f))*0.001f;
    }
    void round4()
    {
        x = ((int)(x*10000.0f))*0.0001f;
        y = ((int)(y*10000.0f))*0.0001f;
    }

	Vector2f sub(const Vector2f& src) const { return Vector2f(x-src.x, y-src.y); }
	Vector2f& sub(const Vector2f& src, Vector2f& out) const
	{
		out.x = x-src.x;
		out.y = y-src.y;
		return out;
	}

	Vector2f add(const Vector2f& src) const { return Vector2f(x+src.x, y+src.y); }
	Vector2f& add(const Vector2f& src, Vector2f& out) const
	{
		out.x = x+src.x;
		out.y = y+src.y;
		return out;
	}

	Vector2f mul(float f) const { return Vector2f(x*f, y*f); }
	Vector2f& mul(float f, Vector2f& out) const
	{
		out.x = x*f;
		out.y = y*f;
		return out;
	}

	Vector2f div(float f) const { return mul(1.0f/f); }
	Vector2f& div(float f, Vector2f& out) const
	{
		out.x = x/f;
		out.y = y/f;
		return out;
	}

	Vector2f& subInPlace(const Vector2f& src) { x -= src.x; y -= src.y; return *this;}
    Vector2f& subInPlace(float xx, float yy) { x -= xx; y -= yy; return *this;}
    Vector2f& subInPlace(float f) { x -= f; y -= f; return *this;}
	Vector2f& addInPlace(const Vector2f& src) { x += src.x; y += src.y; return *this;}
    Vector2f& addInPlace(float f) { x += f; y += f; return *this;}
	Vector2f& mulInPlace(float f) { x *= f; y *= f; return *this;}
	Vector2f& divInPlace(float f) { x /= f; y /= f; return *this;}

	bool Equals(const Vector2f& src) { return (x == src.x && y == src.y); }
	Vector2f& Normalize();
	void FastNormalize() { float f = InvSqrt(x*x + y*y); x*=f; y*=f; }
	
	float Length2() { return x*x + y*y; }
	float Length() { return sqrtf(x*x + y*y); }
    float Distance(Vector2f& to) { return sqrtf(Distance2(to)); }
    float Distance2(Vector2f& to);
	
	bool IsZeroVector() { return (Length2() == 0.0f); }
	Vector2f UnitNormal();
	void UnitNormal(Vector2f& outVec);
    void UnitNormalInPlace();

	Vector2f FastUnitNormal()		
	{
		float f = InvSqrt(x*x + y*y); 
		return Vector2f(-y*f, x*f);
	}

	Vector2f UnitVector();
	Vector2f& UnitVector(Vector2f& outVec);

	Vector2f FastUnitVector()				
	{
		float f = InvSqrt(x*x + y*y); 
		return Vector2f(x*f, y*f);
	}

	Vector2f Rotate2D(float ang);
    Vector2f Rotate2D(float ang, Vector2f& out) const;
	Vector2f& Rotate2DInPlace(float ang);
    Vector2f& Rotate2DInPlace(float sina, float cosa);
    
	Vector2f& Translate2DInPlace(float xx, float yy)
	{
		x += xx; y += yy;
		return *this;
	}
	void Shrink(float factor)
	{
		float f = factor/sqrtf(x*x + y*y);	//factor times reciprocal of length
		x -= (x*f);
		y -= (y*f);
	}
	void Shrink(float factor, float len)
	{
		float f = factor/len;
		x -= (x*f);
		y -= (y*f);
	}
    void Shrink(float factor, float len, Vector2f& dir)
    {
        //shrink along direction specified
        float f = factor/len;
        x -= (x*f*dir.x);
        y -= (y*f*dir.y);
    }
	void FastShrink(float factor);
    Vector2f& ApplyHyperbola(float threshold);

    static int LinesIntersect(const Vector2f& P1, const Vector2f& P2, const Vector2f& P3, const Vector2f& P4, Vector2f& ptIntersection);
    static int LinesIntersect(const Vector2f& P1, const Vector2f& P2, const Vector2f& P3, const Vector2f& P4);
    static bool PointInTriangle(const Vector2f& P, const Vector2f& A, const Vector2f& B, const Vector2f& C);
    
private:
	float InvSqrt(float x);

	////////////////data
public:
	float x;
	float y;
};

#endif	//__VECTOR2F_H_
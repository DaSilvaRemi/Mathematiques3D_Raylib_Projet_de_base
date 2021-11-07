/*******************************************************************************************
*
*   raylib [core] example - Basic window
*
*   Welcome to raylib!
*
*   To test examples, just press F6 and execute raylib_compile_execute script
*   Note that compiled executable is placed in the same folder as .c file
*
*   You can find all basic examples on C:\raylib\raylib\examples folder or
*   raylib official webpage: www.raylib.com
*
*   Enjoy using raylib. :)
*
*   This example has been created using raylib 1.0 (www.raylib.com)
*   raylib is licensed under an unmodified zlib/libpng license (View raylib.h for details)
*
*   Copyright (c) 2014 Ramon Santamaria (@raysan5)
*
********************************************************************************************/

#include "raylib.h"
#include <raymath.h>
#include "rlgl.h"
#include <math.h>
#include <float.h>
#include <vector>
#include <iostream>
#include "tools.h"

#if defined(PLATFORM_DESKTOP)
#define GLSL_VERSION            330
#else   // PLATFORM_RPI, PLATFORM_ANDROID, PLATFORM_WEB
#define GLSL_VERSION            100
#endif

#define EPSILON 1.e-6f


template <typename T> int sgn(T val) {
	return (T(0) < val) - (val < T(0));
}



/**
*
* 
*/
Cylindrical CartesianToCylindrical(Vector3 cart)
{
	Cylindrical cyl;
	cyl.rho = sqrtf(cart.x * cart.x + cart.z * cart.z);
	cyl.y = cart.y;

	if (cyl.rho < EPSILON)cyl.theta = 0;
	else
	{
		cyl.theta = atan2f(cart.x, cart.z);
		if (cyl.theta < 0)cyl.theta += PI * 2;
	}
	return cyl;
}

/**
*
* 
*/
Vector3 CylindricalToCartesian(Cylindrical cyl)
{
	return Vector3{ cyl.rho * sinf(cyl.theta),cyl.y,cyl.rho * cosf(cyl.theta) };
}

/**
*
* 
*/
Vector3 SphericalToCartesian(Spherical sph)
{
	return Vector3{ sph.rho * sinf(sph.phi) * sinf(sph.theta),
	sph.rho * cosf(sph.phi),
	sph.rho * sinf(sph.phi) * cosf(sph.theta) };
}

/**
*
* 
*/
bool InterSegPlane(Segment seg, Plane plane, Vector3& interPt, Vector3& interNormal) {
	Vector3 AB = Vector3Subtract(seg.pt2, seg.pt1);
	float dotABn = Vector3DotProduct(AB, plane.normal);

	if (fabs(dotABn) < EPSILON) return false;
	float t = (plane.d - Vector3DotProduct(seg.pt1, plane.normal)) / dotABn;

	if (t < 0 || t > 1) return false;
	interPt = Vector3Add(seg.pt1, Vector3Scale(AB, t));

	if (dotABn < 0) interNormal = plane.normal;
	else interNormal = Vector3Negate(plane.normal);

	return true;
}

/**
*
* 
*/
bool InterSegSphere(Segment seg, Sphere sphere, Vector3& interPt, Vector3& interNormal) {
	Vector3 AB = Vector3Subtract(seg.pt2, seg.pt1);
	Vector3 OmegaA = Vector3Subtract(seg.pt1, sphere.omega);

	
	float a = Vector3DotProduct(AB, AB);
	float b = 2 * Vector3DotProduct(AB, OmegaA);
	float c = Vector3DotProduct(OmegaA, OmegaA) - powf(sphere.rayon, 2);

	float discrimin = b * b - 4  * a * c;
	if (discrimin < 0) return false;

	float t = 0.0f;

	if (discrimin < EPSILON) {
		t = -(b / (2 * a));
		interPt = Vector3Add(seg.pt1, Vector3Scale(AB, t));
	}
	else {
		discrimin = sqrtf(discrimin);
		float t1 = (-b + discrimin) / (2 * a);
		float t2 = (-b - discrimin) / (2 * a);
		t = t1 < t2 ? t1 : t2;

		//interPt = Vector3Add(seg.pt1, Vector3Scale(AB, t));
		interPt = Vector3Scale(seg.pt1, t);
	
		interNormal = Vector3Normalize({ -AB.z, 0, AB.x });
	}

	if (t >= 0 && t <= 1) return true;
	else return false;
}

/**
*
* 
*/
bool InterSegmentInfiniteCylinder(Segment seg, Cylinder cyl, Vector3& interPt, Vector3& interNormal){
	Vector3 AB = Vector3Subtract(seg.pt2, seg.pt1);
	Vector3 PQ = Vector3Subtract(cyl.pt2, cyl.pt1);
	Vector3 PA = Vector3Subtract(seg.pt1, cyl.pt1);

	Vector3 tmp = Vector3Negate( Vector3Scale(PQ, Vector3DotProduct(AB, PQ) / Vector3DotProduct(PQ, PQ)));
	Vector3 tmpAB = Vector3Add(tmp, AB);
	
	tmp = Vector3Scale(PQ, Vector3DotProduct(Vector3Negate(AB), PQ) / Vector3DotProduct(PQ, PQ));
	Vector3 tmpPA = Vector3Subtract(PA, Vector3Divide(Vector3Scale(PQ, Vector3DotProduct(PA, PQ)), Vector3Multiply(PQ, PQ)));
	
	float a = Vector3DotProduct(tmpAB, tmpAB);
	//float b = Vector3Scale(Vector3CrossProduct(tmp, tmpPA), 2);
	float c = Vector3DotProduct(tmpPA, tmpPA) - powf(cyl.radius, 2);

	float discrimin = powf(b, 2) - 4 * a * c;

	float t = 0;
	if (discrimin < EPSILON) {
		t = -b / (2 * a);
		interPt = Vector3Add(seg.pt1, Vector3Scale(AB, t));
	}
	else {
		float t1 = (-b + sqrtf(discrimin)) / (2 * a);
		float t2 = (-b - sqrtf(discrimin)) / (2 * a);
		t = t1 < t2 ? t1 : t2;

		interPt = Vector3Scale(seg.pt1, t);
		interNormal = Vector3Normalize({ -AB.z, 0, AB.x });
	}

	if (t >= 0 && t <= 1) return true;
	else return false;
}

/**
*
* 
*/
/*void MyDrawCylinderPortion(Quaternion q, Cylinder cyl, float startTheta, float endTheta, int nSegmentsTheta, bool drawCaps, Color color) {
	if (nSegmentsTheta < 3) nSegmentsTheta = 3;

	int numVertex = nSegmentsTheta * 6;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	rlPushMatrix();
	rlTranslatef(cyl.pt1.x, cyl.pt1.y, cyl.pt1.z);

	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

	float height = Vector3Length(Vector3Subtract(cyl.pt2, cyl.pt1));
	rlScalef(cyl.radius, height, cyl.radius); // norme

	rlBegin(RL_TRIANGLES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	float deltaTheta = endTheta / nSegmentsTheta;

	if (cyl.radius > 0)
	{
		// Draw Body -------------------------------------------------------------------------------------
		for (float i = startTheta; i < endTheta; i += deltaTheta)
		{
			rlVertex3f(sinf(DEG2RAD * i) * cyl.radius, 0, cosf(DEG2RAD * i) * cyl.radius); //Bottom Left
			rlVertex3f(sinf(DEG2RAD * (i + deltaTheta)) * cyl.radius, 0, cosf(DEG2RAD * (i + deltaTheta)) * cyl.radius); //Bottom Right
			rlVertex3f(sinf(DEG2RAD * (i + deltaTheta)) * cyl.radius, height, cosf(DEG2RAD * (i + deltaTheta)) * cyl.radius); //Top Right

			rlVertex3f(sinf(DEG2RAD * i) * cyl.radius, height, cosf(DEG2RAD * i) * cyl.radius); //Top Left
			rlVertex3f(sinf(DEG2RAD * i) * cyl.radius, 0, cosf(DEG2RAD * i) * cyl.radius); //Bottom Left
			rlVertex3f(sinf(DEG2RAD * (i + deltaTheta)) * cyl.radius, height, cosf(DEG2RAD * (i + deltaTheta)) * cyl.radius); //Top Right
		}

		// Draw Cap --------------------------------------------------------------------------------------
		for (float i = startTheta; i < endTheta; i += deltaTheta)
		{
			rlVertex3f(0, height, 0);
			rlVertex3f(sinf(DEG2RAD * i) * cyl.radius, height, cosf(DEG2RAD * i) * cyl.radius);
			rlVertex3f(sinf(DEG2RAD * (i + deltaTheta)) * cyl.radius, height, cosf(DEG2RAD * (i + deltaTheta)) * cyl.radius);
		}
	}
	else
	{
		// Draw Cone -------------------------------------------------------------------------------------
		for (float i = startTheta; i < endTheta; i += deltaTheta)
		{
			rlVertex3f(0, height, 0);
			rlVertex3f(sinf(DEG2RAD * i) * nSegmentsTheta, 0, cosf(DEG2RAD * i) * cyl.radius);
			rlVertex3f(sinf(DEG2RAD * (i + deltaTheta)) * cyl.radius, 0, cosf(DEG2RAD * (i + deltaTheta)) * cyl.radius);
		}
	}

	// Draw Base -----------------------------------------------------------------------------------------
	for (float i = startTheta; i < endTheta; i += deltaTheta)
	{
		rlVertex3f(0, 0, 0);
		rlVertex3f(sinf(DEG2RAD * (i + deltaTheta)) * cyl.radius, 0, cosf(DEG2RAD * (i + deltaTheta)) * cyl.radius);
		rlVertex3f(sinf(DEG2RAD * i) * cyl.radius, 0, cosf(DEG2RAD * i) * cyl.radius);
	}
	rlEnd();
	rlPopMatrix();
}*/

void MyDrawCylinderPortion(Quaternion q, Cylinder cyl, float startTheta, float endTheta, int nSegmentsTheta, bool drawCaps, Color color) {
	if (nSegmentsTheta < 3) return;

	int numVertex = nSegmentsTheta * 6;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	rlPushMatrix();
	rlTranslatef(cyl.pt1.x, cyl.pt1.y, cyl.pt1.z);

	//ROTATION
	Vector3 AB = Vector3Subtract(cyl.pt2, cyl.pt1);
	Quaternion qVector = QuaternionFromVector3ToVector3({0, 1, 0}, Vector3Normalize(AB));
	Quaternion qMult = QuaternionMultiply(q, qVector);

	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	rlScalef(cyl.radius, Vector3Length(AB), cyl.radius); // norme

	rlBegin(RL_TRIANGLES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	float deltaTheta = (endTheta - startTheta) / nSegmentsTheta;
	float theta = startTheta;
	Vector3 tmpBottomLeft = CylindricalToCartesian(Cylindrical{1, theta, 0});

	for (int i = 0; i < nSegmentsTheta; i++)
	{
		float nextTheta = theta + deltaTheta;
		Vector3 bottomLeft = tmpBottomLeft;
		Vector3 topLeft = { bottomLeft.x, 1, bottomLeft.z};
		Vector3 bottomRight = CylindricalToCartesian({1, nextTheta, 0});
		Vector3 topRight = {bottomRight.x, 1, bottomRight.z};
	}
	
	rlEnd();
	rlPopMatrix();
}

/**
*
* 
*/
/*void MyDrawCylinderWiresPortion(Quaternion q, Cylinder cyl, float startTheta, float endTheta, int nSegmentsTheta, bool drawCaps, Color color) {
	if (nSegmentsTheta < 3) nSegmentsTheta = 3;

	int numVertex = nSegmentsTheta * 8;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	rlPushMatrix();
	rlTranslatef(cyl.pt1.x, cyl.pt1.y, cyl.pt1.z);

	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

	float height = Vector3Length(Vector3Subtract(cyl.pt2, cyl.pt1));
	rlScalef(cyl.radius, height, cyl.radius); // norme

	rlBegin(RL_LINES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	float deltaTheta = endTheta / nSegmentsTheta;

	for (float i = startTheta; i < endTheta; i += deltaTheta)
	{
		rlVertex3f(sinf(DEG2RAD * i) * cyl.radius, 0, cosf(DEG2RAD * i) * cyl.radius);
		rlVertex3f(sinf(DEG2RAD * (i + deltaTheta)) * cyl.radius, 0, cosf(DEG2RAD * (i + deltaTheta)) * cyl.radius);

		rlVertex3f(sinf(DEG2RAD * (i + deltaTheta)) * cyl.radius, 0, cosf(DEG2RAD * (i + deltaTheta)) * cyl.radius);
		rlVertex3f(sinf(DEG2RAD * (i + deltaTheta)) * cyl.radius, height, cosf(DEG2RAD * (i + deltaTheta)) * cyl.radius);

		rlVertex3f(sinf(DEG2RAD * (i + deltaTheta)) * cyl.radius, height, cosf(DEG2RAD * (i + deltaTheta)) * cyl.radius);
		rlVertex3f(sinf(DEG2RAD * i) * cyl.radius, height, cosf(DEG2RAD * i) * cyl.radius);

		rlVertex3f(sinf(DEG2RAD * i) * cyl.radius, height, cosf(DEG2RAD * i) * cyl.radius);
		rlVertex3f(sinf(DEG2RAD * i) * cyl.radius, 0, cosf(DEG2RAD * i) * cyl.radius);
	}
	rlEnd();
	rlPopMatrix();
}*/

/*void MyDrawCylinderWiresPortion(Quaternion q, Cylinder cyl, float startTheta, float endTheta, int nSegmentsTheta, bool drawCaps, Color color) {
	if (nSegmentsTheta < 3) return;

	int numVertex = nSegmentsTheta * 6;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	rlPushMatrix();
	rlTranslatef(cyl.pt1.x, cyl.pt1.y, cyl.pt1.z);

	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

	float height = Vector3Length(Vector3Subtract(cyl.pt2, cyl.pt1));
	rlScalef(cyl.radius, height, cyl.radius); // norme

	rlBegin(RL_LINES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	float deltaTheta = endTheta / nSegmentsTheta;

	for (float i = startTheta; i < endTheta; i += deltaTheta)
	{
		rlVertex3f(sinf(DEG2RAD * i) * cyl.radius, 0, cosf(DEG2RAD * i) * cyl.radius);
		rlVertex3f(sinf(DEG2RAD * (i + deltaTheta)) * cyl.radius, 0, cosf(DEG2RAD * (i + deltaTheta)) * cyl.radius);

		rlVertex3f(sinf(DEG2RAD * (i + deltaTheta)) * cyl.radius, 0, cosf(DEG2RAD * (i + deltaTheta)) * cyl.radius);
		rlVertex3f(sinf(DEG2RAD * (i + deltaTheta)) * cyl.radius, height, cosf(DEG2RAD * (i + deltaTheta)) * cyl.radius);

		rlVertex3f(sinf(DEG2RAD * (i + deltaTheta)) * cyl.radius, height, cosf(DEG2RAD * (i + deltaTheta)) * cyl.radius);
		rlVertex3f(sinf(DEG2RAD * i) * cyl.radius, height, cosf(DEG2RAD * i) * cyl.radius);

		rlVertex3f(sinf(DEG2RAD * i) * cyl.radius, height, cosf(DEG2RAD * i) * cyl.radius);
		rlVertex3f(sinf(DEG2RAD * i) * cyl.radius, 0, cosf(DEG2RAD * i) * cyl.radius);
	}
	rlEnd();
	rlPopMatrix();
}*/


/**
*
* 
*/
void MyDrawDiskPortion(Quaternion q, Vector3 center, float radius, float startTheta, float endTheta, int nSegmentsTheta, Color color)
{
	if (nSegmentsTheta < 3) nSegmentsTheta = 3;
	
	rlPushMatrix();
	rlTranslatef(center.x, center.y, center.z);
	rlScalef(radius, radius, radius);

	rlBegin(RL_LINES);
	rlColor4ub(color.r, color.g, color.b, color.a);
	
	float deltaTheta = endTheta / nSegmentsTheta;
	// Draw Base -----------------------------------------------------------------------------------------
	for (float i = startTheta; i < endTheta; i += deltaTheta)
	{
		rlVertex3f(0, 0, 0);
		rlVertex3f(sinf(DEG2RAD * (i + deltaTheta)) * radius, 0, cosf(DEG2RAD * (i + deltaTheta)) * radius);
		rlVertex3f(sinf(DEG2RAD * i) * radius, 0, cosf(DEG2RAD * i) * radius);
	}
	
	rlEnd();
	rlPopMatrix();
}

/**
*
* 
*/
void MyDrawDiskWiresPortion(Quaternion q, Vector3 center, float radius, float startTheta, float endTheta, int nSegmentsTheta, Color color)
{
	if (nSegmentsTheta < 3) nSegmentsTheta = 3;

	int numVertex = nSegmentsTheta * 8;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	rlBegin(RL_LINES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	float deltaTheta = endTheta / nSegmentsTheta;

	// Draw Base -----------------------------------------------------------------------------------------
	for (float i = startTheta; i < endTheta; i += deltaTheta)
	{
		rlVertex3f(0, 0, 0);
		rlVertex3f(sinf(DEG2RAD * (i + deltaTheta)) * radius, 0, cosf(DEG2RAD * (i + deltaTheta)) * radius);
		rlVertex3f(sinf(DEG2RAD * i) * radius, 0, cosf(DEG2RAD * i) * radius);
	}

	rlEnd();
	rlPopMatrix();
}

/**
*
* 
*/
void MyDrawSpherePortion(Quaternion q, Sphere sph, float startTheta, float endTheta, float startPhi, float endPhi, int nSegmentsTheta, int nSegmentsPhi, Color color) {
	if (nSegmentsTheta < 3 || nSegmentsPhi < 2) return;

	std::vector<Vector3> vertexBufferTheta(nSegmentsTheta + 1);
	std::fill(vertexBufferTheta.begin(), vertexBufferTheta.end(), Vector3{ 0,1,0 });

	int numVertex = nSegmentsTheta * nSegmentsPhi * 6;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	rlPushMatrix();

	// NOTE: Transformation is applied in inverse order (scale -> translate)
	rlTranslatef(sph.omega.x, sph.omega.y, sph.omega.z);

	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

	rlScalef(sph.rayon, sph.rayon, sph.rayon);


	rlBegin(RL_TRIANGLES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	float deltaPhi = PI / nSegmentsPhi;
	float deltaTheta = 2 * PI / nSegmentsTheta;

	float phi = 0;
	for (int i = 0; i < nSegmentsPhi; i++)
	{
		float theta = 0;
		Vector3 tmpBottomLeft = SphericalToCartesian(Spherical{ 1,theta,phi + deltaPhi });

		for (int j = 0; j < nSegmentsTheta; j++)
		{
			Vector3 topLeft = vertexBufferTheta[j];
			Vector3 bottomLeft = tmpBottomLeft;
			Vector3 topRight = vertexBufferTheta[j + 1];
			Vector3 bottomRight = SphericalToCartesian(Spherical{ 1,theta + deltaTheta,phi + deltaPhi });

			rlVertex3f(bottomLeft.x, bottomLeft.y, bottomLeft.z);
			rlVertex3f(topRight.x, topRight.y, topRight.z);
			rlVertex3f(topLeft.x, topLeft.y, topLeft.z);

			rlVertex3f(bottomLeft.x, bottomLeft.y, bottomLeft.z);
			rlVertex3f(bottomRight.x, bottomRight.y, bottomRight.z);
			rlVertex3f(topRight.x, topRight.y, topRight.z);

			theta += deltaTheta;

			vertexBufferTheta[j] = tmpBottomLeft;
			tmpBottomLeft = bottomRight;
		}
		vertexBufferTheta[vertexBufferTheta.size() - 1] = vertexBufferTheta[0];
		phi += deltaPhi;
	}
	rlEnd();
	rlPopMatrix();

}

/**
*
* 
*/
void MyDrawSphereWiresPortion(Quaternion q, Sphere sph, float startTheta, float endTheta, float startPhi, float endPhi, int nSegmentsTheta, int nSegmentsPhi, Color color) {
	if (nSegmentsTheta < 3 || nSegmentsPhi < 2) return;

	std::vector<Vector3> vertexBufferTheta(nSegmentsTheta + 1);
	std::fill(vertexBufferTheta.begin(), vertexBufferTheta.end(), Vector3{ 0,1,0 });

	int numVertex = nSegmentsTheta * nSegmentsPhi * 4;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	rlPushMatrix();
	// NOTE: Transformation is applied in inverse order (scale -> translate)
	rlTranslatef(sph.omega.x, sph.omega.y, sph.omega.z);

	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	//

	rlScalef(sph.rayon, sph.rayon, sph.rayon);

	rlBegin(RL_LINES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	float deltaPhi = endPhi / nSegmentsPhi;
	float deltaTheta = endTheta / nSegmentsTheta;

	float phi = startPhi;
	for (int i = 0; i < nSegmentsPhi; i++)
	{
		float theta = startTheta;

		for (int j = 0; j < nSegmentsTheta; j++)
		{
			Vector3 topLeft = vertexBufferTheta[j];
			Vector3 bottomLeft = SphericalToCartesian(Spherical{ 1,theta,phi + deltaPhi });
			Vector3 topRight = vertexBufferTheta[j + 1];

			rlVertex3f(topLeft.x, topLeft.y, topLeft.z);
			rlVertex3f(topRight.x, topRight.y, topRight.z);

			rlVertex3f(topLeft.x, topLeft.y, topLeft.z);
			rlVertex3f(bottomLeft.x, bottomLeft.y, bottomLeft.z);

			theta += deltaTheta;

			vertexBufferTheta[j] = bottomLeft;
		}
		vertexBufferTheta[vertexBufferTheta.size() - 1] = vertexBufferTheta[0];
		phi += deltaPhi;
	}
	rlEnd();
	rlPopMatrix();
}

/**
*
* 
*/
void MyDrawSphereEx2(Quaternion q, Sphere sph, int nSegmentsTheta, int nSegmentsPhi, Color color)
{
	MyDrawSpherePortion(q, sph, 0, 2 * PI, 0, PI, nSegmentsTheta, nSegmentsPhi, color);
}

/**
*
* 
*/
void MyDrawSphereWiresEx2(Quaternion q, Sphere sph, int nSegmentsTheta, int nSegmentsPhi, Color color)
{
	MyDrawSphereWiresPortion(q, sph, 0, 2 * PI, 0, PI, nSegmentsTheta, nSegmentsPhi, color);
}

/* Use rlVertex3f method*/
void MyDrawQuad2(Quaternion q, Vector3 center, Vector2 size, Color color) {
	//Center - Hauteur / 2
	Vector3 point1 = Vector3SubtractValue(center, size.y / 2); // -z
	//Center + Hauteur / 2
	Vector3 point2 = Vector3AddValue(center, size.y / 2); //+z
	//Center - Largeur / 2
	Vector3 point3 = Vector3SubtractValue(center, size.x / 2); //-x
	//Center + Largeur / 2
	Vector3 point4 = Vector3AddValue(center, size.x / 2); //+x

	rlPushMatrix();
	// NOTE: Transformation is applied in inverse order (scale -> translate)
	rlTranslatef(center.x, center.y, center.z);
	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	//rlScalef(size.x, size.y, center.y);

	rlBegin(RL_TRIANGLES);
	//

	rlColor4ub(color.r, color.g, color.b, color.a);

	//Left
	rlVertex3f(point4.x, center.y, point1.z);
	rlVertex3f(point3.x, center.y, point1.z);
	rlVertex3f(point4.x, center.y, point2.z);

	//Right
	rlVertex3f(point3.x, center.y, point1.z);
	rlVertex3f(point3.x, center.y, point2.z);
	rlVertex3f(point4.x, center.y, point2.z);
	rlEnd();
	rlPopMatrix();
}

/* Use rlVertex3f method*/
void MyDrawQuadWire2(Quaternion q, Vector3 center, Vector2 size, Color color) {
	//Center - Hauteur / 2
	Vector3 point1 = Vector3SubtractValue(center, size.y / 2); // -z
	//Center + Hauteur / 2
	Vector3 point2 = Vector3AddValue(center, size.y / 2); //+z
	//Center - Largeur / 2
	Vector3 point3 = Vector3SubtractValue(center, size.x / 2); //-x
	//Center + Largeur / 2
	Vector3 point4 = Vector3AddValue(center, size.x / 2); //+x

	rlPushMatrix();
	// NOTE: Transformation is applied in inverse order (scale -> translate)
	rlTranslatef(center.x, center.y, center.z);
	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

	rlBegin(RL_LINES);


	rlColor4ub(color.r, color.g, color.b, color.a);

	//Left
	rlColor4ub(color.r, color.g, color.b, color.a);
	rlVertex3f(point3.x, center.y, point1.z);
	rlVertex3f(point3.x, center.y, point2.z);
	//Right
	rlVertex3f(point4.x, center.y, point1.z);
	rlVertex3f(point4.x, center.y, point2.z);
	//Up
	rlVertex3f(point3.x, center.y, point2.z);
	rlVertex3f(point4.x, center.y, point2.z);
	//Down
	rlVertex3f(point3.x, center.y, point1.z);
	rlVertex3f(point4.x, center.y, point1.z);

	//The Intersec Line
	rlVertex3f(point3.x, center.y, point1.z);
	rlVertex3f(point4.x, center.y, point2.z);
	rlEnd();
	rlPopMatrix();
}

/*void MyDrawCylinder(Quaternion q, Cylinder cyl, int nSegmentsTheta, bool drawCaps, Color color)
{	
	if (nSegmentsTheta < 3)  nSegmentsTheta = 3;

	int numVertex = nSegmentsTheta * 6;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	rlPushMatrix();

	// NOTE: Transformation is applied in inverse order (scale -> translate)
	rlTranslatef(cyl.pt1.x, cyl.pt1.y, cyl.pt1.z);

	//ROTATION
	//Vector3 vect;
	//float angle;
	//QuaternionToAxisAngle(q, &vect, &angle);
	//rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

	float height = Vector3Length(Vector3Subtract(cyl.pt2, cyl.pt1));
	rlScalef(cyl.radius, height, cyl.radius); // norme


	rlBegin(RL_TRIANGLES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	float deltaTheta = 2 * PI / nSegmentsTheta;

	float theta = 0;
	Vector3 tmpBottomLeft = SphericalToCartesian(Spherical{ cyl.pt1.x + cyl.radius, theta,  cyl.pt1.z + cyl.radius });

	if (cyl.radius > 0) {
		for (int j = 0; j < nSegmentsTheta; j++)
		{
			//std::fill(vertexBufferTheta.begin(), vertexBufferTheta.end(), Vector3{ 0, (float) j, 0 });

			Vector3 topLeft = SphericalToCartesian(Spherical{ cyl.pt2.x + cyl.radius, height,  cyl.pt2.z + cyl.radius });
			Vector3 bottomLeft = tmpBottomLeft;
			Vector3 topRight = SphericalToCartesian(Spherical{ cyl.pt1.x + cyl.radius, height,  cyl.pt1.z + cyl.radius });
			Vector3 bottomRight = SphericalToCartesian(Spherical{ cyl.pt1.x + cyl.radius, theta + deltaTheta,  cyl.pt1.z + cyl.radius });


			rlVertex3f(bottomLeft.x, bottomLeft.y, bottomLeft.z);
			rlVertex3f(bottomRight.x, bottomRight.y, bottomRight.y);
			rlVertex3f(topRight.x, topRight.y, topRight.z);

			rlVertex3f(topLeft.x, height, topLeft.z);
			rlVertex3f(bottomLeft.x, bottomLeft.y, bottomLeft.y);
			rlVertex3f(topRight.x, height, topRight.z);

			theta += deltaTheta;
		}
	}
	else {
		theta = 0;
		for (int j = 0; j < nSegmentsTheta; j++)
		{
			Vector3 topLeft = SphericalToCartesian(Spherical{ cyl.pt2.x + cyl.radius, height,  cyl.pt2.z + cyl.radius });
			Vector3 topRight = SphericalToCartesian(Spherical{ cyl.pt2.x + cyl.radius, height,  cyl.pt2.z + cyl.radius });

			rlVertex3f(0, height, 0);
			rlVertex3f(topLeft.x, 0, topLeft.z);
			rlVertex3f(topRight.x, 0, topRight.z);

			theta += deltaTheta;
		}
	}
	theta = 0;
	// Draw Base
	for (int j = 0; j < nSegmentsTheta; j ++)
	{
		Vector3 bottomLeft = tmpBottomLeft;
		Vector3 bottomRight = SphericalToCartesian(Spherical{ cyl.pt1.x * cyl.radius, theta + deltaTheta,  cyl.pt1.z * cyl.radius });

		rlVertex3f(0, 0, 0);
		rlVertex3f(bottomRight.x, 0, bottomRight.y);
		rlVertex3f(bottomLeft.x, 0, bottomLeft.z);

		theta += deltaTheta;
	}

	rlEnd();
	rlPopMatrix();
}*/

/**
*
* 
*/
void MyDrawCylinder(Quaternion q, Cylinder cyl, int nSegmentsTheta, bool drawCaps, Color color)
{
	MyDrawCylinderPortion(q, cyl, 0.0f, 2 * PI, nSegmentsTheta, drawCaps, color);
}

/**
 *
 * 
 */
void MyDrawCylinderWires(Quaternion q, Cylinder cyl, int nSegmentsTheta, bool drawCaps, Color color) {
	MyDrawCylinderWiresPortion(q, cyl, 0, 2 * PI, nSegmentsTheta, drawCaps, color);
}

void MyDrawDisk(Quaternion q, Vector3 center, float radius, int nSegmentsTheta, Color color) {
	MyDrawDiskPortion(q, center, radius, 0, 2 * PI, nSegmentsTheta, color);
}

void MyDrawDiskWires(Quaternion q, Vector3 center, float radius, int nSegmentsTheta, Color color) {
	MyDrawDiskWiresPortion(q, center, radius, 0, 2 * PI, nSegmentsTheta, color);
}

void MyUpdateOrbitalCamera(Camera* camera, float deltaTime)
{
	static Spherical sphPos = { 20,PI / 4.f,PI / 4.f };
	static Spherical sphSpeed = { 10,.4f,.4f };
	float rhoMin = 4;
	float rhoMax = 40;

	static Vector2 prevMousePos = { 0,0 };
	Vector2 mousePos = GetMousePosition();
	Vector2 mouseVect = Vector2Subtract(mousePos, prevMousePos);
	prevMousePos = mousePos;

	Spherical sphDelta = { -GetMouseWheelMove() * sphSpeed.rho * deltaTime,
	IsMouseButtonDown(MOUSE_RIGHT_BUTTON) ? mouseVect.x * sphSpeed.theta * deltaTime  * - 1 : 0,
	IsMouseButtonDown(MOUSE_RIGHT_BUTTON) ? mouseVect.y * sphSpeed.phi * deltaTime * -1 : 0 };

	Spherical newSphPos = sphPos + sphDelta;
	newSphPos = { Clamp(newSphPos.rho,rhoMin,rhoMax),
	newSphPos.theta,
	Clamp(newSphPos.phi,PI / 100.f,.99f * PI) };

	sphPos = newSphPos;

	camera->position = SphericalToCartesian(sphPos);

}

int main(int argc, char* argv[])
{
	// Initialization
	//--------------------------------------------------------------------------------------
	float screenSizeCoef = 0.7f;
	const int screenWidth = 1920 * screenSizeCoef;
	const int screenHeight = 1080 * screenSizeCoef;

	InitWindow(screenWidth, screenHeight, "Bouncy Sphere");

	SetTargetFPS(60);

	//CAMERA
	Vector3 cameraPos = { 8.0f, 15.0f, 14.0f };
	Camera camera = { 0 };
	camera.position = cameraPos;
	camera.target = { 0.0f, 0.0f, 0.0f };
	camera.up = { 0.0f, 1.0f, 0.0f };
	camera.fovy = 45.0f;
	camera.type = CAMERA_PERSPECTIVE;
	SetCameraMode(camera, CAMERA_CUSTOM);  // Set an orbital camera mode

	//TEST CONVERSION CARTESIAN->CYLINDRICAL
	Vector3 pos = { 1,1,1 };
	Cylindrical cyl = CartesianToCylindrical(pos);
	cyl = cyl + cyl;

	/*	TEST INTERSECTIONS	*/

	Quaternion qOrient = QuaternionFromAxisAngle({ 0,0,1 }, PI * .2f);

	//TEST INTERSECTION SEGMENT PLANE
	Segment segment = { { -4, 0, -5 } , { 4, 0, 5 } };
	Plane plane = { Vector3RotateByQuaternion({ 0,1,0 }, qOrient), 2 };

	Vector3 interSectPt = { 0, 0, 0 };
	Vector3 interSecNormal = { 0, 0, 0 };
	bool planeHaveIntersec = InterSegPlane(segment, plane, interSectPt, interSecNormal);

	//TEST INTERSECTION SEGMENT PLANE
	Sphere sphere = { {0, 0, 0}, 2 };
	bool sphereHaveIntersec = InterSegSphere(segment, sphere, interSectPt, interSecNormal);

	//TEST INTERSECTION SEGMENT CYLYNDRE
	Cylinder cylinder = { {0, 0, 0}, {0, 4, 0}, 2};

	// Main game loop
	while (!WindowShouldClose())    // Detect window close button or ESC key
	{
		// Update
		//----------------------------------------------------------------------------------
		// TODO: Update your variables here
		//----------------------------------------------------------------------------------

		float deltaTime = GetFrameTime();
		float time = (float)GetTime();

		MyUpdateOrbitalCamera(&camera, deltaTime);

		// Draw
		//----------------------------------------------------------------------------------
		BeginDrawing();

		ClearBackground(RAYWHITE);

		BeginMode3D(camera);
		{
			//
			//3D REFERENTIAL
			DrawGrid(20, 1.0f);        // Draw a grid
			DrawLine3D({ 0 }, { 0,10,0 }, DARKGRAY);
			DrawSphere({ 10,0,0 }, .2f, RED);
			DrawSphere({ 0,10,0 }, .2f, GREEN);
			DrawSphere({ 0,0,10 }, .2f, BLUE);

			//INTERSEC BETWEEN SEGMENT AND PLANE
			/*MyDrawQuadWire2(qOrient, Vector3Scale(plane.normal, plane.d), {2, 2}, WHITE);
			MyDrawQuad2(qOrient, Vector3Scale(plane.normal, plane.d), { 2, 2 }, BLUE);

			if (planeHaveIntersec) {
				DrawSphere(interSectPt, .2f, DARKBROWN);
			}
			DrawLine3D(segment.pt1, segment.pt2, DARKGREEN);*/

			// INTERSEC BETWEEN SEGMENT AND SPHERE
			Quaternion qOrient = QuaternionFromAxisAngle({1,0,0}, PI * .5f);
			MyDrawSphereEx2(qOrient, sphere, 40, 20, BLUE);
			MyDrawSphereWiresEx2(qOrient, sphere, 40, 20, WHITE);

			if (sphereHaveIntersec) {
				DrawLine3D(interSectPt, interSecNormal, DARKPURPLE);
				DrawSphere(interSectPt, .2f, DARKBROWN);
			}
			DrawLine3D(segment.pt1, segment.pt2, DARKGREEN);

			/*MyDrawCylinder(qOrient, cylinder, 25, false, BLUE);
			MyDrawCylinderWires(qOrient, cylinder, 25, false, WHITE);*/
		}
		EndMode3D();

		EndDrawing();
		//----------------------------------------------------------------------------------
	}

	// De-Initialization
	//--------------------------------------------------------------------------------------  
	CloseWindow();        // Close window and OpenGL context
	//--------------------------------------------------------------------------------------

	return 0;
}
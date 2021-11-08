#include "tools.h"

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
/*bool InterSegmentInfiniteCylinder(Segment seg, Cylinder cyl, Vector3& interPt, Vector3& interNormal) {
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

	//float discrimin = powf(b, 2) - 4 * a * c;

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
}*/

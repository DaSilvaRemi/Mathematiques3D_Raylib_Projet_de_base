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

Vector3 GlobalToLocalPos(Vector3 posGlobal, Referential localRef) {
	Vector3 globalOrigin = { 0,0,0 };
	Vector3 globalVect = Vector3Subtract(Vector3Subtract(posGlobal, globalOrigin), Vector3Subtract(localRef.origin, globalOrigin));
	return GlobalToLocalVect(globalVect, localRef);
}

Vector3 GlobalToLocalVect(Vector3 vectGlobal, Referential localRef) {
	return { Vector3DotProduct(vectGlobal, localRef.i), Vector3DotProduct(vectGlobal, localRef.j) , Vector3DotProduct(vectGlobal, localRef.k) };
}

Vector3 LocalToGlobalPos(Vector3 localPos, Referential localRef) {
	Vector3 globalOrigin = { 0,0,0 };
	return Vector3Add(Vector3Subtract(localRef.origin, globalOrigin), LocalToGlobalVect(localPos, localRef));
}

Vector3 LocalToGlobalVect(Vector3 localVect, Referential localRef) {
	return Vector3Add(Vector3Add(Vector3Scale(localRef.i, localVect.x), Vector3Scale(localRef.j, localVect.y)), Vector3Scale(localRef.k, localVect.z));
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

bool InterSegQuad(Segment seg, Quad quad, Vector3& interPt, Vector3& interNormal) {
	bool isIntersec = InterSegPlane(seg, Plane(quad.referential.j, quad.referential.origin), interPt, interNormal);
	if (!isIntersec) return false;

	Vector3 localPos = GlobalToLocalPos(interPt, quad.referential);
	return ((fabsf(localPos.x) <= quad.extension.x / 2) && (fabsf(localPos.z) <= quad.extension.z / 2));
}

bool InterSegDisk(Segment seg, Disk disk, Vector3& interPt, Vector3& interNormal) {
	bool isIntersec = InterSegPlane(seg, Plane(GlobalToLocalPos(disk.referential.i, disk.referential), disk.referential.origin), interPt, interNormal);
	if (!isIntersec) return false;

	Vector3 localPos = GlobalToLocalPos(interPt, disk.referential);
	return (fabsf(localPos.x) <= disk.radius && fabsf(localPos.y) <= 1 && fabsf(localPos.z) <= disk.radius);
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

	float discrimin = b * b - 4 * a * c;
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

		interPt = Vector3Add(seg.pt1, Vector3Scale(AB, t));
		//interPt = Vector3Scale(seg.pt1, t);

		interNormal = Vector3Normalize(Vector3Subtract(interPt, sphere.omega));
	}

	if (t >= 0 && t <= 1) return true;
	else return false;
}

/**
*
*
*/
bool InterSegmentInfiniteCylinder(Segment seg, Cylinder cyl, Vector3& interPt, Vector3& interNormal) {
	Vector3 ptA = seg.pt1;
	Vector3 ptB = seg.pt2;
	Vector3 ptP = cyl.pt1;
	Vector3 ptQ = cyl.pt2;
	float r = cyl.radius;

	Vector3 AB = Vector3Subtract(ptB, ptA);
	Vector3 PQ = Vector3Subtract(ptQ, ptP);
	Vector3 PA = Vector3Subtract(ptA, ptP);

	float ABdotPQ = Vector3DotProduct(AB, PQ);
	float PAdotPQ = Vector3DotProduct(PA, PQ);
	float ABdotPA = Vector3DotProduct(AB, PA);

	float ABcarre = Vector3DotProduct(AB, AB);
	float PQcarre = Vector3DotProduct(PQ, PQ);
	float PAcarre = Vector3DotProduct(PA, PA);
	float ABPQcarre = powf(ABdotPQ, 2);

	float PAdotPQcarre = powf(PAdotPQ, 2);

	float ABPQsurPQcarre = ABdotPQ / PQcarre;
	float PAPQsurPQcarre = PAdotPQ / PQcarre;
	float PAPQcarreeSurPQcaree = PAdotPQcarre / PQcarre;

	// Soit at² + bt + c = 0
	float a = ABcarre - 2 * ABPQcarre / PQcarre + powf(ABPQsurPQcarre, 2) * PQcarre;
	float b = 2 * (ABdotPA - ABPQsurPQcarre * PAdotPQ - PAPQsurPQcarre * ABdotPQ + ABPQsurPQcarre * PAPQsurPQcarre * PQcarre);
	float c = PAcarre - 2 * PAPQcarreeSurPQcaree + powf(PAPQsurPQcarre, 2) * PQcarre - powf(r, 2);

	if (a < EPSILON) return false;

	// Discriminant
	float discriminant = powf(b, 2) - 4 * a * c;

	// test sur discriminant
	if (discriminant < 0) return false;

	// calculs racines
	float t = 0;
	float racineDiscriminant = sqrtf(discriminant);
	float t1 = (-b - racineDiscriminant) / (2 * a);
	float t2 = (-b + racineDiscriminant) / (2 * a);
	t = t1 < t2 ? t1 : t2;

	// test sur t
	if (t < 0 || t > 1) return false;

	// interPt
	interPt = Vector3Add(ptA, Vector3Scale(AB, t));

	// Valeurs initiales
	Vector3 OInter = interPt;//Vector3Subtract(interPt, {0,0,0});
	Vector3 OP = ptP;//Vector3Subtract(ptP, {0,0,0});
	Vector3 PInter = Vector3Subtract(interPt, ptP);
	float PQdotPInter = Vector3DotProduct(PQ, PInter);

	// ON = OP + (PQdotPInter / PQcarre) * PQ
	float PQdotPInterSurPQCarre = PQdotPInter / PQcarre;
	Vector3 ON = Vector3Add(OP, Vector3Scale(PQ, PQdotPInterSurPQCarre));
	Vector3 NI = Vector3Subtract(OInter, ON);

	// set interNormal
	interNormal = Vector3Normalize(NI);

	return true;
}

bool InterSegmentFiniteCylinder(Segment seg, Cylinder cyl, Vector3& interPt, Vector3& interNormal) {
	bool cylinderIsIntersec = InterSegmentInfiniteCylinder(seg, cyl, interPt, interNormal);
	bool isIntersec = false;
	
	Vector3 PInter = Vector3Subtract(interPt, cyl.pt1);
	Vector3 PQ = Vector3Subtract(cyl.pt2, cyl.pt1);
	float PInterdotPQ = Vector3DotProduct(PInter, PQ);
	float PQcarre = Vector3DotProduct(PQ, PQ);

	if (cylinderIsIntersec) {
		if (PInterdotPQ < 0 || PInterdotPQ > PQcarre) {
			Vector3 maxInterPt = { FLT_MAX };
			Vector3 tmpInterPt;
			Vector3 tmpInterNormal;


			bool isDiskIntersec = InterSegDisk(seg, { Referential(cyl.pt1), cyl.radius }, tmpInterPt, tmpInterNormal);

			if (isDiskIntersec && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(maxInterPt, seg.pt1)) {
				interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
				interNormal = { tmpInterNormal.x, tmpInterNormal.y, tmpInterNormal.z };
				maxInterPt = interPt;
				isIntersec = true;
			}

			bool isDiskIntersec2 = InterSegDisk(seg, { Referential(cyl.pt2), cyl.radius }, tmpInterPt, tmpInterNormal);
			if (isDiskIntersec2 && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(maxInterPt, seg.pt1)) {
				interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
				interNormal = { tmpInterNormal.x, tmpInterNormal.y, tmpInterNormal.z };
				isIntersec = true;
			}
		}
		else {
			isIntersec = true;
		}
	}
	else {
		Vector3 maxInterPt = { FLT_MAX };
		Vector3 tmpInterPt;
		Vector3 tmpInterNormal;

		bool isDiskIntersec = InterSegDisk(seg, { Referential(cyl.pt1), cyl.radius }, tmpInterPt, tmpInterNormal);

		if (isDiskIntersec && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(maxInterPt, seg.pt1)) {
			interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
			interNormal = { tmpInterNormal.x, tmpInterNormal.y, tmpInterNormal.z };
			maxInterPt = interPt;
			isIntersec = true;
		}

		bool isDiskIntersec2 = InterSegDisk(seg, { Referential(cyl.pt2), cyl.radius }, tmpInterPt, tmpInterNormal);
		if (isDiskIntersec2 && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(maxInterPt, seg.pt1)) {
			interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
			interNormal = { tmpInterNormal.x, tmpInterNormal.y, tmpInterNormal.z };
			isIntersec = true;
		}
	}
	
	return isIntersec;
}

bool InterSegmentCapsule(Segment seg, Capsule caps, Vector3& interPt, Vector3& interNormal) {
	bool isIntersec = false;
	bool tmpIsIntersec = false;
	interPt = { FLT_MAX };
	interNormal = { FLT_MAX };

	Vector3 up = LocalToGlobalPos({0, caps.height, 0}, caps.referential);
	Vector3 down = LocalToGlobalPos({ 0, 0, 0 }, caps.referential);

	Quaternion qUp = QuaternionFromAxisAngle({ 0, 0, 1 }, 0.5 * PI);
	Quaternion qDown = QuaternionFromAxisAngle({ 0, 0, 1 }, 1.5 * PI);
	Quaternion qIdentity = QuaternionIdentity();

	Referential ref = Referential(down);
	ref.RotateByQuaternion(caps.referential.q);
	Cylinder cylinder = Cylinder(ref, caps.radius, caps.height);
	cylinder.UpdateCylinder();

	Sphere sphereUp = { up, caps.radius };
	Sphere sphereDown = { down, caps.radius };

	Vector3 tmpInterPt;
	Vector3 tmpInterNormal;
	tmpIsIntersec = InterSegmentFiniteCylinder(seg, cylinder, tmpInterPt, tmpInterNormal);

	if (tmpIsIntersec && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(interPt, seg.pt1)) {
		interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
		interNormal = { tmpInterNormal.x, tmpInterNormal.y, tmpInterNormal.z };
		isIntersec = true;
	}

	tmpIsIntersec = InterSegSphere(seg, sphereUp, tmpInterPt, tmpInterNormal);

	if (tmpIsIntersec && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(interPt, seg.pt1)) {
		interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
		interNormal = { tmpInterNormal.x, tmpInterNormal.y, tmpInterNormal.z };
		isIntersec = true;
	}

	tmpIsIntersec = InterSegSphere(seg, sphereDown, tmpInterPt, tmpInterNormal);

	if (tmpIsIntersec && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(interPt, seg.pt1)) {
		interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
		interNormal = { tmpInterNormal.x, tmpInterNormal.y, tmpInterNormal.z };
		isIntersec = true;
	}

	//interPt = LocalToGlobalPos(interPt, caps.referential);

	return isIntersec;
}

bool IntersecSegRoundedBox(Segment seg, RoundedBox roundedBox, Vector3 &interPt, Vector3 &interNormal) {
	bool tmpIsIntersec = false;
	bool isIntersec = false;
	interPt = { FLT_MAX };
	interNormal = { FLT_MAX };
	

	// référentiel de capsules
	Quaternion qLeft = QuaternionFromAxisAngle({ 1, 0, 0 }, -PI * 0.5f);
	Capsule capsLeftBottom = { Referential(roundedBox.ref.origin, qLeft), roundedBox.radius, roundedBox.extension.z };

	Vector3 tmpInterPt;
	Vector3 tmpInterNormal;
	
	tmpIsIntersec = InterSegmentCapsule(seg, capsLeftBottom, tmpInterPt, tmpInterNormal);
	if (tmpIsIntersec && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(interPt, seg.pt1)) {
		interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
		interNormal = { tmpInterNormal.x, tmpInterNormal.y, tmpInterNormal.z };
		isIntersec = true;
	}

	Quaternion qFront = QuaternionFromAxisAngle({ 0, 0, 1 }, -PI * 0.5f);
	Capsule capsFrontBottom = { Referential(roundedBox.ref.origin, qFront), roundedBox.radius, roundedBox.extension.x };

	tmpIsIntersec = InterSegmentCapsule(seg, capsFrontBottom, tmpInterPt, tmpInterNormal);
	if (tmpIsIntersec && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(interPt, seg.pt1)) {
		interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
		interNormal = { tmpInterNormal.x, tmpInterNormal.y, tmpInterNormal.z };
		isIntersec = true;
	}

	Quaternion qUp = QuaternionIdentity();
	Capsule capsFrontLeft = { Referential(roundedBox.ref.origin, qUp), roundedBox.radius, roundedBox.extension.y };

	tmpIsIntersec = InterSegmentCapsule(seg, capsFrontLeft, tmpInterPt, tmpInterNormal);
	if (tmpIsIntersec && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(interPt, seg.pt1)) {
		interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
		interNormal = { tmpInterNormal.x, tmpInterNormal.y, tmpInterNormal.z };
		isIntersec = true;
	}

	// toutes les autres capsules
	Capsule capsFrontTop = { Referential(Vector3Add(roundedBox.ref.origin, { 0, roundedBox.extension.y, 0 }), qFront), roundedBox.radius, roundedBox.extension.x };
	
	tmpIsIntersec = InterSegmentCapsule(seg, capsFrontTop, tmpInterPt, tmpInterNormal);
	if (tmpIsIntersec && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(interPt, seg.pt1)) {
		interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
		interNormal = { tmpInterNormal.x, tmpInterNormal.y, tmpInterNormal.z };
		isIntersec = true;
	}

	Capsule capsFrontRight = { Referential(Vector3Add(roundedBox.ref.origin, {roundedBox.extension.x, 0, 0}), qUp), roundedBox.radius, roundedBox.extension.y };
	
	tmpIsIntersec = InterSegmentCapsule(seg, capsFrontRight, tmpInterPt, tmpInterNormal);
	if (tmpIsIntersec && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(interPt, seg.pt1)) {
		interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
		interNormal = { tmpInterNormal.x, tmpInterNormal.y, tmpInterNormal.z };
		isIntersec = true;
	}

	Capsule capsRightBottom = { Referential(Vector3Add(roundedBox.ref.origin, { roundedBox.extension.x, 0, 0 }), qLeft), roundedBox.radius, roundedBox.extension.z };
	
	tmpIsIntersec = InterSegmentCapsule(seg, capsRightBottom, tmpInterPt, tmpInterNormal);
	if (tmpIsIntersec && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(interPt, seg.pt1)) {
		interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
		interNormal = { tmpInterNormal.x, tmpInterNormal.y, tmpInterNormal.z };
		isIntersec = true;
	}

	Capsule capsRightTop = { Referential(Vector3Add(roundedBox.ref.origin, { roundedBox.extension.x, roundedBox.extension.y, 0 }), qLeft), roundedBox.radius, roundedBox.extension.z };
	
	tmpIsIntersec = InterSegmentCapsule(seg, capsRightTop, tmpInterPt, tmpInterNormal);
	if (tmpIsIntersec && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(interPt, seg.pt1)) {
		interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
		interNormal = { tmpInterNormal.x, tmpInterNormal.y, tmpInterNormal.z };
		isIntersec = true;
	}

	Capsule capsLeftTop = { Referential(Vector3Add(roundedBox.ref.origin, { 0, roundedBox.extension.y, 0 }), qLeft), roundedBox.radius, roundedBox.extension.z };
	
	tmpIsIntersec = InterSegmentCapsule(seg, capsLeftTop, tmpInterPt, tmpInterNormal);
	if (tmpIsIntersec && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(interPt, seg.pt1)) {
		interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
		interNormal = { tmpInterNormal.x, tmpInterNormal.y, tmpInterNormal.z };
		isIntersec = true;
	}

	Capsule capsBackBottom = { Referential(Vector3Add(roundedBox.ref.origin, {0, 0, -roundedBox.extension.z}), qFront), roundedBox.radius, roundedBox.extension.x };
	
	tmpIsIntersec = InterSegmentCapsule(seg, capsBackBottom, tmpInterPt, tmpInterNormal);
	if (tmpIsIntersec && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(interPt, seg.pt1)) {
		interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
		interNormal = { tmpInterNormal.x, tmpInterNormal.y, tmpInterNormal.z };
		isIntersec = true;
	}

	Capsule capsBackTop = { Referential(Vector3Add(roundedBox.ref.origin, { 0, roundedBox.extension.y, -roundedBox.extension.z }), qFront), roundedBox.radius, roundedBox.extension.x };
	
	tmpIsIntersec = InterSegmentCapsule(seg, capsBackTop, tmpInterPt, tmpInterNormal);
	if (tmpIsIntersec && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(interPt, seg.pt1)) {
		interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
		interNormal = { tmpInterNormal.x, tmpInterNormal.y, tmpInterNormal.z };
		isIntersec = true;
	}

	Capsule capsBackLeft = { Referential(Vector3Add(roundedBox.ref.origin, { 0, 0, -roundedBox.extension.z }), qUp), roundedBox.radius, roundedBox.extension.y };
	
	tmpIsIntersec = InterSegmentCapsule(seg, capsBackLeft, tmpInterPt, tmpInterNormal);
	if (tmpIsIntersec && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(interPt, seg.pt1)) {
		interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
		interNormal = { tmpInterNormal.x, tmpInterNormal.y, tmpInterNormal.z };
		isIntersec = true;
	}

	Capsule capsBackRight = { Referential(Vector3Add(roundedBox.ref.origin, { roundedBox.extension.x, 0, -roundedBox.extension.z }), qUp), roundedBox.radius, roundedBox.extension.y };
	
	tmpIsIntersec = InterSegmentCapsule(seg, capsBackRight, tmpInterPt, tmpInterNormal);
	if (tmpIsIntersec && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(interPt, seg.pt1)) {
		interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
		interNormal = { tmpInterNormal.x, tmpInterNormal.y, tmpInterNormal.z };
		isIntersec = true;
	}

	Quaternion qFrontQuad = QuaternionFromAxisAngle({ 1, 0, 0 }, PI * 0.5f);
	Quad quadFront = { Referential(Vector3Add(roundedBox.ref.origin, {roundedBox.extension.x / 2, roundedBox.extension.y / 2, roundedBox.radius}), qFrontQuad), {roundedBox.extension.x, roundedBox.extension.z, roundedBox.extension.y } };
	
	tmpIsIntersec = InterSegQuad(seg, quadFront, tmpInterPt, tmpInterNormal);
	if (tmpIsIntersec && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(interPt, seg.pt1)) {
		interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
		interNormal = { tmpInterNormal.x, tmpInterNormal.y, tmpInterNormal.z };
		isIntersec = true;
	}

	Quaternion qBackQuad = QuaternionFromAxisAngle({ 1, 0, 0 }, -PI * 0.5f);
	Quad quadBack = { Referential(Vector3Add(roundedBox.ref.origin, {roundedBox.extension.x / 2,  roundedBox.extension.y / 2, -(roundedBox.extension.z + roundedBox.radius)}), qBackQuad), {roundedBox.extension.x, roundedBox.extension.z, roundedBox.extension.y} };
	
	tmpIsIntersec = InterSegQuad(seg, quadBack, tmpInterPt, tmpInterNormal);
	if (tmpIsIntersec && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(interPt, seg.pt1)) {
		interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
		interNormal = { tmpInterNormal.x, tmpInterNormal.y, tmpInterNormal.z };
		isIntersec = true;
	}

	Quaternion qRightQuad = QuaternionFromAxisAngle({ 0, 0, 1 }, -PI * 0.5f);
	Quad quadRight = { Referential(Vector3Add(roundedBox.ref.origin, {roundedBox.extension.x + roundedBox.radius, roundedBox.extension.y / 2, -roundedBox.extension.z / 2}), qRightQuad), {roundedBox.extension.y, roundedBox.extension.x, roundedBox.extension.z} };
	
	tmpIsIntersec = InterSegQuad(seg, quadRight, tmpInterPt, tmpInterNormal);
	if (tmpIsIntersec && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(interPt, seg.pt1)) {
		interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
		interNormal = { tmpInterNormal.x, tmpInterNormal.y, tmpInterNormal.z };
		isIntersec = true;
	}

	Quaternion qLeftQuad = QuaternionFromAxisAngle({ 0, 0, 1 }, PI * 0.5f);
	Quad quadLeft = { Referential(Vector3Add(roundedBox.ref.origin, {-roundedBox.radius, roundedBox.extension.y / 2, -roundedBox.extension.z / 2}), qLeftQuad), {roundedBox.extension.y, roundedBox.extension.x, roundedBox.extension.z } };
	
	tmpIsIntersec = InterSegQuad(seg, quadLeft, tmpInterPt, tmpInterNormal);
	if (tmpIsIntersec && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(interPt, seg.pt1)) {
		interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
		interNormal = { tmpInterNormal.x, tmpInterNormal.y, tmpInterNormal.z };
		isIntersec = true;
	}

	Quaternion qTopQuad = QuaternionIdentity();
	Quad quadTop = { Referential(Vector3Add(roundedBox.ref.origin, {roundedBox.extension.x / 2, roundedBox.extension.y + roundedBox.radius, -roundedBox.extension.z / 2}), qTopQuad), {roundedBox.extension.x, roundedBox.extension.y, roundedBox.extension.z } };
	
	tmpIsIntersec = InterSegQuad(seg, quadTop, tmpInterPt, tmpInterNormal);
	if (tmpIsIntersec && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(interPt, seg.pt1)) {
		interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
		interNormal = { tmpInterNormal.x, tmpInterNormal.y, tmpInterNormal.z };
		isIntersec = true;
	}

	Quaternion qBottomQuad = QuaternionFromAxisAngle({ 0, 0, 1 }, PI);
	Quad quadBottom = { Referential(Vector3Add(roundedBox.ref.origin, {roundedBox.extension.x / 2, -roundedBox.radius, -roundedBox.extension.z / 2}), qBottomQuad), {roundedBox.extension.x, roundedBox.extension.y, roundedBox.extension.z } };
	tmpIsIntersec = InterSegQuad(seg, quadBottom, tmpInterPt, tmpInterNormal);
	if (tmpIsIntersec && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(interPt, seg.pt1)) {
		interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
		interNormal = { tmpInterNormal.x, tmpInterNormal.y, tmpInterNormal.z };
		isIntersec = true;
	}

	return isIntersec;
}

bool IntersecSegBox(Segment seg, Box box, Vector3& interPt, Vector3& interNormal) {
	bool tmpIsIntersec = false;
	bool isIntersec = false;
	interPt = { FLT_MAX };
	interNormal = { FLT_MAX };

	Quaternion qUp = QuaternionIdentity();
	Quaternion qLeft = QuaternionFromAxisAngle({ 1, 0, 0 }, PI * 0.5f);
	Quaternion qRight = QuaternionFromAxisAngle({ 1, 0, 0 }, PI * -0.5f);
	Quaternion qFront = QuaternionFromAxisAngle({ 0, 0, 1 }, PI * -0.5f);
	Quaternion qBack = QuaternionFromAxisAngle({ 0, 0, 1 }, PI * 0.5f);
	//Quaternion qDown = QuaternionFromAxisAngle({ 1, 0, 0 }, PI * 1);

	Quaternion qX = QuaternionFromAxisAngle({ 1, 0, 0 }, PI);
	Quaternion qY = QuaternionFromAxisAngle({ 0, 1, 0 }, 0);
	Quaternion qDown = QuaternionMultiply(qX, qY);

	Referential referentialQuadUp = Referential(GlobalToLocalPos({ -0.5f, 0.5f, 0 }, box.ref));
	Referential referentialQuadFront = Referential(GlobalToLocalPos({ 0, 0, 0 }, box.ref));
	Referential referentialQuadBack = Referential(GlobalToLocalPos({ -1, 0, 0 }, box.ref));
	Referential referentialQuadLeft = Referential(GlobalToLocalPos({ -0.5f, 0, 0.5f }, box.ref));
	Referential referentialQuadRight = Referential(GlobalToLocalPos({ -0.5f, 0, -0.5f }, box.ref));
	Referential referentialQuadDown = Referential(GlobalToLocalPos({ -0.5f, -0.5f, 0 }, box.ref));

	referentialQuadUp.RotateByQuaternion(qUp);
	referentialQuadFront.RotateByQuaternion(qFront);
	referentialQuadBack.RotateByQuaternion(qBack);
	referentialQuadLeft.RotateByQuaternion(qLeft);
	referentialQuadDown.RotateByQuaternion(qDown);

	Quad quadUp = { referentialQuadUp, box.extension };
	Quad quadFront = { referentialQuadFront, box.extension };
	Quad quadBack = { referentialQuadBack, box.extension };
	Quad quadLeft = { referentialQuadLeft, box.extension };
	Quad quadRight = { referentialQuadRight, box.extension };
	Quad quadDown = { referentialQuadDown, box.extension };

	Vector3 tmpInterPt;
	Vector3 tmpInterNormal;
	tmpIsIntersec = InterSegQuad(seg, quadUp, tmpInterPt, tmpInterNormal);

	if (tmpIsIntersec && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(interPt, seg.pt1)) {
		interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
		interNormal = { tmpInterNormal.x, tmpInterNormal.y, tmpInterNormal.z };
		isIntersec = true;
	}

	tmpIsIntersec = InterSegQuad(seg, quadFront, tmpInterPt, tmpInterNormal);

	if (tmpIsIntersec && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(interPt, seg.pt1)) {
		interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
		interNormal = { tmpInterNormal.x, tmpInterNormal.y, tmpInterNormal.z };
		isIntersec = true;
	}

	tmpIsIntersec = InterSegQuad(seg, quadBack, tmpInterPt, tmpInterNormal);

	if (tmpIsIntersec && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(interPt, seg.pt1)) {
		interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
		interNormal = { tmpInterNormal.x, tmpInterNormal.y, tmpInterNormal.z };
		isIntersec = true;
	}

	tmpIsIntersec = InterSegQuad(seg, quadLeft, tmpInterPt, tmpInterNormal);

	if (tmpIsIntersec && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(interPt, seg.pt1)) {
		interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
		interNormal = { tmpInterNormal.x, tmpInterNormal.y, tmpInterNormal.z };
		isIntersec = true;
	}

	tmpIsIntersec = InterSegQuad(seg, quadRight, tmpInterPt, tmpInterNormal);

	if (tmpIsIntersec && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(interPt, seg.pt1)) {
		interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
	}

	tmpIsIntersec = InterSegQuad(seg, quadDown, tmpInterPt, tmpInterNormal);

	if (tmpIsIntersec && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(interPt, seg.pt1)) {
		interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
		interNormal = { tmpInterNormal.x, tmpInterNormal.y, tmpInterNormal.z };
		isIntersec = true;
	}

	return isIntersec;
}

float random_float(float a, float b) {
	return a + ((((float)std::rand()) / (float)RAND_MAX) * (b - a));
}
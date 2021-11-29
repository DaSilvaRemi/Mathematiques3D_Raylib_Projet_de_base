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
	bool isIntersec = InterSegPlane(seg, Plane(quad.referential.j, LocalToGlobalPos(quad.referential.origin, quad.referential), { 0,0,0 }), interPt, interNormal);

	return isIntersec && (fabsf(interPt.x) <= quad.extension.x && fabsf(interPt.z) <= quad.extension.z);
}

bool InterSegDisk(Segment seg, Disk disk, Vector3& interPt, Vector3& interNormal) {
	bool isIntersec = InterSegPlane(seg, Plane(disk.referential.i, LocalToGlobalPos(disk.referential.origin, disk.referential), { 0,0,0 }), interPt, interNormal);

	return isIntersec && (fabsf(interPt.x) <= disk.radius && fabsf(interPt.z) <= disk.radius);
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
	bool isCylinderIntersec = InterSegmentInfiniteCylinder(seg, cyl, interPt, interNormal);
	
	Vector3 PInter = Vector3Subtract(interPt, cyl.pt1);
	Vector3 PQ = Vector3Subtract(cyl.pt2, cyl.pt1);
	float PInterdotPQ = Vector3DotProduct(PInter, PQ);
	float PQcarre = Vector3DotProduct(PQ, PQ);

	if (isCylinderIntersec) {
		if (PInterdotPQ < 0 || PInterdotPQ > PQcarre) {
			bool isDiskIntersec = InterSegDisk(seg, { cyl.pt1, cyl.radius }, interPt, interNormal);
			if (isDiskIntersec) {
				return isDiskIntersec;
			}

			bool isDiskIntersec2 = InterSegDisk(seg, { cyl.pt2, cyl.radius }, interPt, interNormal);
			if (isDiskIntersec2) {
				return isDiskIntersec2;
			}
		}
		else {
			return true;
		}
	}
	else {
		bool isDiskIntersec = InterSegDisk(seg, { cyl.pt1, cyl.radius }, interPt, interNormal);
		if (isDiskIntersec) {
			return isDiskIntersec;
		}

		bool isDiskIntersec2 = InterSegDisk(seg, { cyl.pt2, cyl.radius }, interPt, interNormal);
		if (isDiskIntersec2) {
			return isDiskIntersec2;
		}
	}
	
	return false;
}

bool InterSegmentFiniteCylinder(Segment seg, Capsule cyl, Vector3& interPt, Vector3& interNormal) {
	
}


bool IntersecRoundedBox(Segment seg, RoundedBox roundedBox) {
	bool isIntersec = false;
	Vector3 interPt = { FLT_MAX };
	Vector3 interNormal = { FLT_MAX };

	Quaternion qUp = QuaternionIdentity();
	Quaternion qLeft = QuaternionFromAxisAngle({ 1, 0, 0 }, PI * 0.5f);
	Quaternion qRight = QuaternionFromAxisAngle({ 1, 0, 0 }, PI * -0.5f);
	Quaternion qFront = QuaternionFromAxisAngle({ 0, 0, 1 }, PI * -0.5f);
	Quaternion qBack = QuaternionFromAxisAngle({ 0, 0, 1 }, PI * 0.5f);
	
	Quaternion qX = QuaternionFromAxisAngle({ 1, 0, 0 }, PI);
	Quaternion qY = QuaternionFromAxisAngle({ 0, 1, 0 }, 0);
	Quaternion qDown = QuaternionMultiply(qX, qY);

	Referential referentialQuadUp = Referential({ -0.5f, 0.75f, 0 }); // ok
	Referential referentialQuadFront = Referential({ 0.25f, 0, 0 }); // ok
	Referential referentialQuadBack = Referential({ -1.25f, 0, 0 }); // ok
	Referential referentialQuadLeft = Referential({ -0.5f, 0, 0.75f });
	Referential referentialQuadRight = Referential({ -0.5f, 0, -0.75f }); // ok
	Referential referentialQuadDown = Referential({ -0.5f, -0.75f, 0 });

	referentialQuadUp.RotateByQuaternion(qUp);
	referentialQuadFront.RotateByQuaternion(qFront);
	referentialQuadBack.RotateByQuaternion(qBack);
	referentialQuadLeft.RotateByQuaternion(qLeft);
	referentialQuadDown.RotateByQuaternion(qDown);

	Quad quadUp = { referentialQuadUp, {1, 1, 1} };
	Quad quadFront = { referentialQuadFront, {1, 1, 1} };
	Quad quadBack = { referentialQuadBack, {1, 1, 1} };
	Quad quadLeft = { referentialQuadLeft, {1, 1, 1} };
	Quad quadRight = { referentialQuadRight, {1, 1, 1} };
	Quad quadDown = { referentialQuadDown, {1, 1, 1} };

	Vector3 tmpInterPt;
	Vector3 tmpInterSeg;
	isIntersec = InterSegQuad(seg, quadUp, tmpInterPt, tmpInterSeg);

	if (isIntersec && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(interPt, seg.pt1)) {
		interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
	}

	isIntersec = InterSegQuad(seg, quadFront, tmpInterPt, tmpInterSeg);

	if (isIntersec && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(interPt, seg.pt1)) {
		interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
	}

	isIntersec = InterSegQuad(seg, quadBack, tmpInterPt, tmpInterSeg);

	if (isIntersec && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(interPt, seg.pt1)) {
		interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
	}

	isIntersec = InterSegQuad(seg, quadLeft, tmpInterPt, tmpInterSeg);

	if (isIntersec && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(interPt, seg.pt1)) {
		interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
	}

	isIntersec = InterSegQuad(seg, quadRight, tmpInterPt, tmpInterSeg);

	if (isIntersec && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(interPt, seg.pt1)) {
		interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
	}

	isIntersec = InterSegQuad(seg, quadDown, tmpInterPt, tmpInterSeg);

	if (isIntersec && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(interPt, seg.pt1)) {
		interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
	}

	Referential referentiaCapsUp = Referential({ 0, 0.5f, -0.5f });
	Referential referentialCapsDown = Referential({ 0, -0.5f, -0.5f });
	Referential referentialCapsRight = Referential({ 0, -0.5f, -0.5f });
	Referential referentialCapsLeft = Referential({ 0, -0.5f, 0.5f });

	Capsule capsuleUp = { referentiaCapsUp, 0.25f };
	Capsule capsuleDown = { referentialCapsDown, 0.25f };
	Capsule capsuleRight = { referentialCapsRight, 0.25f };
	Capsule capsuleLeft = { referentialCapsLeft, 0.25f };

	referentiaCapsUp.RotateByQuaternion(qLeft);
	referentialCapsRight.RotateByQuaternion(qUp);
	referentialCapsLeft.RotateByQuaternion(qUp);
	referentialCapsDown.RotateByQuaternion(qLeft);

	isIntersec = InterSegCas(seg, quadUp, tmpInterPt, tmpInterSeg);

	if (isIntersec && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(interPt, seg.pt1)) {
		interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
	}

	isIntersec = InterSegQuad(seg, quadFront, tmpInterPt, tmpInterSeg);

	if (isIntersec && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(interPt, seg.pt1)) {
		interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
	}

	isIntersec = InterSegQuad(seg, quadBack, tmpInterPt, tmpInterSeg);

	if (isIntersec && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(interPt, seg.pt1)) {
		interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
	}

	isIntersec = InterSegQuad(seg, quadLeft, tmpInterPt, tmpInterSeg);

	if (isIntersec && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(interPt, seg.pt1)) {
		interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
	}

	isIntersec = InterSegQuad(seg, quadRight, tmpInterPt, tmpInterSeg);

	if (isIntersec && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(interPt, seg.pt1)) {
		interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
	}

	isIntersec = InterSegQuad(seg, quadDown, tmpInterPt, tmpInterSeg);

	if (isIntersec && Vector3Distance(tmpInterPt, seg.pt1) < Vector3Distance(interPt, seg.pt1)) {
		interPt = { tmpInterPt.x, tmpInterPt.y, tmpInterPt.z };
	}

	Referential referentiaCapsUpBack = Referential({ -1, 0.5f, -0.5f });
	Referential referentialCapsDownBack = Referential({ -1, -0.5f, -0.5f });
	Referential referentialCapsLeftBack = Referential({ -1, -0.5f, -0.5f });
	Referential referentialCapsRightBack = Referential({ -1, -0.5f, 0.5f });

	Capsule capsuleUpBack = { referentiaCapsUpBack, 0.25f };
	Capsule capsuleDownBack = { referentialCapsDownBack, 0.25f };
	Capsule capsuleLeftBack = { referentialCapsLeftBack, 0.25f };
	Capsule capsuleRightBack = { referentialCapsRightBack, 0.25f };

	referentiaCapsUpBack.RotateByQuaternion(qLeft);
	referentialCapsDownBack.RotateByQuaternion(qLeft);
	referentialCapsLeftBack.RotateByQuaternion(qUp);
	referentialCapsRightBack.RotateByQuaternion(qUp);

	Referential referentiaCapsUpLeft = Referential({ -1, 0.5f, -0.5f });
	Referential referentialCapsDownLeft = Referential({ -1, -0.5f, -0.5f });

	Capsule capsuleUpRight = { referentiaCapsUpLeft, 0.25f };
	Capsule capsuleDownRight = { referentialCapsDownLeft, 0.25f };

	referentiaCapsUpLeft.RotateByQuaternion(qFront);
	referentialCapsDownLeft.RotateByQuaternion(qFront);

	Referential referentiaCapsUpRight = Referential({ -1, 0.5f, 0.5f });
	Referential referentialCapsDownRight = Referential({ -1, -0.5f, 0.5f });

	Capsule capsuleUpLeft = { referentiaCapsUpRight, 0.25f };
	Capsule capsuleDownLeft = { referentialCapsDownRight, 0.25f };

	referentiaCapsUpRight.RotateByQuaternion(qFront);
	referentialCapsDownRight.RotateByQuaternion(qFront);
}
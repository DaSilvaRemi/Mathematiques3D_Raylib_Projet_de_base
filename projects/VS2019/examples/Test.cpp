#include "Test.h"

void TestDisplaySegment() {
	TestDisplaySegment({ { -4, 0, -5 } , { 4, 0, 5 } } , BLUE);
}

void TestDisplaySegment(Segment segment, Color primaryColor) {
	Quaternion qOrient = QuaternionFromAxisAngle({ 1, 0, 0 }, PI * .5f);
	TestDisplaySegment(qOrient, segment, primaryColor);
}

void TestDisplaySegment(Quaternion q, Segment segment, Color primaryColor) {
	MyDrawSegment(q, segment, primaryColor);
}

void TestDisplaySphere() {
	Quaternion qOrient = QuaternionFromAxisAngle({ 1,0,0 }, PI * .5f);
	Sphere sphere = { {0, 0, 0}, 2 };
	TestDisplaySphere(qOrient, sphere);
}

void TestDisplaySphere(Quaternion q, Sphere sphere) {
	TestDisplaySphere(q, sphere, 40, 20, BLUE, WHITE);
}

void TestDisplaySphere(Quaternion q, Sphere sphere, int nSegmentTheta, int nSegmentPhi, Color primaryColor, Color secondaryColor) {
	MyDrawSphereEx2(q, sphere, nSegmentTheta, nSegmentPhi, primaryColor);
	MyDrawSphereWiresEx2(q, sphere, nSegmentTheta, nSegmentPhi, secondaryColor);
}

void TestDisplayQuad() {
	Quaternion qOrient = QuaternionFromAxisAngle({ 0,0, 0 }, PI * .2f);
	Quad quad = { Referential({0, 0, 0}), {4, 1, 4} };
	TestDisplayQuad(qOrient, quad);
}

void TestDisplayQuad(Quaternion q, Quad quad) {
	TestDisplayQuad(q, quad, BLUE, BLUE);
}

void TestDisplayQuad(Quaternion q, Quad quad, Color primaryColor, Color secondaryColor) {
	MyDrawQuad2(q, quad, primaryColor);
	MyDrawQuadWire2(q, quad, secondaryColor);
}

void TestDisplayPlane() {
	Quaternion qOrient = QuaternionFromAxisAngle({ 0,0, 0 }, PI * .2f);
	Plane plane = { Vector3RotateByQuaternion({ 0,1,0 }, qOrient), 2.0f };
	TestDisplayPlane(qOrient, plane);
}

void TestDisplayPlane(Quaternion q, Plane plane) {
	TestDisplayPlane(q, plane, { 2, 2 }, BLUE, WHITE);
}

void TestDisplayPlane(Quaternion q, Plane plane, Vector2 size, Color primaryColor, Color secondaryColor) {
	MyDrawQuad2(q, { Referential(Vector3Scale(plane.normal, plane.d)), { size.x, size.y, 1 } }, primaryColor);
	MyDrawQuadWire2(q, { Referential(Vector3Scale(plane.normal, plane.d)), { size.x, size.y, 1 } }, secondaryColor);
}

void TestDisplayCylinder() {
	Quaternion qOrient = QuaternionFromAxisAngle({ 0,0, 0 }, PI * .2f);
	Referential ref = Referential({ 0, 0, 0 });
	ref.RotateByQuaternion(qOrient);
	Cylinder cylinder = Cylinder(ref, 2, 4);
	TestDisplayCylinder(qOrient, cylinder);
}

void TestDisplayCylinder(Quaternion q, Cylinder cylinder) {
	TestDisplayCylinder(q, cylinder, 25, true, BLUE, RED);
}

void TestDisplayCylinder(Quaternion q, Cylinder cylinder, int nSegmensTheta, bool drawCaps, Color primaryColor, Color secondaryColor) {
	MyDrawCylinder(q, cylinder, nSegmensTheta, drawCaps, primaryColor);
	MyDrawCylinderWires(q, cylinder, nSegmensTheta, drawCaps, secondaryColor);
}

void TestDisplayDisk() {
	Quaternion qOrient = QuaternionFromAxisAngle({ 0, 0, 0 }, PI * .2f);
	Referential diskReferential = Referential({ 1, 1, 1 });
	Disk disk = { diskReferential, 2.0f };
	TestDisplayDisk(qOrient, disk);
}

void TestDisplayDisk(Quaternion q, Disk disk) {
	TestDisplayDisk(q, disk, 25, BLUE, WHITE);
}

void TestDisplayDisk(Quaternion q, Disk disk, int nSegmensTheta, Color primaryColor, Color secondaryColor) {
	MyDrawDisk(q, disk.referential.origin, disk.radius, nSegmensTheta, primaryColor);
	MyDrawDiskWires(q, disk.referential.origin, disk.radius, nSegmensTheta, secondaryColor);
}

void TestDisplayCaps() {
	Quaternion q = QuaternionIdentity();
	Referential referential = Referential({ 0,1,0 });
	Capsule capsule = { referential, 3 };
	TestDisplayCaps(q, capsule, BLUE, WHITE);
}

void TestDisplayCaps(Quaternion q, Capsule capsule, Color primaryColor, Color secondaryColor) {
	MyDrawCapsule(q, capsule, primaryColor);
	MyDrawCapsuleWires(q, capsule, secondaryColor);
}

void TestDisplayRoundBox() {
	Quaternion q = QuaternionFromAxisAngle({ 0, 0, 0 }, PI * .2f);
	RoundedBox roundedBox = { Referential({ 0, 0, 0 }), {5, 5, 5}, 1 };
	TestDisplayRoundBox(q, roundedBox);
}

void TestDisplayRoundBox(Quaternion q, RoundedBox roundedBox) {
	TestDisplayRoundBox(q, roundedBox, BLUE, WHITE);
}

void TestDisplayRoundBox(Quaternion q, RoundedBox roundedBox, Color primaryColor, Color secondaryColor) {
	MyDrawRoundBox(q, roundedBox, primaryColor);
	//MyDrawRoundBoxWires(q, roundedBox, secondaryColor);
}


void TestDisplayBox() {
	Quaternion q = QuaternionFromAxisAngle({ 0, 0, 0 }, PI * .2f);
	Box box = { Referential({ 0, 5, 0 }), { 2, 2, 2 } };
	TestDisplayBox(q, box);
}

void TestDisplayBox(Quaternion q, Box box) {
	TestDisplayBox(q, box, BLUE, WHITE);
}

void TestDisplayBox(Quaternion q, Box box, Color primaryColor, Color secondaryColor) {
	MyDrawBox(q, box, primaryColor);
	MyDrawBoxWires(q, box, secondaryColor);
}

//TEST INTERSECTION SEGMENT PLANE
void TestIntersecSegmentPlane(float time) {
	Quaternion qOrient = QuaternionFromAxisAngle({ 1, 0, 0 }, PI * .2f * time);
	Segment segment = { { -4, 0, -5 } , { 4, 0, 5 } };
	Plane plane = { Vector3RotateByQuaternion({ 0,1,0 }, qOrient), 2.0f };

	Vector3 interSecPt = { 0, 0, 0 };
	Vector3 interSecNormal = { 0, 0, 0 };
	bool planeHaveIntersec = InterSegPlane(segment, plane, interSecPt, interSecNormal);

	TestDisplaySegment(qOrient, segment, PURPLE);
	TestDisplayPlane(qOrient, plane);

	if (planeHaveIntersec) {
		TestDisplaySegment(qOrient, { interSecPt, interSecNormal }, DARKPURPLE);
		DrawSphere(interSecPt, .2f, DARKBROWN);
	}	
}

//TEST INTERSECTION SEGMENT SPHERE
void TestIntersecSegmentSphere(float time) {
	Quaternion qOrient = QuaternionFromAxisAngle({ 1, 0, 0 }, PI * .5f * time);
	Segment segment = { { -4, 0, -5 } , { 4, 0, 5 } };
	Sphere sphere = { {0, 0, 0}, 2 };
	

	Vector3 interSecPt = { 0, 0, 0 };
	Vector3 interSecNormal = { 0, 0, 0 };
	bool sphereHaveIntersec = InterSegSphere(segment, sphere, interSecPt, interSecNormal);

	TestDisplaySegment(qOrient, segment, PURPLE);
	TestDisplaySphere(qOrient, sphere);

	if (sphereHaveIntersec) {
		TestDisplaySegment(qOrient, { interSecPt, interSecNormal }, DARKPURPLE);
		DrawSphere(interSecPt, .2f, DARKBROWN);
	}
}

//TEST INTERSECTION SEGMENT QUAD
void TestIntersecInterSegmentQuad(float time) {

	Quaternion qOrient = QuaternionFromAxisAngle({ 1, 0,1 }, PI * .5f * time);
	Segment segment = { { -4, 0, -5 } , { 4, 0, 5 } };
	Quad quad = { Referential({0, 0, 0}), {4, 1, 4} };

	Vector3 interSecPt = { 0, 0, 0 };
	Vector3 interSecNormal = { 0, 0, 0 };
	bool quadHaveIntersec = InterSegQuad(segment, quad, interSecPt, interSecNormal);

	TestDisplaySegment(segment, PURPLE);
	TestDisplayQuad(qOrient, quad);

	if (quadHaveIntersec) {
		TestDisplaySegment(qOrient, { interSecPt, interSecNormal }, DARKPURPLE);
		DrawSphere(interSecPt, .2f, DARKBROWN);
	}
}

//TEST INTERSECTION SEGMENT DISK
void TestIntersecInterSegmentDisk(float time) {

	Quaternion qOrient = QuaternionFromAxisAngle({ 1, 0,0 }, PI * .5f * time);
	Segment segment = { { -4, 0, -5 } , { 4, 0, 5 } };
	Disk disk = { Referential({0, 0, 0}), 5 };

	Vector3 interSecPt = { 0, 0, 0 };
	Vector3 interSecNormal = { 0, 0, 0 };
	bool diskHaveIntersec = InterSegDisk(segment, disk, interSecPt, interSecNormal);

	TestDisplaySegment(segment, PURPLE);
	TestDisplayDisk(qOrient, disk);

	if (diskHaveIntersec) {
		TestDisplaySegment(qOrient, { interSecPt, interSecNormal }, DARKPURPLE);
		DrawSphere(interSecPt, .2f, DARKBROWN);
	}
	}

void TestIntersecSegmentInfiniteCylinder(float time) {

	Quaternion qOrient = QuaternionFromAxisAngle({ 1, 0, 0 }, PI * 0.5f);
	//Quaternion qOrient = QuaternionIdentity();
	Segment segment = { { -4, 3, -5 } , { 4, 0, 5 } };

	Referential ref = Referential({ 0, 0, 0 });
	ref.RotateByQuaternion(qOrient);
	Cylinder cylinder = Cylinder(ref, 2, 4);

	Vector3 interSecPt = { 0, 0, 0 };
	Vector3 interSecNormal = { 0, 0, 0 };
	bool cylinderHaveIntersec = InterSegmentInfiniteCylinder(segment, cylinder, interSecPt, interSecNormal);

	DrawLine3D(segment.pt1, segment.pt2, RED);
	TestDisplayCylinder(qOrient, cylinder);

	if (cylinderHaveIntersec) {
		//TestDisplaySegment(qOrient, { Vector3Subtract(interSecNormal, interSecPt), interSecPt }, DARKPURPLE);
		DrawLine3D(interSecPt, Vector3Add(interSecPt, interSecNormal), RED);
		DrawSphere(interSecPt, .2f, DARKBROWN);
	}
}

void TestIntersecParalleleSegmentFiniteCylinder(float time) {

	Quaternion qOrient = QuaternionFromAxisAngle({1, 0, 0 }, PI * .5f);
	Segment segment = { { -8, 3, 0 } , { 8, 3, 0} };

	Referential ref = Referential({ 0, 0, 0 });
	ref.RotateByQuaternion(qOrient);
	Cylinder cylinder = Cylinder(ref, 2, 4);

	Vector3 interSecPt = { 0, 0, 0 };
	Vector3 interSecNormal = { 0, 0, 0 };
	bool cylinderHaveIntersec = InterSegmentFiniteCylinder(segment, cylinder, interSecPt, interSecNormal);

	TestDisplaySegment(qOrient, segment, PURPLE);
	//DrawLine3D(segment.pt1, segment.pt2, PURPLE);
	TestDisplayCylinder(qOrient, cylinder);

	if (cylinderHaveIntersec) {
		DrawLine3D(interSecPt, Vector3Add(interSecPt, interSecNormal), RED);
		DrawSphere(interSecPt, .2f, DARKBROWN);
	}
}

void TestIntersecSegmentFiniteCylinderDisk(float time) {

	Quaternion qOrient = QuaternionFromAxisAngle({ 1, 0,0 }, PI * .5f * time);
	Segment segment = { { -4, 0, -5 } , { -4, 0, 5 } };

	Referential ref = Referential({ 0, 0, 0 });
	ref.RotateByQuaternion(qOrient);
	Cylinder cylinder = Cylinder(ref, 2, 4);

	Vector3 interSecPt = { 0, 0, 0 };
	Vector3 interSecNormal = { 0, 0, 0 };
	bool sphereHaveIntersec = InterSegmentFiniteCylinder(segment, cylinder, interSecPt, interSecNormal);

	TestDisplaySegment(segment, PURPLE);
	TestDisplayCylinder(qOrient, cylinder);

	if (sphereHaveIntersec) {
		DrawLine3D(interSecPt, Vector3Add(interSecPt, interSecNormal), RED);
		DrawSphere(interSecPt, .2f, DARKBROWN);
	}
}

void TestIntersecSegmentDiskNoInfiniteCylinder(float time) {

	Quaternion qOrient = QuaternionFromAxisAngle({ 1, 0,0 }, PI * .5f * time);
	Segment segment = { { -4, 0, -5 } , { -4, 0, 5 } };

	Referential ref = Referential({ 0, 0, 0 });
	ref.RotateByQuaternion(qOrient);
	Cylinder cylinder = Cylinder(ref, 2, 4);

	Vector3 interSecPt = { 0, 0, 0 };
	Vector3 interSecNormal = { 0, 0, 0 };
	bool sphereHaveIntersec = InterSegmentInfiniteCylinder(segment, cylinder, interSecPt, interSecNormal);

	TestDisplaySegment(segment, PURPLE);
	TestDisplayCylinder(qOrient, cylinder);

	if (sphereHaveIntersec) {
		DrawLine3D(interSecPt, Vector3Add(interSecPt, interSecNormal), RED);
		DrawSphere(interSecPt, .2f, DARKBROWN);
	}
}

void TestIntersecSegmentFiniteCylinderNoDisk(float time) {

	Quaternion qOrient = QuaternionFromAxisAngle({ 1, 0,0 }, PI * .5f * time);
	Segment segment = { {2, 0, -5 } , { 2, 0, 5 } };

	Referential ref = Referential({ 0, 0, 0 });
	ref.RotateByQuaternion(qOrient);
	Cylinder cylinder = Cylinder(ref, 2, 4);

	Vector3 interSecPt = { 0, 0, 0 };
	Vector3 interSecNormal = { 0, 0, 0 };
	bool sphereHaveIntersec = InterSegmentFiniteCylinder(segment, cylinder, interSecPt, interSecNormal);

	TestDisplaySegment(segment, PURPLE);
	TestDisplayCylinder(qOrient, cylinder, 20, false, BLUE, WHITE);

	if (sphereHaveIntersec) {
		DrawLine3D(interSecPt, Vector3Add(interSecPt, interSecNormal), RED);
		DrawSphere(interSecPt, .2f, DARKBROWN);
	}
}

void TestIntersecSegmentFiniteCylinder(float time) {
	TestIntersecSegmentFiniteCylinderNoDisk(time);
	TestIntersecSegmentDiskNoInfiniteCylinder(time);
	TestIntersecSegmentFiniteCylinderDisk(time);
	TestIntersecParalleleSegmentFiniteCylinder(time);
}

void TestAllDisplay() {
	TestDisplaySegment();
	TestDisplaySphere();
	TestDisplayPlane();
	TestDisplayCylinder();
	TestDisplayDisk();
	TestDisplayQuad();
	TestDisplayCaps();
	TestDisplayRoundBox();
}

void TestAllIntersec() {
	TestIntersecSegmentPlane();
	TestIntersecSegmentSphere();
	TestIntersecSegmentInfiniteCylinder();
	TestIntersecSegmentFiniteCylinder();
	TestIntersecInterSegmentDisk();
	TestIntersecInterSegmentQuad();
	TestIntersecSegmentFiniteCylinder();
}

void TestAll() {
	TestAllDisplay();
	TestAllIntersec();
}
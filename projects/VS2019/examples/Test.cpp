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
	MyDrawQuad2(q, quad.referential.origin, { quad.extension.x, quad.extension.z }, primaryColor);
	MyDrawQuadWire2(q, quad.referential.origin, { quad.extension.x, quad.extension.z }, secondaryColor);
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
	MyDrawQuad2(q, Vector3Scale(plane.normal, plane.d), size, primaryColor);
	MyDrawQuadWire2(q, Vector3Scale(plane.normal, plane.d), size, secondaryColor);
}

void TestDisplayCylinder() {
	Quaternion qOrient = QuaternionFromAxisAngle({ 0,0, 0 }, PI * .2f);
	Cylinder cylinder = { {0, 0, 0}, {0, 4, 0}, 2 };
	TestDisplayCylinder(qOrient, cylinder);
}

void TestDisplayCylinder(Quaternion q, Cylinder cylinder) {
	TestDisplayCylinder(q, cylinder, 25, true, BLUE, WHITE);
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

void TestDisplayBox() {
	Quaternion q = QuaternionFromAxisAngle({ 0, 0, 0 }, PI * .2f);
	Vector3 center = { 0,1,0 };
	Vector3 size = { 2, 2, 2 };
	TestDisplayBox(q, center, size);
}

void TestDisplayBox(Quaternion q, Vector3 center, Vector3 size) {
	TestDisplayBox(q, center, size, BLUE, WHITE);
}

void TestDisplayBox(Quaternion q, Vector3 center, Vector3 size, Color primaryColor, Color secondaryColor) {
	MyDrawBox(q, center, size, primaryColor);
	MyDrawBoxWires(q, center, size, secondaryColor);
}

//TEST INTERSECTION SEGMENT PLANE
void TestIntersecSegmentPlane(float time) {
	Quaternion qOrient = QuaternionFromAxisAngle({ 1, 0, 0 }, PI * .2f * time);
	Segment segment = { { -4, 0, -5 } , { 4, 0, 5 } };
	Plane plane = { Vector3RotateByQuaternion({ 0,1,0 }, qOrient), 2.0f };

	Vector3 interSectPt = { 0, 0, 0 };
	Vector3 interSecNormal = { 0, 0, 0 };
	bool planeHaveIntersec = InterSegPlane(segment, plane, interSectPt, interSecNormal);

	TestDisplaySegment(qOrient, segment, PURPLE);
	TestDisplayPlane(qOrient, plane);

	if (planeHaveIntersec) {
		TestDisplaySegment(qOrient, { interSectPt, interSecNormal }, DARKPURPLE);
		DrawSphere(interSectPt, .2f, DARKBROWN);
	}	
}

//TEST INTERSECTION SEGMENT SPHERE
void TestIntersecSegmentSphere(float time) {
	Quaternion qOrient = QuaternionFromAxisAngle({ 1, 0, 0 }, PI * .5f * time);
	Segment segment = { { -4, 0, -5 } , { 4, 0, 5 } };
	Sphere sphere = { {0, 0, 0}, 2 };
	

	Vector3 interSectPt = { 0, 0, 0 };
	Vector3 interSecNormal = { 0, 0, 0 };
	bool sphereHaveIntersec = InterSegSphere(segment, sphere, interSectPt, interSecNormal);

	TestDisplaySegment(qOrient, segment, PURPLE);
	TestDisplaySphere(qOrient, sphere);

	if (sphereHaveIntersec) {
		TestDisplaySegment(qOrient, { interSectPt, interSecNormal }, DARKPURPLE);
		DrawSphere(interSectPt, .2f, DARKBROWN);
	}
}

//TEST INTERSECTION SEGMENT QUAD
void TestIntersecInterSegmentQuad(float time) {

	Quaternion qOrient = QuaternionFromAxisAngle({ 1, 0,1 }, PI * .5f * time);
	Segment segment = { { -4, 0, -5 } , { 4, 0, 5 } };
	Quad quad = { Referential({0, 0, 0}), {4, 1, 4} };

	Vector3 interSectPt = { 0, 0, 0 };
	Vector3 interSecNormal = { 0, 0, 0 };
	bool quadHaveIntersec = InterSegQuad(segment, quad, interSectPt, interSecNormal);

	TestDisplaySegment(segment, PURPLE);
	TestDisplayQuad(qOrient, quad);

	if (quadHaveIntersec) {
		TestDisplaySegment(qOrient, { interSectPt, interSecNormal }, DARKPURPLE);
		DrawSphere(interSectPt, .2f, DARKBROWN);
	}
}

//TEST INTERSECTION SEGMENT DISK
void TestIntersecInterSegmentDisk(float time) {

	Quaternion qOrient = QuaternionFromAxisAngle({ 1, 0,0 }, PI * .5f * time);
	Segment segment = { { -4, 0, -5 } , { 4, 0, 5 } };
	Disk disk = { Referential({0, 0, 0}), 5 };

	Vector3 interSectPt = { 0, 0, 0 };
	Vector3 interSecNormal = { 0, 0, 0 };
	bool diskHaveIntersec = InterSegDisk(segment, disk, interSectPt, interSecNormal);

	TestDisplaySegment(segment, PURPLE);
	TestDisplayDisk(qOrient, disk);

	if (diskHaveIntersec) {
		TestDisplaySegment(qOrient, { interSectPt, interSecNormal }, DARKPURPLE);
		DrawSphere(interSectPt, .2f, DARKBROWN);
	}
	}

void TestIntersecSegmentCylinder(float time) {

	Quaternion qOrient = QuaternionFromAxisAngle({ 1, 0,0 }, PI * .5f * time);
	Segment segment = { { -4, 0, -5 } , { 4, 0, 5 } };
	Cylinder cylinder = { {0, 0, 0}, {0, 4, 0}, 2 };

	Vector3 interSectPt = { 0, 0, 0 };
	Vector3 interSecNormal = { 0, 0, 0 };
	bool sphereHaveIntersec = InterSegmentInfiniteCylinder(segment, cylinder, interSectPt, interSecNormal);

	TestDisplaySegment(segment, PURPLE);
	TestDisplayCylinder(qOrient, cylinder);

	if (sphereHaveIntersec) {
		TestDisplaySegment(qOrient, { interSectPt, interSecNormal }, DARKPURPLE);
		DrawSphere(interSectPt, .2f, DARKBROWN);
	}
}

void TestAllDisplay() {
	TestDisplaySegment();
	TestDisplaySphere();
	TestDisplayPlane();
	TestDisplayCylinder();
	TestDisplayDisk();
	TestDisplayQuad();
	TestDisplayBox();
}

void TestAllIntersec() {
	TestIntersecSegmentPlane();
	TestIntersecSegmentSphere();
	TestIntersecSegmentCylinder();
	TestIntersecInterSegmentDisk();
	TestIntersecInterSegmentQuad();
}

void TestAll() {
	TestAllDisplay();
	TestAllIntersec();
}
#pragma once

#include "raylib.h"
#include <raymath.h>
#include "rlgl.h"
#include <math.h>
#include <float.h>
#include <vector>
#include <iostream>
#include "Struct.h"
#include "Tools.h"
#include "MyDrawMethods.h"

//TestDrawing

void TestDisplaySegment();
void TestDisplaySegment(Segment segment, Color primaryColor);
void TestDisplaySegment(Quaternion q, Segment segment, Color primaryColor);

void TestDisplaySphere();
void TestDisplaySphere(Quaternion q, Sphere sphere);
void TestDisplaySphere(Quaternion q, Sphere sphere, int nSegmentTheta, int nSegmentPhi, Color primaryColor, Color secondaryColor);

void TestDisplayQuad();
void TestDisplayQuad(Quaternion q, Quad quad);
void TestDisplayQuad(Quaternion q, Quad quad, Color primaryColor, Color secondaryColor);

void TestDisplayPlane();
void TestDisplayPlane(Quaternion q, Plane plane);
void TestDisplayPlane(Quaternion q, Plane plane, Vector2 size, Color primaryColor, Color secondaryColor);

void TestDisplayCylinder();
void TestDisplayCylinder(Quaternion q, Cylinder cylinder);
void TestDisplayCylinder(Quaternion q, Cylinder cylinder, int nSegmensTheta, bool drawCaps, Color primaryColor, Color secondaryColor);

void TestDisplayDisk();
void TestDisplayDisk(Quaternion q, Disk disk);
void TestDisplayDisk(Quaternion q, Disk disk, int nSegmensTheta, Color primaryColor, Color secondaryColor);

void TestDisplayBox();
void TestDisplayBox(Quaternion q, Vector3 center, Vector3 size);
void TestDisplayBox(Quaternion q, Vector3 center, Vector3 size, Color primaryColor, Color secondaryColor);

//Test Intersection

void TestIntersecSegmentPlane(float time = 1);
void TestIntersecSegmentSphere(float time = 1);
void TestIntersecInterSegmentQuad(float time = 1);
void TestIntersecInterSegmentDisk(float time = 1);
void TestIntersecSegmentInfiniteCylinder(float time = 1);

void TestIntersecParalleleSegmentFiniteCylinder(float time = 1);
void TestIntersecSegmentInfiniteCylinderNoDisk(float time = 1);
void TestIntersecSegmentFiniteCylinderDisk(float time = 1);
void TestIntersecSegmentDiskNoInfiniteCylinder(float time = 1);
void TestIntersecSegmentFiniteCylinderNoDisk(float time = 1);
void TestIntersecSegmentFiniteCylinder(float time = 1);

void TestAllDisplay();
void TestAllIntersec();
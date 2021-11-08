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

//Test Intersection

void TestIntersecSegmentPlane(float time = 1);
void TestIntersecSegmentSphere(float time = 1);
void TestIntersecSegmentCylinder(float time = 1);

void TestAllDisplay();
void TestAllIntersec();
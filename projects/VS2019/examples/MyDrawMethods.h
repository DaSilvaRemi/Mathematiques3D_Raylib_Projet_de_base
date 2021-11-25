#include "raylib.h"
#include <raymath.h>
#include "rlgl.h"
#include <math.h>
#include <float.h>
#include <vector>
#include <iostream>
#include "Tools.h"
#include "Struct.h"

void MyDrawSegment(Quaternion q, Segment seg, Color color);

void MyDrawCylinderPortion(Quaternion q, Cylinder cyl, float startTheta, float endTheta, int nSegmentsTheta, bool drawCaps, Color color);

void MyDrawDiskPortion(Quaternion q, Vector3 center, float radius, float startTheta, float endTheta, int nSegmentsTheta, Color color);

void MyDrawDiskWiresPortion(Quaternion q, Vector3 center, float radius, float startTheta, float endTheta, int nSegmentsTheta, Color color);

void MyDrawSpherePortion(Quaternion q, Sphere sph, float startTheta, float endTheta, float startPhi, float endPhi, int nSegmentsTheta, int nSegmentsPhi, Color color);

void MyDrawSphereWiresPortion(Quaternion q, Sphere sph, float startTheta, float endTheta, float startPhi, float endPhi, int nSegmentsTheta, int nSegmentsPhi, Color color);

void MyDrawSphereEx2(Quaternion q, Sphere sph, int nSegmentsTheta, int nSegmentsPhi, Color color);

void MyDrawSphereWiresEx2(Quaternion q, Sphere sph, int nSegmentsTheta, int nSegmentsPhi, Color color);

void MyDrawQuad2(Quaternion q, Vector3 center, Vector2 size, Color color);

void MyDrawQuadWire2(Quaternion q, Vector3 center, Vector2 size, Color color);

void MyDrawCylinder(Quaternion q, Cylinder cyl, int nSegmentsTheta, bool drawCaps, Color color);

void MyDrawCylinderWires(Quaternion q, Cylinder cyl, int nSegmentsTheta, bool drawCaps, Color color);

void MyDrawCylinderPortion(Quaternion q, Cylinder cyl, float startTheta, float endTheta, int nSegmentsTheta, bool drawCaps, Color color);

void MyDrawCylinderWiresPortion(Quaternion q, Cylinder cyl, float startTheta, float endTheta, int nSegmentsTheta, bool drawCaps, Color color);

void MyDrawDisk(Quaternion q, Vector3 center, float radius, int nSegmentsTheta, Color color);

void MyDrawDiskWires(Quaternion q, Vector3 center, float radius, int nSegmentsTheta, Color color);

void MyDrawDiskPortion(Quaternion q, Vector3 center, float radius, float startTheta, float endTheta, int nSegmentsTheta, Color color);

void MyDrawDiskWiresPortion(Quaternion q, Vector3 center, float radius, float startTheta, float endTheta, int nSegmentsTheta, Color color);

void MyDrawCapsule(Quaternion q, Capsule capsule, Color color);

void MyDrawCapsuleWires(Quaternion q, Capsule capsule, Color color);

void MyDrawRoundBox(Quaternion q, Vector3 center, Vector3 size, Color color);

void MyDrawRoundBoxWires(Quaternion q, Vector3 center, Vector3 size, Color color);

void MyDrawBox(Quaternion q, Vector3 center, Vector3 size, Color color);

void MyDrawBoxWires(Quaternion q, Vector3 center, Vector3 size, Color color);
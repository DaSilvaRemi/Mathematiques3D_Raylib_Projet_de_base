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

void MyDrawCylinderPortion(Cylinder cyl, float startTheta, float endTheta, int nSegmentsTheta, bool drawCaps, Color color);

void MyDrawDiskPortion(Disk disk, float startTheta, float endTheta, int nSegmentsTheta, Color color);

void MyDrawDiskWiresPortion(Disk disk, float startTheta, float endTheta, int nSegmentsTheta, Color color);

void MyDrawSpherePortion(Quaternion q, Sphere sph, float startTheta, float endTheta, float startPhi, float endPhi, int nSegmentsTheta, int nSegmentsPhi, Color color);

void MyDrawSphereWiresPortion(Quaternion q, Sphere sph, float startTheta, float endTheta, float startPhi, float endPhi, int nSegmentsTheta, int nSegmentsPhi, Color color);

void MyDrawSphereEx2(Quaternion q, Sphere sph, int nSegmentsTheta, int nSegmentsPhi, Color color);

void MyDrawSphereWiresEx2(Quaternion q, Sphere sph, int nSegmentsTheta, int nSegmentsPhi, Color color);

void MyDrawQuad2(Quad quad, Color color);

void MyDrawQuadWire2(Quad quad, Color color);

void MyDrawCylinder(Cylinder cyl, int nSegmentsTheta, bool drawCaps, Color color);

void MyDrawCylinderWires(Cylinder cyl, int nSegmentsTheta, bool drawCaps, Color color);

void MyDrawCylinderPortion(Cylinder cyl, float startTheta, float endTheta, int nSegmentsTheta, bool drawCaps, Color color);

void MyDrawCylinderWiresPortion(Cylinder cyl, float startTheta, float endTheta, int nSegmentsTheta, bool drawCaps, Color color);

void MyDrawDisk(Disk disk, int nSegmentsTheta, Color color);

void MyDrawDiskWires(Disk disk, int nSegmentsTheta, Color color);

void MyDrawDiskPortion(Disk disk, float startTheta, float endTheta, int nSegmentsTheta, Color color);

void MyDrawDiskWiresPortion(Disk disk, float startTheta, float endTheta, int nSegmentsTheta, Color color);

void MyDrawCapsule(Capsule capsule, Color color);

void MyDrawCapsuleWires(Capsule capsule, Color color);

void MyDrawRoundedBoxV2(RoundedBox roundedBox, Color color);

void MyDrawRoundBoxWiresV2(RoundedBox roundedBox, Color color);

void MyDrawBox(Box box, Color color);

void MyDrawBoxWires(Box box, Color color);
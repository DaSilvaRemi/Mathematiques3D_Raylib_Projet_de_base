#include "raylib.h"
#include <raymath.h>
#include "rlgl.h"
#include <math.h>
#include <float.h>
#include <vector>
#include <iostream>
#include "Struct.h"

#define EPSILON 1.e-6f

Cylindrical CartesianToCylindrical(Vector3 cart);

Vector3 CylindricalToCartesian(Cylindrical cyl);

Vector3 SphericalToCartesian(Spherical sph);

Vector3 GlobalToLocalPos(Vector3 posGlobal, Referential localRef);

Vector3 GlobalToLocalVect(Vector3 vectGlobal, Referential localRef);

Vector3 LocalToGlobalPos(Vector3 localPos, Referential localRef);

Vector3 LocalToGlobalVect(Vector3 localVect, Referential localRef);

bool InterSegPlane(Segment seg, Plane plane, Vector3& interPt, Vector3& interNormal);

bool InterSegSphere(Segment seg, Sphere sphere, Vector3& interPt, Vector3& interNormal);

bool InterSegmentInfiniteCylinder(Segment seg, Cylinder cyl, Vector3& interPt, Vector3& interNormal);
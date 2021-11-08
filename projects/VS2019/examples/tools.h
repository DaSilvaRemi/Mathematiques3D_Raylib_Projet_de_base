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

bool InterSegPlane(Segment seg, Plane plane, Vector3& interPt, Vector3& interNormal);

bool InterSegSphere(Segment seg, Sphere sphere, Vector3& interPt, Vector3& interNormal);

bool InterSegmentInfiniteCylinder(Segment seg, Cylinder cyl, Vector3& interPt, Vector3& interNormal);
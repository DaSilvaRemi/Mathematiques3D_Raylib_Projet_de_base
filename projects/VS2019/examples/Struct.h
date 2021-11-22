#pragma once

#include "raylib.h"
#include <raymath.h>
#include "rlgl.h"
#include <math.h>
#include <float.h>
#include <vector>
#include <iostream>

struct Cylindrical {
	float rho;
	float theta;
	float y;

	inline Cylindrical operator+(Cylindrical a) {
		return { a.rho + rho,a.theta + theta,a.y + y };
	}
};

struct Spherical {
	float rho;
	float theta;
	float phi;

	inline Spherical operator+(Spherical a) {
		return { a.rho + rho,a.theta + theta,a.phi + phi };
	}

};

struct Referential {
	Vector3 origin;
	Vector3 i;
	Vector3 j;
	Vector3 k;

	Referential(Vector3 pt) {
		origin = pt;
		i = { 1, 0, 0 };
		j = { 0, 1, 0 };
		k = { 0, 0, 1 };
	};

	void RotateByQuaternion(Quaternion q) {
		i = Vector3RotateByQuaternion(i, q);
		j = Vector3RotateByQuaternion(j, q);
		k = Vector3RotateByQuaternion(k, q);
	}
};

struct Line {
	Vector3 pt;
	Vector3 dir;
};

struct Segment {
	Vector3 pt1;
	Vector3 pt2;
};

struct Plane {
	Vector3 normal;
	float d;

	Plane(Vector3 _normal, Vector3 O, Vector3 H) {
		Plane(_normal, Vector3Subtract(H, O));
	}

	Plane(Vector3 _normal, Vector3 OH) {
		Plane(_normal, Vector3DotProduct(OH, normal));
	}

	Plane(Vector3 _normal, float _d) {
		normal = _normal;
		d = _d;
	}
};

struct Quad {
	Referential referential;
	Vector3 extension;
};

struct Sphere {
	Vector3 omega;
	float rayon;
};

struct Cylinder {
	Vector3 pt1;
	Vector3 pt2;
	float radius;
};

struct Disk {
	Referential referential;
	float radius;
};

struct Capsule {
	Referential referential;
	float radius;
};

struct RoundexBox{
	Quad quad;
};

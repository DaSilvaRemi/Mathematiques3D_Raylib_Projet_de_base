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

	Referential() {
		origin = { 0, 0, 0 };
		i = { 1, 0, 0 };
		j = { 0, 1, 0 };
		k = { 0, 0, 1 };
	}

	Referential(Vector3 pt) {
		origin = pt;
		i = { 1, 0, 0 };
		j = { 0, 1, 0 };
		k = { 0, 0, 1 };
	}

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
	Referential ref;
	Vector3 pt1;
	Vector3 pt2;
	float radius;
	float height;

	Cylinder(Referential ref, float radius, float height) {
		this->ref = ref;
		this->radius = radius;
		this->height = height;
		this->pt1 = ref.origin;
		this->pt2 = Vector3Add(this->pt1, Vector3Scale(ref.j, height));
	}
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
	std::vector<Quad> quads;

	RoundexBox() {
		Quaternion qUp = QuaternionIdentity();
		Quaternion qLeft = QuaternionFromAxisAngle({ 1, 0, 0 }, PI * 0.5f);
		Quaternion qRight = QuaternionFromAxisAngle({ 1, 0, 0 }, PI * -0.5f);
		Quaternion qFront = QuaternionFromAxisAngle({ 0, 0, 1 }, PI * -0.5f);
		Quaternion qBack = QuaternionFromAxisAngle({ 0, 0, 1 }, PI * 0.5f);
		//Quaternion qDown = QuaternionFromAxisAngle({ 1, 0, 0 }, PI * 1);

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

		quads.push_back({ referentialQuadUp, {1, 1, 1} });
		quads.push_back({ referentialQuadFront, {1, 1, 1} });
		quads.push_back({ referentialQuadBack, {1, 1, 1} });
		quads.push_back({ referentialQuadLeft, {1, 1, 1} });
		quads.push_back({ referentialQuadRight, {1, 1, 1} });
		quads.push_back({ referentialQuadDown, {1, 1, 1} });
	}
};

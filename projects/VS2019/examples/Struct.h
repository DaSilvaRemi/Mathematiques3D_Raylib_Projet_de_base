#pragma once
#ifndef _MY_STRUCT_H
#define _MY_STRUCT_H

#include "raylib.h"
#include <raymath.h>

/**
 * @brief Cylindrical coordinates 
 */
struct Cylindrical
{
    float rho;
    float theta;
    float y;

    /**
     * @brief overload of the + operator to add two coordinates
     *
     * @param a The coordinates to add
     */
    Cylindrical operator+(Cylindrical a)
    {
        return {a.rho + rho, a.theta + theta, a.y + y};
    }
};

/**
 * @brief Spherical  coordinates 
 */
struct Spherical
{
    float rho;
    float theta;
    float phi;

    /**
     * @brief overload of the + operator to add two coordinates
     *
     * @param a The coordinates to add
     */
    Spherical operator+(Spherical a)
    {
        return {a.rho + rho, a.theta + theta, a.phi + phi};
    }
};

/**
 * @brief A referential with origin, a quaternion and 3 unitary vector
 */
struct Referential
{
    Vector3 origin;
    Vector3 i;
    Vector3 j;
    Vector3 k;
    Quaternion q;

    /**
     * @brief Default constructor of the struct we set origin to 0 and a quaternion identity
     * 
     */
    Referential()
    {
        origin = {0, 0, 0};
        i = {1, 0, 0};
        j = {0, 1, 0};
        k = {0, 0, 1};
        q = QuaternionIdentity();
    }

    /**
    * @brief Constructor of Referential
    *
    * @param pt A point in 3D space set to the origin of the referential
    */
    Referential(Vector3 pt)
    {
        origin = pt;
        i = {1, 0, 0};
        j = {0, 1, 0};
        k = {0, 0, 1};
        q = QuaternionIdentity();
    }

    /**
    * @brief Constructor of Referential
    *
    * @param pt A point in 3D space set to the origin of the referential
    * @param q A quaternion to rotate the referential
    */
    Referential(Vector3 pt, Quaternion q)
    {
        origin = pt;
        i = {1, 0, 0};
        j = {0, 1, 0};
        k = {0, 0, 1};
        this->q = QuaternionIdentity();
        this->RotateByQuaternion(q);
    }

    /**
    * @brief Rotate the referential with a Quaternion
    * @remarks It multiply the save Quaternion by the param Quaternion
    *
    * @param q A quaternion to rotate the referential
    */
    void RotateByQuaternion(Quaternion q)
    {
        i = Vector3RotateByQuaternion(i, q);
        j = Vector3RotateByQuaternion(j, q);
        k = Vector3RotateByQuaternion(k, q);
        this->q = QuaternionMultiply(this->q, q);
    }
};

/**
 * @brief A line struct with a point and a direction
 */
struct Line
{
    Vector3 pt;
    Vector3 dir;
};

/**
 * @brief A segment represented with two points 
 */
struct Segment
{
    Vector3 pt1;
    Vector3 pt2;
};

/**
 * @brief A Plane represented with a normal and a distance
 */
struct Plane
{
    Vector3 normal;
    float d;

    /**
     * @brief Default constructor of Plane
     *
     * @param _normal The normal of the plane
     * @param pt A point in the space to calculate the distance
     */
    Plane(Vector3 _normal, Vector3 pt)
    {
        normal = _normal;
        d = Vector3DotProduct(_normal, pt);
    }

    /**
     * @brief Constructor of Plane
     *
     * @param _normal The normal of the plane
     * @param _d The distance
     */
    Plane(Vector3 _normal, float _d)
    {
        normal = _normal;
        d = _d;
    }
};

/**
 * @brief A Quad have a referential to have its position and extension to have its size. 
 */
struct Quad
{
    Referential referential;
    Vector3 extension;
};

/**
 * @brief A Sphere have an omega wich represent the sphere, and it have a radius
 */
struct Sphere
{
    Vector3 omega;
    float radius;
};

/**
 * @brief A Cylinder is two point, a ref for position, a radius and height.
 * @remarks The two point are calculated with the referential
 */
struct Cylinder
{
    Referential ref;
    Vector3 pt1;
    Vector3 pt2;
    float radius;
    float height;

    /**
     * @brief Default constructor of Cylinder
     * @remarks If you rotate the referential after this constructor you have to call @see UpdateCylinder to update the two point
     *
     * @param ref The referential of the cylinder
     * @param radius The radius of the cylinder
     * @param height The height of the cylinder
     */
    Cylinder(Referential ref, float radius, float height)
    {
        this->ref = ref;
        this->radius = radius;
        this->height = height;
        this->pt1 = ref.origin;
        this->pt2 = Vector3Add(this->pt1, Vector3Scale(ref.j, height));
    }

    /**
     * @brief It update the position of the two point according to the referential.
     * @remarks After rotated the referential of the cylinder the two point don't have the good coordinates, so this method update coordinates 
     * 
     */
    void UpdateCylinder()
    {
        this->pt1 = ref.origin;
        this->pt2 = Vector3Add(this->pt1, Vector3Scale(ref.j, height));
    }
};

/**
 * @brief A Disk with a referential for the pos and a radius
 */
struct Disk
{
    Referential referential;
    float radius;
};

/**
 * @brief A Capsule with a referential, radius and height
 */
struct Capsule
{
    Referential referential;
    float radius;
    float height;
};

/**
 * @brief A RoundedBox with a referential, extensions and radius.
 */
struct RoundedBox
{
    Referential ref;
    Vector3 extension;
    float radius;
};

/**
 * @brief A RoundedBox with a referential, extensions
 */
struct Box
{
    Referential ref;
    Vector3 extension;
};

#endif

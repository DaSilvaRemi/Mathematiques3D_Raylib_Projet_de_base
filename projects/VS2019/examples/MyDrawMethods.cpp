#include "MyDrawMethods.h"

/**
 * @brief Draw segment
 *
 * @param q The quaternion to turn the segment
 * @param seg The maths segment
 * @param color The color will be display for the shape
 */
void MyDrawSegment(Quaternion q, Segment seg, Color color)
{
    rlPushMatrix();

    //ROTATION
    Vector3 vect;
    float angle;
    QuaternionToAxisAngle(q, &vect, &angle);
    rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

    rlBegin(RL_LINES);

    //Set color and draw a line from PT1 to PT2
    rlColor4ub(color.r, color.g, color.b, color.a);
    rlVertex3f(seg.pt1.x, seg.pt1.y, seg.pt1.z);
    rlVertex3f(seg.pt2.x, seg.pt2.y, seg.pt2.z);

    rlEnd();
    rlPopMatrix();
}

/**
 * @brief Draw a portion of the sphere
 * @remarks end theta set to 2 PI and end Phi set to PI
 *
 * @param q The quaternion to rotate the shape
 * @param sph The maths sphere
 * @param nSegmentsTheta The number of segment display on Theta
 * @param nSegmentsPhi The number of segment display on Phi
 * @param color The color will be display for the shape
 */
void MyDrawSphereEx2(Quaternion q, Sphere sph, int nSegmentsTheta, int nSegmentsPhi, Color color)
{
    MyDrawSpherePortion(q, sph, 0, 2 * PI, 0, PI, nSegmentsTheta, nSegmentsPhi, color);
}

/**
 * @brief Draw a portion of the sphere wires
 * @remarks end theta set to 2 PI and end Phi set to PI
 *
 * @param q The quaternion to rotate the shape
 * @param sph The maths sphere
 * @param nSegmentsTheta The number of segment display on Theta
 * @param nSegmentsPhi The number of segment display on Phi
 * @param color The color will be display for the shape
 */
void MyDrawSphereWiresEx2(Quaternion q, Sphere sph, int nSegmentsTheta, int nSegmentsPhi, Color color)
{
    MyDrawSphereWiresPortion(q, sph, 0, 2 * PI, 0, PI, nSegmentsTheta, nSegmentsPhi, color);
}

/**
 * @brief Draw a portion of the sphere
 *
 * @param q The quaternion to rotate the shape
 * @param sph The maths sphere
 * @param startTheta The start of Theta
 * @param endTheta The end of Theta
 * @param startPhi The start of Phi
 * @param endPhi The end of Phi
 * @param nSegmentsTheta The number of segment display on Theta
 * @param nSegmentsPhi The number of segment display on Phi
 * @param color The color will be display for the shape
 */
void MyDrawSpherePortion(Quaternion q, Sphere sph, float startTheta, float endTheta, float startPhi, float endPhi, int nSegmentsTheta, int nSegmentsPhi, Color color)
{
    if (nSegmentsTheta < 3 || nSegmentsPhi < 2) return;

    //Calculate the offset for Phi and Theta
    float deltaPhi = (endPhi - startPhi) / nSegmentsPhi;
    float deltaTheta = (endTheta - startTheta) / nSegmentsTheta;

    /* Because we have missing segments we need to add a segment in the vertexBuffer */
    std::vector<Vector3> vertexBufferTheta(nSegmentsTheta + 1);
    for (int i = 0; i < nSegmentsTheta + 1; i++)
        vertexBufferTheta[i] = SphericalToCartesian(
            {1, Lerp(startTheta, endTheta, i / static_cast<float>(nSegmentsTheta)), startPhi});

    //Num vertex is the number of the segment will be display multiply by the number of rlVertex3d
    int numVertex = nSegmentsTheta * nSegmentsPhi * 6;
    if (rlCheckBufferLimit(numVertex)) rlglDraw();

    rlPushMatrix();

    // NOTE: Transformation is applied in inverse order (scale -> translate)
    rlTranslatef(sph.omega.x, sph.omega.y, sph.omega.z);

    //ROTATION
    Vector3 vect;
    float angle;
    QuaternionToAxisAngle(q, &vect, &angle);
    rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

    rlScalef(sph.radius, sph.radius, sph.radius);

    rlBegin(RL_TRIANGLES);
    rlColor4ub(color.r, color.g, color.b, color.a);

    float phi = startPhi;
    for (int i = 0; i < nSegmentsPhi; i++)
    {
        //For each segment on phi we calculate the bottom left and the theta
        float theta = startTheta;
        Vector3 tmpBottomLeft = SphericalToCartesian(Spherical{1, theta, phi + deltaPhi});

        for (int j = 0; j < nSegmentsTheta; j++)
        {
            //For each segment on theta we calculate the position of the triangles
            Vector3 topLeft = vertexBufferTheta[j];
            Vector3 bottomLeft = tmpBottomLeft;
            Vector3 topRight = vertexBufferTheta[j + 1];
            Vector3 bottomRight = SphericalToCartesian(Spherical{1, theta + deltaTheta, phi + deltaPhi});

            //Draw the first triangle
            rlVertex3f(bottomLeft.x, bottomLeft.y, bottomLeft.z);
            rlVertex3f(topRight.x, topRight.y, topRight.z);
            rlVertex3f(topLeft.x, topLeft.y, topLeft.z);

            //Draw the second triangle
            rlVertex3f(bottomLeft.x, bottomLeft.y, bottomLeft.z);
            rlVertex3f(bottomRight.x, bottomRight.y, bottomRight.z);
            rlVertex3f(topRight.x, topRight.y, topRight.z);

            //Set the offset in the theta
            theta += deltaTheta;

            //Inversed the position of the bottom left and the bottomRight for each segment
            vertexBufferTheta[j] = tmpBottomLeft;
            tmpBottomLeft = bottomRight;
        }

        //We calculate the last buffer of theta for each segment on theta
        vertexBufferTheta[vertexBufferTheta.size() - 1] = SphericalToCartesian(Spherical{
            1, theta + deltaTheta, phi + deltaPhi
        });

        //Set the offset to phi
        phi += deltaPhi;
    }

    rlEnd();
    rlPopMatrix();
}

/**
 * @brief Draw a portion of the sphere wires
 *
 * @param q The quaternion to rotate the shape
 * @param sph The maths sphere
 * @param startTheta The start of Theta
 * @param endTheta The end of Theta
 * @param startPhi The start of Phi
 * @param endPhi The end of Phi
 * @param nSegmentsTheta The number of segment display on Theta
 * @param nSegmentsPhi The number of segment display on Phi
 * @param color The color will be display for the shape
 */
void MyDrawSphereWiresPortion(Quaternion q, Sphere sph, float startTheta, float endTheta, float startPhi, float endPhi, int nSegmentsTheta, int nSegmentsPhi, Color color)
{
    if (nSegmentsTheta < 3 || nSegmentsPhi < 2) return;

    //Calculate the offset for Phi and Theta
    float deltaPhi = (endPhi - startPhi) / nSegmentsPhi;
    float deltaTheta = (endTheta - startTheta) / nSegmentsTheta;

    /* Because we have missing segments we need to add a segment in the vertexBuffer */
    std::vector<Vector3> vertexBufferTheta(nSegmentsTheta + 1);
    for (size_t i = 0; i < nSegmentsTheta + 1; i++)
        vertexBufferTheta[i] = SphericalToCartesian(
            {1, Lerp(startTheta, endTheta, i / static_cast<float>(nSegmentsTheta)), startPhi});

    //Fill the buffer with points
    std::fill(vertexBufferTheta.begin(), vertexBufferTheta.end(), Vector3{0, 1, 0});

    //Num vertex is the number of the segment will be display multiply by the number of rlVertex3d
    int numVertex = nSegmentsTheta * nSegmentsPhi * 4;
    if (rlCheckBufferLimit(numVertex)) rlglDraw();

    rlPushMatrix();
    // NOTE: Transformation is applied in inverse order (scale -> translate)
    rlTranslatef(sph.omega.x, sph.omega.y, sph.omega.z);

    //ROTATION
    Vector3 vect;
    float angle;
    QuaternionToAxisAngle(q, &vect, &angle);
    rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

    rlScalef(sph.radius, sph.radius, sph.radius);

    rlBegin(RL_LINES);
    rlColor4ub(color.r, color.g, color.b, color.a);

    float phi = startPhi;
    for (int i = 0; i < nSegmentsPhi; i++)
    {
        //For each segment on phi we calculate the theta
        float theta = startTheta;

        for (int j = 0; j < nSegmentsTheta; j++)
        {
            Vector3 topLeft = SphericalToCartesian(Spherical{1, theta, phi}); // vertexBufferTheta[j];
            Vector3 bottomLeft = SphericalToCartesian(Spherical{1, theta, phi + deltaPhi}); //tmpBottomLeft;
            Vector3 topRight = SphericalToCartesian(Spherical{1, theta + deltaTheta, phi}); //vertexBufferTheta[j + 1];

            //Draw the horizontal line
            rlVertex3f(topLeft.x, topLeft.y, topLeft.z);
            rlVertex3f(topRight.x, topRight.y, topRight.z);

            //Draw tbe vertical line
            rlVertex3f(topLeft.x, topLeft.y, topLeft.z);
            rlVertex3f(bottomLeft.x, bottomLeft.y, bottomLeft.z);

            theta += deltaTheta;
            vertexBufferTheta[j] = bottomLeft;
        }

        Vector3 topRight = SphericalToCartesian(Spherical{1, theta + deltaTheta, phi}); //vertexBufferTheta[j + 1];
        Vector3 bottomRight = SphericalToCartesian(Spherical{1, endTheta, phi + deltaPhi});

        //Draw the vertical right segment
        rlVertex3f(topRight.x, topRight.y, topRight.z);
        rlVertex3f(bottomRight.x, bottomRight.y, bottomRight.z);

        //We calculate the last buffer of theta for each segment on theta
        vertexBufferTheta[vertexBufferTheta.size() - 1] = vertexBufferTheta[0];
        phi += deltaPhi;
    }
    rlEnd();
    rlPopMatrix();
}

/**
 * @brief Draw the quad
 * @remarks The quaternion need to be set in the quad referential and the extension are defined in our case like size
 *
 * @param quad The quad
 * @param color The color will be display for the shape
 */
void MyDrawQuad2(Quad quad, Color color)
{
    //Calculate the cordonates of the each point
    float zNeg = -quad.extension.z / 2;
    float z = quad.extension.z / 2;
    float xNeg = -quad.extension.x / 2;
    float x = quad.extension.x / 2;

    rlPushMatrix();
    // NOTE: Transformation is applied in inverse order (scale -> translate)
    rlTranslatef(quad.referential.origin.x, quad.referential.origin.y, quad.referential.origin.z);

    //ROTATION
    Vector3 vect;
    float angle;
    QuaternionToAxisAngle(quad.referential.q, &vect, &angle);
    rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

    rlBegin(RL_TRIANGLES);
    rlColor4ub(color.r, color.g, color.b, color.a);

    //Left Triangle
    rlVertex3f(x, 0, zNeg);
    rlVertex3f(xNeg, 0, zNeg);
    rlVertex3f(x, 0, z);

    //Right Triangle
    rlVertex3f(xNeg, 0, zNeg);
    rlVertex3f(xNeg, 0, z);
    rlVertex3f(x, 0, z);

    rlEnd();
    rlPopMatrix();
}

/**
 * @brief Draw the quad wires
 * @remarks The quaternion need to be set in the quad referential and the extension are defined in our case like size
 *
 * @param quad The quad
 * @param color The color will be display for the shape
 */
void MyDrawQuadWire2(Quad quad, Color color)
{
    //Calculate the cordonates of the each point
    float zNeg = -quad.extension.z / 2;
    float z = quad.extension.z / 2;
    float xNeg = -quad.extension.x / 2;
    float x = quad.extension.x / 2;

    rlPushMatrix();
    // NOTE: Transformation is applied in inverse order (scale -> translate)
    rlTranslatef(quad.referential.origin.x, quad.referential.origin.y, quad.referential.origin.z);
    //ROTATION
    Vector3 vect;
    float angle;
    QuaternionToAxisAngle(quad.referential.q, &vect, &angle);
    rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

    rlBegin(RL_LINES);
    rlColor4ub(color.r, color.g, color.b, color.a);

    //Left Line
    rlVertex3f(xNeg, 0, zNeg);
    rlVertex3f(xNeg, 0, z);

    //Right Line
    rlVertex3f(x, 0, zNeg);
    rlVertex3f(x, 0, z);

    //Top Line
    rlVertex3f(xNeg, 0, z);
    rlVertex3f(x, 0, z);

    //Bottom Line
    rlVertex3f(xNeg, 0, zNeg);
    rlVertex3f(x, 0, zNeg);

    //Cross Line
    rlVertex3f(xNeg, 0, zNeg);
    rlVertex3f(x, 0, z);

    rlEnd();
    rlPopMatrix();
}

/**
 * @brief Draw the cylinder
 * @remarks The quaternion need to be set in the referential and the end of theta set to 2 PI
 *
 * @param cyl The maths cylinder
 * @param nSegmentsTheta The number of segment display on Theta
 * @param drawCaps If draw the "caps" / disk to mean it's finite cylinder or not
 * @param color The color will be display for the shape
 */
void MyDrawCylinder(Cylinder cyl, int nSegmentsTheta, bool drawCaps, Color color)
{
    MyDrawCylinderPortion(cyl, 0, 2 * PI, nSegmentsTheta, drawCaps, color);
}

/**
 * @brief Draw the cylinder wires
 * @remarks The quaternion need to be set in the referential and the end of theta set to 2 PI
 *
 * @param cyl The maths cylinder
 * @param nSegmentsTheta The number of segment display on Theta
 * @param drawCaps If draw the "caps" / disk to mean it's finite cylinder or not
 * @param color The color will be display for the shape
 */
void MyDrawCylinderWires(Cylinder cyl, int nSegmentsTheta, bool drawCaps, Color color)
{
    MyDrawCylinderWiresPortion(cyl, 0, 2 * PI, nSegmentsTheta, drawCaps, color);
}

/**
 * @brief Draw a portion of the cylinder
 * @remarks The quaternion need to be set in the referential
 *
 * @param cyl The maths cylinder
 * @param startTheta The start of Theta
 * @param endTheta The end of Theta
 * @param nSegmentsTheta The number of segment display on Theta
 * @param drawCaps If draw the "caps" / disk to mean it's finite cylinder or not
 * @param color The color will be display for the shape
 */
void MyDrawCylinderPortion(Cylinder cyl, float startTheta, float endTheta, int nSegmentsTheta, bool drawCaps, Color color)
{
    if (nSegmentsTheta < 3) return;

    //Num vertex is the number of the segment will be display multiply by the number of rlVertex3d
    int numVertex = nSegmentsTheta * 6;
    if (rlCheckBufferLimit(numVertex)) rlglDraw();

    //Make the push of the matrix and translate the shape
    rlPushMatrix();
    rlTranslatef(cyl.pt1.x, cyl.pt1.y, cyl.pt1.z);

    //Calculate PQ to have the "height" of the cylinder
    Vector3 PQ = Vector3Subtract(cyl.pt2, cyl.pt1);

    //Rotate and scale the shape
    Vector3 vect;
    float angle;
    QuaternionToAxisAngle(cyl.ref.q, &vect, &angle);
    rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
    rlScalef(cyl.radius, Vector3Length(PQ), cyl.radius);

    rlBegin(RL_TRIANGLES);
    rlColor4ub(color.r, color.g, color.b, color.a);

    //Calculate the offset of the theta
    float deltaTheta = (endTheta - startTheta) / nSegmentsTheta;
    float theta = startTheta;

    //Set the bottom left to change it after
    Vector3 tmpBottomLeft = CylindricalToCartesian(Cylindrical{1, theta, 0});

    for (int i = 0; i < nSegmentsTheta; i++)
    {
        //For each rectangle will be draw, we calculate the bottom Left/Right and top Left/Right of the current rectangle
        Vector3 bottomLeft = tmpBottomLeft;
        Vector3 topLeft = {bottomLeft.x, 1, bottomLeft.z};
        Vector3 bottomRight = CylindricalToCartesian({1, theta + deltaTheta, 0});
        Vector3 topRight = {bottomRight.x, 1, bottomRight.z};

        //We draw the left triangle
        rlVertex3f(bottomLeft.x, bottomLeft.y, bottomLeft.z);
        rlVertex3f(topRight.x, topRight.y, topRight.z);
        rlVertex3f(topLeft.x, topLeft.y, topLeft.z);

        //We draw the right rectangle
        rlVertex3f(bottomLeft.x, bottomLeft.y, bottomLeft.z);
        rlVertex3f(bottomRight.x, bottomRight.y, bottomRight.z);
        rlVertex3f(topRight.x, topRight.y, topRight.z);

        //Set the offet to the theta and set the new tmpBottomLeft
        theta += deltaTheta;
        tmpBottomLeft = bottomRight; //For each rectangle the bottom left and the bottom right are inversed.
    }

    //If drawCaps it's set to true we draw the disk to mean it's a finite cylinder
    if (drawCaps)
    {
        Quaternion qX = QuaternionFromAxisAngle({1, 0, 0}, PI);
        Quaternion qY = QuaternionFromAxisAngle({0, 1, 0}, startTheta);

        //Draw the bottom disk (set in the P point)
        MyDrawDiskPortion({Referential({0, 0, 0}, QuaternionMultiply(qX, qY)), 1}, startTheta, endTheta, nSegmentsTheta,
                          color);
        //Draw the top disk (set in the Q point)
        MyDrawDiskPortion({Referential({0, 1, 0}, QuaternionIdentity()), 1}, startTheta, endTheta, nSegmentsTheta,
                          color);
    }

    rlEnd();
    rlPopMatrix();
}

/**
 * @brief Draw a portion of the cylinder wires
 * @remarks If drawCaps set to false the number of segment will be higher, because we need to draw the missed corner wires
 * @remarks The quaternion need to be set in the referential
 *
 * @param cyl The maths cylinder
 * @param startTheta The start of Theta
 * @param endTheta The end of Theta
 * @param nSegmentsTheta The number of segment display on Theta
 * @param drawCaps If draw the "caps" / disk to mean it's finite cylinder or not
 * @param color The color will be display for the shape
 */
void MyDrawCylinderWiresPortion(Cylinder cyl, float startTheta, float endTheta, int nSegmentsTheta, bool drawCaps, Color color)
{
    if (nSegmentsTheta < 3) return;

    /* Num vertex is the number of the segment will be display multiply by the number of rlVertex3
     * Because we draw wire in this case we need to check if we draw or not the caps because the number of rlVertex will not be the same.
     */
    int numVertex = nSegmentsTheta * (drawCaps ? 0 : 4) * 2;
    if (rlCheckBufferLimit(numVertex)) rlglDraw();

    //Make the push of the matrix and translate the shape
    rlPushMatrix();
    rlTranslatef(cyl.pt1.x, cyl.pt1.y, cyl.pt1.z);

    //Calculate PQ to have the "height" of the cylinder
    Vector3 PQ = Vector3Subtract(cyl.pt2, cyl.pt1);

    //Rotate and scale the shape
    Vector3 vect;
    float angle;
    QuaternionToAxisAngle(cyl.ref.q, &vect, &angle);
    rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
    rlScalef(cyl.radius, Vector3Length(PQ), cyl.radius); // norme

    rlBegin(RL_LINES);
    rlColor4ub(color.r, color.g, color.b, color.a);

    float deltaTheta = (endTheta - startTheta) / nSegmentsTheta;
    float theta = startTheta;
    Vector3 tmpBottomLeft = CylindricalToCartesian(Cylindrical{1, theta, 0});

    for (int i = 0; i < nSegmentsTheta; i++)
    {
        //For each rectangle will be draw, we calculate the bottom Left/Right and top Left/Right of the current rectangle
        Vector3 bottomLeft = tmpBottomLeft;
        Vector3 topLeft = {bottomLeft.x, 1, bottomLeft.z};
        Vector3 bottomRight = CylindricalToCartesian({1, theta + deltaTheta, 0});
        Vector3 topRight = {bottomRight.x, 1, bottomRight.z};

        //We draw the left line
        rlVertex3f(bottomLeft.x, bottomLeft.y, bottomLeft.z);
        rlVertex3f(topLeft.x, topLeft.y, topLeft.z);

        //Draw caps set to false we need to draw the corner
        if (!drawCaps)
        {
            //We draw the bottom corner
            rlVertex3f(bottomLeft.x, bottomLeft.y, bottomLeft.z);
            rlVertex3f(bottomRight.x, bottomRight.y, bottomRight.z);

            //We draw the top corner
            rlVertex3f(topLeft.x, topLeft.y, topLeft.z);
            rlVertex3f(topRight.x, topRight.y, topRight.z);
        }

        //Set the offet to the theta and set the new tmpBottomLeft
        theta += deltaTheta;
        tmpBottomLeft = bottomRight;
    }

    if (drawCaps)
    {
        Quaternion qX = QuaternionFromAxisAngle({1, 0, 0}, PI);
        Quaternion qY = QuaternionFromAxisAngle({0, 1, 0}, startTheta);

        //Draw the bottom disk (set in the P point)
        MyDrawDiskWiresPortion({Referential({0, 0, 0}, QuaternionMultiply(qX, qY)), 1}, startTheta, endTheta,
                               nSegmentsTheta, color);
        //Draw the top disk (set in the Q point)
        MyDrawDiskWiresPortion({Referential({0, 1, 0}, QuaternionIdentity()), 1}, startTheta, endTheta, nSegmentsTheta,
                               color);
    }

    rlEnd();
    rlPopMatrix();
}

/**
 * @brief Draw Disk
 * @remarks The end the set to 2 PI and The quaternion need to be set in the referential
 *
 * @param disk The maths disk
 * @param nSegmentsTheta The number of segment display on Theta
 * @param color The color will be display for the shape
 */
void MyDrawDisk(Disk disk, int nSegmentsTheta, Color color)
{
    MyDrawDiskPortion(disk, 0, 2 * PI, nSegmentsTheta, color);
}

/**
 * @brief Draw Disk wires
 * @remarks The end the set to 2 PI and The quaternion need to be set in the referential
 *
 * @param disk The maths disk
 * @param nSegmentsTheta The number of segment display on Theta
 * @param color The color will be display for the shape
 */
void MyDrawDiskWires(Disk disk, int nSegmentsTheta, Color color)
{
    MyDrawDiskWiresPortion(disk, 0, 2 * PI, nSegmentsTheta, color);
}

/**
 * @brief Draw Disk portion
 * @remarks The quaternion need to be set in the referential
 *
 * @param disk The maths disk
 * @param startTheta The start of Theta
 * @param endTheta The end of Theta
 * @param nSegmentsTheta The number of segment display on Theta
 * @param color The color will be display for the shape
 */
void MyDrawDiskPortion(Disk disk, float startTheta, float endTheta, int nSegmentsTheta, Color color)
{
    if (nSegmentsTheta < 3) return;
    //Num vertex is the number of the segment will be display multiply by the number of rlVertex3d
    int numVertex = nSegmentsTheta * 3;
    if (rlCheckBufferLimit(numVertex)) rlglDraw();

    rlPushMatrix();
    rlTranslatef(disk.referential.origin.x, disk.referential.origin.y, disk.referential.origin.z);

    //ROTATION
    Vector3 vect;
    float angle;
    QuaternionToAxisAngle(disk.referential.q, &vect, &angle);
    rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
    rlScalef(disk.radius, 0, disk.radius);

    rlBegin(RL_TRIANGLES);
    rlColor4ub(color.r, color.g, color.b, color.a);

    //Calculate the offset of the theta
    float deltaTheta = (endTheta - startTheta) / nSegmentsTheta;
    float theta = startTheta;

    //Set the bottom left to change it after
    Vector3 tmpBottomLeft = CylindricalToCartesian(Cylindrical{1, theta, 0});

    for (int i = 0; i < nSegmentsTheta; i++)
    {
        //For each segment we calculate the bottom left and right
        Vector3 bottomLeft = tmpBottomLeft;
        Vector3 bottomRight = CylindricalToCartesian({1, theta + deltaTheta, 0});

        //Draw a triangle
        rlVertex3f(0, 0, 0);
        rlVertex3f(bottomLeft.x, bottomLeft.y, bottomLeft.z);
        rlVertex3f(bottomRight.x, bottomRight.y, bottomRight.z);

        //Set the offet to the theta and inversed bottom left and bottom right
        theta += deltaTheta;
        tmpBottomLeft = bottomRight;
    }

    rlEnd();
    rlPopMatrix();
}

/**
 * @brief Draw Disk wires portion
 * @remarks The quaternion need to be set in the referential
 *
 * @param disk The maths disk
 * @param startTheta The start of Theta
 * @param endTheta The end of Theta
 * @param nSegmentsTheta The number of segment display on Theta
 * @param color The color will be display for the shape
 */
void MyDrawDiskWiresPortion(Disk disk, float startTheta, float endTheta, int nSegmentsTheta, Color color)
{
    if (nSegmentsTheta < 3) return;

    rlPushMatrix();
    rlTranslatef(disk.referential.origin.x, disk.referential.origin.y, disk.referential.origin.z);

    //ROTATION
    Vector3 vect;
    float angle;
    QuaternionToAxisAngle(disk.referential.q, &vect, &angle);
    rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
    rlScalef(disk.radius, 0, disk.radius);

    rlBegin(RL_LINES);
    rlColor4ub(color.r, color.g, color.b, color.a);

    //Calculate the offset of the theta
    float deltaTheta = (endTheta - startTheta) / nSegmentsTheta;
    float theta = startTheta;

    //Set the bottom left to change it after
    Vector3 tmpBottomLeft = CylindricalToCartesian(Cylindrical{1, theta, 0});

    for (int i = 0; i < nSegmentsTheta; i++)
    {
        //For each segment we calculate the bottom left and right
        Vector3 bottomLeft = tmpBottomLeft;
        Vector3 bottomRight = CylindricalToCartesian({1, theta + deltaTheta, 0});

        //Draw line from center to the bottom left
        rlVertex3f(0, 0, 0);
        rlVertex3f(bottomLeft.x, bottomLeft.y, bottomLeft.z);

        //Draw the corner of the disk
        rlVertex3f(bottomRight.x, bottomRight.y, bottomRight.z);
        rlVertex3f(bottomLeft.x, bottomLeft.y, bottomLeft.z);

        //Set the offet to the theta and inversed bottom left and bottom right
        theta += deltaTheta;
        tmpBottomLeft = bottomRight;
    }

    //Draw line from center to the bottom right after all segment for finish the disk
    Vector3 bottomRight = tmpBottomLeft;
    rlVertex3f(0, 0, 0);
    rlVertex3f(bottomRight.x, bottomRight.y, bottomRight.z);

    rlEnd();
    rlPopMatrix();
}

/**
 * @brief Draw Capsule
 * @remarks The quaternion need to be set in the referential
 *
 * @param capsule The maths capsule
 * @param color The color will be display for the shape
 */
void MyDrawCapsule(Capsule capsule, Color color)
{
    rlPushMatrix();
    rlTranslatef(capsule.referential.origin.x, capsule.referential.origin.y, capsule.referential.origin.z);

    //ROTATION
    Vector3 vect;
    float angle;
    QuaternionToAxisAngle(capsule.referential.q, &vect, &angle);
    rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

    //Calculate the position of P and Q point
    Vector3 P = {0, 0, 0};
    Vector3 Q = {0, capsule.height, 0};

    //Calculate the quaternions
    Quaternion qTop = QuaternionFromAxisAngle({0, 0, 1}, 0.5 * PI);
    Quaternion qBottom = QuaternionFromAxisAngle({0, 0, 1}, -0.5 * PI);
    Quaternion qIdentity = QuaternionIdentity();

    //Set the referential from P and set the capsule with ref, the radius and height
    auto ref = Referential(P, qIdentity);
    auto cylinder = Cylinder(ref, capsule.radius, capsule.height);

    //Set the top and the bottom sphere
    Sphere sphereTop = {Q, capsule.radius};
    Sphere sphereBottom = {P, capsule.radius};

    //Draw complete cylinder with caps
    MyDrawCylinder(cylinder, 25, true, color);
    //Draw half sphere (From 0 to PI)
    MyDrawSpherePortion(qTop, sphereTop, 0, PI, 0, PI, 30, 30, color);
    MyDrawSpherePortion(qBottom, sphereBottom, 0, PI, 0, PI, 25, 25, color);

    rlPopMatrix();
}

/**
 * @brief Draw Capsule wires
 * @remarks The quaternion need to be set in the referential
 *
 * @param capsule The maths capsule
 * @param color The color will be display for the shape
 */
void MyDrawCapsuleWires(Capsule capsule, Color color)
{
    rlPushMatrix();
    rlTranslatef(capsule.referential.origin.x, capsule.referential.origin.y, capsule.referential.origin.z);

    //ROTATION
    Vector3 vect;
    float angle;
    QuaternionToAxisAngle(capsule.referential.q, &vect, &angle);
    rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

    //Calculate the position of P and Q point
    Vector3 P = {0, 0, 0};
    Vector3 Q = {0, capsule.height, 0};

    Quaternion qUp = QuaternionFromAxisAngle({0, 0, 1}, 0.5 * PI);
    Quaternion qDown = QuaternionFromAxisAngle({0, 0, 1}, 1.5 * PI);
    Quaternion qIdentity = QuaternionIdentity();

    //Set the referential from P and set the capsule with ref, the radius and height
    auto ref = Referential(P, qIdentity);
    auto cylinder = Cylinder(ref, capsule.radius, capsule.height);

    //Set the top and the bottom sphere
    Sphere sphereUp = {Q, capsule.radius};
    Sphere sphereDown = {P, capsule.radius};

    //Draw complete cylinder with caps
    MyDrawCylinderWires(cylinder, 25, true, color);
    //Draw half sphere (From 0 to PI)
    MyDrawSphereWiresPortion(qUp, sphereUp, 0, PI, 0, PI, 30, 30, color);
    MyDrawSphereWiresPortion(qDown, sphereDown, 0, PI, 0, PI, 25, 25, color);

    rlPopMatrix();
}

/**
 * @brief Draw Rounded Box
 * @remarks The quaternion need to be set in the referential, a rounded box it's 4 quad and 8 cylinder 
 *
 * @param roundedBox The maths roundedBox
 * @param color The color will be display for the shape
 */
void MyDrawRoundedBoxV2(RoundedBox roundedBox, Color color)
{
    rlPushMatrix();

    //ROTATION
    Vector3 vect;
    float angle;
    QuaternionToAxisAngle(roundedBox.ref.q, &vect, &angle);
    rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

    // The origin of the RoundedBox is originally set in share sphere set in the intersection of third priamry capsules
    // We need to translate it of one negative extension in X and Y, and one positive extension in Z
    // Example : We have a roundedBox.extension which is the length of the cylinder of the caps (so x2 to real extension)
    Vector3 posRef = Vector3Add(roundedBox.ref.origin, {-roundedBox.extension.x / 2, -roundedBox.extension.y / 2, roundedBox.extension.z / 2 });

    // We take 3 caps to our referential and move all the shapes from this referential
    //Move to the translation referential and resize to the Z to have the left corner (Z local ref)
    Quaternion qLeft = QuaternionFromAxisAngle({1, 0, 0}, -PI * 0.5f);
    Capsule capsLeftBottom = {Referential(posRef, qLeft), roundedBox.radius, roundedBox.extension.z};
    MyDrawCapsule(capsLeftBottom, color);

    //Move to the translation referential and resize to the X  to have the front corner (X local ref)
    Quaternion qFront = QuaternionFromAxisAngle({0, 0, 1}, -PI * 0.5f);
    Capsule capsFrontBottom = {Referential(posRef, qFront), roundedBox.radius, roundedBox.extension.x};
    MyDrawCapsule(capsFrontBottom, color);

    //Move to the translation referential and resize to the Y to have the front up corner (Y local ref)
    Quaternion qUp = QuaternionIdentity();
    Capsule capsFrontLeft = {Referential(posRef, qUp), roundedBox.radius, roundedBox.extension.y};
    MyDrawCapsule(capsFrontLeft, color);
    // end referential of capsules

    // All other capsules calculated from the translated referential
    Capsule capsFrontTop = {
        Referential(Vector3Add(posRef, {0, roundedBox.extension.y, 0}), qFront), roundedBox.radius,
        roundedBox.extension.x
    };
    MyDrawCapsule(capsFrontTop, color);

    Capsule capsFrontRight = {
        Referential(Vector3Add(posRef, {roundedBox.extension.x, 0, 0}), qUp), roundedBox.radius, roundedBox.extension.y
    };
    MyDrawCapsule(capsFrontRight, color);

    Capsule capsRightBottom = {
        Referential(Vector3Add(posRef, {roundedBox.extension.x, 0, 0}), qLeft), roundedBox.radius,
        roundedBox.extension.z
    };
    MyDrawCapsule(capsRightBottom, color);

    Capsule capsRightTop = {
        Referential(Vector3Add(posRef, {roundedBox.extension.x, roundedBox.extension.y, 0}), qLeft), roundedBox.radius,
        roundedBox.extension.z
    };
    MyDrawCapsule(capsRightTop, color);

    Capsule capsLeftTop = {
        Referential(Vector3Add(posRef, {0, roundedBox.extension.y, 0}), qLeft), roundedBox.radius,
        roundedBox.extension.z
    };
    MyDrawCapsule(capsLeftTop, color);

    Capsule capsBackBottom = {
        Referential(Vector3Add(posRef, {0, 0, -roundedBox.extension.z}), qFront), roundedBox.radius,
        roundedBox.extension.x
    };
    MyDrawCapsule(capsBackBottom, color);

    Capsule capsBackTop = {
        Referential(Vector3Add(posRef, {0, roundedBox.extension.y, -roundedBox.extension.z}), qFront),
        roundedBox.radius, roundedBox.extension.x
    };
    MyDrawCapsule(capsBackTop, color);

    Capsule capsBackLeft = {
        Referential(Vector3Add(posRef, {0, 0, -roundedBox.extension.z}), qUp), roundedBox.radius, roundedBox.extension.y
    };
    MyDrawCapsule(capsBackLeft, color);

    Capsule capsBackRight = {
        Referential(Vector3Add(posRef, {roundedBox.extension.x, 0, -roundedBox.extension.z}), qUp), roundedBox.radius,
        roundedBox.extension.y
    };
    MyDrawCapsule(capsBackRight, color);
    //End of the capsule

    //All the quad calculate from translated referential.
    //Because of the radius of the caps we need to take part of the radius in the coordinates of the quads

    //The front quad set to be push in Z local axes
    Quaternion qFrontQuad = QuaternionFromAxisAngle({1, 0, 0}, PI * 0.5f);
    Referential refQuadFront = Referential(Vector3Add(posRef, {
                                                   roundedBox.extension.x / 2, roundedBox.extension.y / 2,
                                                   roundedBox.radius
                                               }), qFrontQuad);
    Quad quadFront = {refQuadFront, {roundedBox.extension.x, roundedBox.extension.z, roundedBox.extension.y}};
    MyDrawQuad2(quadFront, color);

    //The front quad set to be push in Z local axes
    Quaternion qBackQuad = QuaternionFromAxisAngle({1, 0, 0}, -PI * 0.5f);
    Referential refQuadBack = Referential(Vector3Add(posRef, {
                                                  roundedBox.extension.x / 2, roundedBox.extension.y / 2,
                                                  -(roundedBox.extension.z + roundedBox.radius)
                                              }), qBackQuad);
    Quad quadBack = {refQuadBack, {roundedBox.extension.x, roundedBox.extension.z, roundedBox.extension.y}};
    MyDrawQuad2(quadBack, color);

    //The front quad set to be push in X local axes
    Quaternion qRightQuad = QuaternionFromAxisAngle({0, 0, 1}, -PI * 0.5f);
    Referential refQuadRight = Referential(Vector3Add(posRef,
                                            {
                                                   roundedBox.extension.x + roundedBox.radius,
                                                   roundedBox.extension.y / 2, -roundedBox.extension.z / 2
                                               }), qRightQuad);
    Quad quadRight = {refQuadRight, {roundedBox.extension.y, roundedBox.extension.x, roundedBox.extension.z}};
    MyDrawQuad2(quadRight, color);

    //The front quad set to be push in X local axes
    Quaternion qLeftQuad = QuaternionFromAxisAngle({0, 0, 1}, PI * 0.5f);
    Referential refQuadLeft = Referential(Vector3Add(posRef,
                                            {
                                                  -roundedBox.radius, roundedBox.extension.y / 2,
                                                  -roundedBox.extension.z / 2
                                              }), qLeftQuad);
    Quad quadLeft = {refQuadLeft, {roundedBox.extension.y, roundedBox.extension.x, roundedBox.extension.z}};
    MyDrawQuad2(quadLeft, color);

    //The front quad set to be push in Y local axes
    Quaternion qTopQuad = QuaternionIdentity();
    Referential refQuadTop = Referential(Vector3Add(posRef,
                                            {
                                                 roundedBox.extension.x / 2,
                                                 roundedBox.extension.y + roundedBox.radius,
                                                 -roundedBox.extension.z / 2
                                             }), qTopQuad);
    Quad quadTop = {refQuadTop, {roundedBox.extension.x, roundedBox.extension.y, roundedBox.extension.z}};
    MyDrawQuad2(quadTop, color);

    //The front quad set to be push in Y local axes
    Quaternion qBottomQuad = QuaternionFromAxisAngle({0, 0, 1}, PI);
    Referential refQuadBottom = Referential(Vector3Add(posRef,
                                                {
                                                    roundedBox.extension.x / 2, -roundedBox.radius,
                                                    -roundedBox.extension.z / 2
                                                }), qBottomQuad);
    Quad quadBottom = {refQuadBottom, {roundedBox.extension.x, roundedBox.extension.y, roundedBox.extension.z}};
    MyDrawQuad2(quadBottom, color);

    rlPopMatrix();
}

/**
 * @brief Draw Rounded Box wires
 * @remarks The quaternion need to be set in the referential, a rounded box it's 4 quad and 8 cylinder 
 *
 * @param roundedBox The maths roundedBox
 * @param color The color will be display for the shape
 */
void MyDrawRoundBoxWiresV2(RoundedBox roundedBox, Color color)
{
    rlPushMatrix();

    //ROTATION
    Vector3 vect;
    float angle;
    QuaternionToAxisAngle(roundedBox.ref.q, &vect, &angle);
    rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

    /* For drawing our RoundedBox we don't take the center referential but we take a translated referential because we start our draw from a sphere referential */
    Vector3 posRef = Vector3Add(roundedBox.ref.origin, {
                                    -roundedBox.extension.x / 2, -roundedBox.extension.y / 2, roundedBox.extension.z / 2
                                });

    // We take 3 caps to our referential and move all the shapes from this referential
    //Move to the translation referential and resize to the Z to have the left corner (Z local ref)
    Quaternion qLeft = QuaternionFromAxisAngle({1, 0, 0}, -PI * 0.5f);
    Capsule capsLeftBottom = {Referential(posRef, qLeft), roundedBox.radius, roundedBox.extension.z};
    MyDrawCapsule(capsLeftBottom, color);

    //Move to the translation referential and resize to the X  to have the front corner (X local ref)
    Quaternion qFront = QuaternionFromAxisAngle({0, 0, 1}, -PI * 0.5f);
    Capsule capsFrontBottom = {Referential(posRef, qFront), roundedBox.radius, roundedBox.extension.x};
    MyDrawCapsule(capsFrontBottom, color);

    //Move to the translation referential and resize to the Y to have the front up corner (Y local ref)
    Quaternion qUp = QuaternionIdentity();
    Capsule capsFrontLeft = {Referential(posRef, qUp), roundedBox.radius, roundedBox.extension.y};
    MyDrawCapsuleWires(capsFrontLeft, color);
    // end referential of capsules

    // All other capsules calculated from the translated referential
    /* We start on the ref plus the extension Y to go top
     * We size the caps with X extension
     * All other capsule are the same working we translate on one axes and we resize it
     */
    Capsule capsFrontTop = {
        Referential(Vector3Add(posRef, {0, roundedBox.extension.y, 0}), qFront), roundedBox.radius,
        roundedBox.extension.x
    };
    MyDrawCapsuleWires(capsFrontTop, color);

    Capsule capsFrontRight = {
        Referential(Vector3Add(posRef, {roundedBox.extension.x, 0, 0}), qUp), roundedBox.radius, roundedBox.extension.y
    };
    MyDrawCapsuleWires(capsFrontRight, color);

    Capsule capsRightBottom = {
        Referential(Vector3Add(posRef, {roundedBox.extension.x, 0, 0}), qLeft), roundedBox.radius,
        roundedBox.extension.z
    };
    MyDrawCapsuleWires(capsRightBottom, color);

    Capsule capsRightTop = {
        Referential(Vector3Add(posRef, {roundedBox.extension.x, roundedBox.extension.y, 0}), qLeft), roundedBox.radius,
        roundedBox.extension.z
    };
    MyDrawCapsuleWires(capsRightTop, color);

    Capsule capsLeftTop = {
        Referential(Vector3Add(posRef, {0, roundedBox.extension.y, 0}), qLeft), roundedBox.radius,
        roundedBox.extension.z
    };
    MyDrawCapsuleWires(capsLeftTop, color);

    Capsule capsBackBottom = {
        Referential(Vector3Add(posRef, {0, 0, -roundedBox.extension.z}), qFront), roundedBox.radius,
        roundedBox.extension.x
    };
    MyDrawCapsuleWires(capsBackBottom, color);

    Capsule capsBackTop = {
        Referential(Vector3Add(posRef, {0, roundedBox.extension.y, -roundedBox.extension.z}), qFront),
        roundedBox.radius, roundedBox.extension.x
    };
    MyDrawCapsuleWires(capsBackTop, color);

    Capsule capsBackLeft = {
        Referential(Vector3Add(posRef, {0, 0, -roundedBox.extension.z}), qUp), roundedBox.radius, roundedBox.extension.y
    };
    MyDrawCapsuleWires(capsBackLeft, color);

    Capsule capsBackRight = {
        Referential(Vector3Add(posRef, {roundedBox.extension.x, 0, -roundedBox.extension.z}), qUp), roundedBox.radius,
        roundedBox.extension.y
    };
    MyDrawCapsuleWires(capsBackRight, color);
    //End of the capsule

    //All the quad calculate from translated referential.
    //Because of the radius of the caps we need to take part of the radius in the coordinates of the quads

    /* We start on the ref plus the extension /2 to have exactly into cylinder and we set the radius to the set to out the quad
    * We size the caps with X extension, and Z to the Y and Y to Z because we turn to PI / 2
    * All other quad are the same working we translate it and resize it according local referential
    */
    Quaternion qFrontQuad = QuaternionFromAxisAngle({1, 0, 0}, PI * 0.5f);
    Referential refQuadFront = Referential(Vector3Add(posRef, {
                                                   roundedBox.extension.x / 2, roundedBox.extension.y / 2,
                                                   roundedBox.radius
                                               }), qFrontQuad);
    Quad quadFront = {refQuadFront, {roundedBox.extension.x, roundedBox.extension.z, roundedBox.extension.y}};
    MyDrawQuadWire2(quadFront, color);

    //The front quad set to be push in Z local axes
    Quaternion qBackQuad = QuaternionFromAxisAngle({1, 0, 0}, -PI * 0.5f);
    Referential refQuadBack = Referential(Vector3Add(posRef, {
                                                  roundedBox.extension.x / 2, roundedBox.extension.y / 2,
                                                  -(roundedBox.extension.z + roundedBox.radius)
                                              }), qBackQuad);
    Quad quadBack = {refQuadBack, {roundedBox.extension.x, roundedBox.extension.z, roundedBox.extension.y}};
    MyDrawQuadWire2(quadBack, color);

    //The front quad set to be push in X local axes
    Quaternion qRightQuad = QuaternionFromAxisAngle({0, 0, 1}, -PI * 0.5f);
    Referential refQuadRight = Referential(Vector3Add(posRef, {
                                                   roundedBox.extension.x + roundedBox.radius,
                                                   roundedBox.extension.y / 2, -roundedBox.extension.z / 2
                                               }), qRightQuad);
    Quad quadRight = {refQuadRight, {roundedBox.extension.y, roundedBox.extension.x, roundedBox.extension.z}};
    MyDrawQuadWire2(quadRight, color);

    //The front quad set to be push in X local axes
    Quaternion qLeftQuad = QuaternionFromAxisAngle({0, 0, 1}, PI * 0.5f);
    Referential refQuadLeft = Referential(Vector3Add(posRef, {
                                                  -roundedBox.radius, roundedBox.extension.y / 2,
                                                  -roundedBox.extension.z / 2
                                              }), qLeftQuad);
    Quad quadLeft = {refQuadLeft, {roundedBox.extension.y, roundedBox.extension.x, roundedBox.extension.z}};
    MyDrawQuadWire2(quadLeft, color);

    //The front quad set to be push in Y local axes
    Quaternion qTopQuad = QuaternionIdentity();
    Referential refQuadTop = Referential(Vector3Add(posRef, {
                                                 roundedBox.extension.x / 2,
                                                 roundedBox.extension.y + roundedBox.radius,
                                                 -roundedBox.extension.z / 2
                                             }), qTopQuad);
    Quad quadTop = {refQuadTop, {roundedBox.extension.x, roundedBox.extension.y, roundedBox.extension.z}};
    MyDrawQuadWire2(quadTop, color);

    //The front quad set to be push in Y local axes
    Quaternion qBottomQuad = QuaternionFromAxisAngle({0, 0, 1}, PI);
    Referential refQuadBottom = Referential(Vector3Add(posRef, {
                                                    roundedBox.extension.x / 2, -roundedBox.radius,
                                                    -roundedBox.extension.z / 2
                                                }), qBottomQuad);
    Quad quadBottom = {refQuadBottom, {roundedBox.extension.x, roundedBox.extension.y, roundedBox.extension.z}};
    MyDrawQuadWire2(quadBottom, color);

    rlPopMatrix();
}

/**
 * @brief Draw Box
 * @remarks The quaternion need to be set in the referential 
 *
 * @param box The maths Box
 * @param color The color will be display for the shape
 */
void MyDrawBox(Box box, Color color)
{
    rlPushMatrix();

    //ROTATION
    Vector3 vect;
    float angle;
    QuaternionToAxisAngle(box.ref.q, &vect, &angle);
    rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

    //Set all quaternions
    Quaternion qTop = QuaternionIdentity();
    Quaternion qFront = QuaternionFromAxisAngle({0, 0, 1}, PI * -0.5f);
    Quaternion qBack = QuaternionFromAxisAngle({0, 0, 1}, PI * 0.5f);
    Quaternion qLeft = QuaternionFromAxisAngle({1, 0, 0}, PI * 0.5f);
    Quaternion qRight = QuaternionFromAxisAngle({1, 0, 0}, PI * -0.5f);

    Quaternion qX = QuaternionFromAxisAngle({1, 0, 0}, PI);
    Quaternion qY = QuaternionFromAxisAngle({0, 1, 0}, 0);
    Quaternion qBottom = QuaternionMultiply(qX, qY);

    //Set the referential of all quad according to the referential of box and it extension
    Referential referentialQuadTop = Referential(Vector3Add(box.ref.origin, {0, box.extension.y / 2, 0}), qTop);
    Referential referentialQuadFront = Referential(Vector3Add(box.ref.origin, {box.extension.x / 2, 0, 0}), qFront);
    Referential referentialQuadBack = Referential(Vector3Add(box.ref.origin, {-box.extension.x / 2, 0, 0}), qBack);
    Referential referentialQuadLeft = Referential(Vector3Add(box.ref.origin, {0, 0, box.extension.z / 2}), qLeft);
    Referential referentialQuadRight = Referential(Vector3Add(box.ref.origin, {0, 0, -box.extension.z / 2}), qRight);
    Referential referentialQuadBottom = Referential(Vector3Add(box.ref.origin, {0, -box.extension.y / 2, 0}), qBottom);

    // Set all the quad with referential and the extension according to local ref of the quad
    // Example if we turn to 90° or PI/2 around X axes the Z extension become Y extension and Y extension become Z extension
    Quad quadTop = {referentialQuadTop, box.extension};
    Quad quadFront = {referentialQuadFront, {box.extension.y, box.extension.x, box.extension.z}};
    Quad quadBack = {referentialQuadBack, {box.extension.y, box.extension.x, box.extension.z}};
    Quad quadLeft = {referentialQuadLeft, {box.extension.x, box.extension.z, box.extension.y}};
    Quad quadRight = {referentialQuadRight, {box.extension.x, box.extension.z, box.extension.y}};
    Quad quadBottom = {referentialQuadBottom, box.extension};

    //Draw all the quads
    MyDrawQuad2(quadTop, color);
    MyDrawQuad2(quadFront, color);
    MyDrawQuad2(quadBack, color);
    MyDrawQuad2(quadLeft, color);
    MyDrawQuad2(quadRight, color);
    MyDrawQuad2(quadBottom, color);

    rlPopMatrix();
}

/**
 * @brief Draw Box wires
 * @remarks The quaternion need to be set in the referential 
 *
 * @param box The maths Box
 * @param color The color will be display for the shape
 */
void MyDrawBoxWires(Box box, Color color)
{
    rlPushMatrix();

    //ROTATION
    Vector3 vect;
    float angle;
    QuaternionToAxisAngle(box.ref.q, &vect, &angle);
    rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

    //Set all quaternions
    Quaternion qTop = QuaternionIdentity();
    Quaternion qFront = QuaternionFromAxisAngle({0, 0, 1}, PI * -0.5f);
    Quaternion qBack = QuaternionFromAxisAngle({0, 0, 1}, PI * 0.5f);
    Quaternion qLeft = QuaternionFromAxisAngle({1, 0, 0}, PI * 0.5f);
    Quaternion qRight = QuaternionFromAxisAngle({1, 0, 0}, PI * -0.5f);

    Quaternion qX = QuaternionFromAxisAngle({1, 0, 0}, PI);
    Quaternion qY = QuaternionFromAxisAngle({0, 1, 0}, 0);
    Quaternion qBottom = QuaternionMultiply(qX, qY);

    //Set the referential of all quad according to the referential of box and it extension
    Referential referentialQuadTop = Referential(Vector3Add(box.ref.origin, {0, box.extension.y / 2, 0}), qTop);
    Referential referentialQuadFront = Referential(Vector3Add(box.ref.origin, {box.extension.x / 2, 0, 0}), qFront);
    Referential referentialQuadBack = Referential(Vector3Add(box.ref.origin, {-box.extension.x / 2, 0, 0}), qBack);
    Referential referentialQuadLeft = Referential(Vector3Add(box.ref.origin, {0, 0, box.extension.z / 2}), qLeft);
    Referential referentialQuadRight = Referential(Vector3Add(box.ref.origin, {0, 0, -box.extension.z / 2}), qRight);
    Referential referentialQuadBottom = Referential(Vector3Add(box.ref.origin, {0, -box.extension.y / 2, 0}), qBottom);

    // Set all the quad with referential and the extension according to local ref of the quad
    // Example if we turn to 90° or PI/2 around X axes the Z extension become Y extension and Y extension become Z extension
    Quad quadTop = {referentialQuadTop, box.extension};
    Quad quadFront = {referentialQuadFront, {box.extension.y, box.extension.x, box.extension.z}};
    Quad quadBack = {referentialQuadBack, {box.extension.y, box.extension.x, box.extension.z}};
    Quad quadLeft = {referentialQuadLeft, {box.extension.x, box.extension.z, box.extension.y}};
    Quad quadRight = {referentialQuadRight, {box.extension.x, box.extension.z, box.extension.y}};
    Quad quadBottom = {referentialQuadBottom, box.extension};

    //Draw all the quads wires
    MyDrawQuadWire2(quadTop, RED);
    MyDrawQuadWire2(quadFront, BROWN);
    MyDrawQuadWire2(quadBack, PURPLE);
    MyDrawQuadWire2(quadLeft, YELLOW);
    MyDrawQuadWire2(quadRight, BLUE);
    MyDrawQuadWire2(quadBottom, DARKGREEN);

    rlPopMatrix();
}

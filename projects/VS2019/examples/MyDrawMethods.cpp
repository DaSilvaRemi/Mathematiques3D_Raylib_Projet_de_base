#include "MyDrawMethods.h"

void MyDrawSegment(Quaternion q, Segment seg, Color color) {
	rlPushMatrix();

	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

	rlBegin(RL_LINES);
	
	rlColor4ub(color.r, color.g, color.b, color.a);
	rlVertex3f(seg.pt1.x, seg.pt1.y, seg.pt1.z);
	rlVertex3f(seg.pt2.x, seg.pt2.y, seg.pt2.z);
	
	rlEnd();
	rlPopMatrix();
}

void MyDrawCylinderPortion(Quaternion q, Cylinder cyl, float startTheta, float endTheta, int nSegmentsTheta, bool drawCaps, Color color) {
	if (nSegmentsTheta < 3) return;

	int numVertex = nSegmentsTheta * 6;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	rlPushMatrix();
	rlTranslatef(cyl.pt1.x, cyl.pt1.y, cyl.pt1.z);

	//ROTATION
	Vector3 AB = Vector3Subtract(cyl.pt2, cyl.pt1);
	Quaternion qVector = QuaternionFromVector3ToVector3({0, 1, 0}, Vector3Normalize(AB));
	Quaternion qMult = QuaternionMultiply(q, qVector);

	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(qMult, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	rlScalef(cyl.radius, Vector3Length(AB), cyl.radius); // norme

	rlBegin(RL_TRIANGLES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	float deltaTheta = (endTheta - startTheta) / nSegmentsTheta;
	float theta = startTheta;
	Vector3 tmpBottomLeft = CylindricalToCartesian(Cylindrical{1, theta, 0});

	for (int i = 0; i < nSegmentsTheta; i++)
	{
		Vector3 bottomLeft = tmpBottomLeft;
		Vector3 topLeft = { bottomLeft.x, 1, bottomLeft.z};
		Vector3 bottomRight = CylindricalToCartesian({1, theta + deltaTheta, 0});
		Vector3 topRight = {bottomRight.x, 1, bottomRight.z};

		rlVertex3f(bottomLeft.x, bottomLeft.y, bottomLeft.z);
		rlVertex3f(topRight.x, topRight.y, topRight.z);
		rlVertex3f(topLeft.x, topLeft.y, topLeft.z);

		rlVertex3f(bottomLeft.x, bottomLeft.y, bottomLeft.z);
		rlVertex3f(bottomRight.x, bottomRight.y, bottomRight.z);
		rlVertex3f(topRight.x, topRight.y, topRight.z);

		theta = theta + deltaTheta;
		tmpBottomLeft = bottomRight;
	}

	if(drawCaps)
	{
		Quaternion qX = QuaternionFromAxisAngle({1, 0, 0}, PI);
		Quaternion qY = QuaternionFromAxisAngle({0, 1, 0}, startTheta);
		MyDrawDiskPortion(QuaternionMultiply(qX, qY), {0, 0, 0}, 1, startTheta, endTheta, nSegmentsTheta, color);
		MyDrawDiskPortion(QuaternionIdentity(), {0, 1, 0}, 1, startTheta, endTheta, nSegmentsTheta, color);
	}
	
	rlEnd();
	rlPopMatrix();
}

void MyDrawCylinderWiresPortion(Quaternion q, Cylinder cyl, float startTheta, float endTheta, int nSegmentsTheta, bool drawCaps, Color color) {
	if (nSegmentsTheta < 3) return;

	int numVertex = nSegmentsTheta * (drawCaps ? 0 : 4) * 2;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	rlPushMatrix();
	rlTranslatef(cyl.pt1.x, cyl.pt1.y, cyl.pt1.z);

	//ROTATION
	Vector3 AB = Vector3Subtract(cyl.pt2, cyl.pt1);
	Quaternion qVector = QuaternionFromVector3ToVector3({0, 1, 0}, Vector3Normalize(AB));
	Quaternion qMult = QuaternionMultiply(q, qVector);

	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(qMult, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	rlScalef(cyl.radius, Vector3Length(AB), cyl.radius); // norme

	rlBegin(RL_LINES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	float deltaTheta = (endTheta - startTheta) / nSegmentsTheta;
	float theta = startTheta;
	Vector3 tmpBottomLeft = CylindricalToCartesian(Cylindrical{1, theta, 0});

	for (int i = 0; i < nSegmentsTheta; i++)
	{
		Vector3 bottomLeft = tmpBottomLeft;
		Vector3 topLeft = { bottomLeft.x, 1, bottomLeft.z};
		Vector3 bottomRight = CylindricalToCartesian({1, theta + deltaTheta, 0});
		Vector3 topRight = {bottomRight.x, 1, bottomRight.z};

		rlVertex3f(bottomLeft.x, bottomLeft.y, bottomLeft.z);
		rlVertex3f(topLeft.x, topLeft.y, topLeft.z);
		
		if(!drawCaps)
		{
			rlVertex3f(bottomLeft.x, bottomLeft.y, bottomLeft.z);
			rlVertex3f(bottomRight.x, bottomRight.y, bottomRight.z);

			rlVertex3f(topLeft.x, topLeft.y, topLeft.z);
			rlVertex3f(topRight.x, topRight.y, topRight.z);
		}

		theta += deltaTheta;
		tmpBottomLeft = bottomRight;
	}

	Vector3 bottomRight = tmpBottomLeft;
	Vector3 topRight = { bottomRight.x, 0, bottomRight.z };
	
	if(drawCaps)
	{
		Quaternion qX = QuaternionFromAxisAngle({1, 0, 0}, PI);
		Quaternion qY = QuaternionFromAxisAngle({0, 1, 0}, startTheta);
		MyDrawDiskWiresPortion(QuaternionMultiply(qX, qY), {0, 0, 0}, 1, startTheta, endTheta, nSegmentsTheta, color);
		MyDrawDiskWiresPortion(QuaternionIdentity(), {0, 1, 0}, 1, startTheta, endTheta, nSegmentsTheta, color);
	}

	rlEnd();
	rlPopMatrix();
}

void MyDrawSpherePortion(Quaternion q, Sphere sph, float startTheta, float endTheta, float startPhi, float endPhi, int nSegmentsTheta, int nSegmentsPhi, Color color) {
	if (nSegmentsTheta < 3 || nSegmentsPhi < 2) return;

	float deltaPhi = (endPhi - startPhi) / nSegmentsPhi;
	float deltaTheta = (endTheta - startTheta) / nSegmentsTheta;

	std::vector<Vector3> vertexBufferTheta(nSegmentsTheta + 1);
	std::fill(vertexBufferTheta.begin(), vertexBufferTheta.end(), Vector3{ 0,1,0 });

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
	
	rlScalef(sph.rayon, sph.rayon, sph.rayon);
	
	rlBegin(RL_TRIANGLES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	float phi = startPhi;
	for (int i = 0; i < nSegmentsPhi; i++)
	{
		float theta = startTheta;
		Vector3 tmpBottomLeft = SphericalToCartesian(Spherical{ 1, theta, phi + deltaPhi });

		for (int j = 0; j < nSegmentsTheta; j++)
		{
			Vector3 topLeft = vertexBufferTheta[j];
			Vector3 bottomLeft = tmpBottomLeft;
			Vector3 topRight = vertexBufferTheta[j + 1];
			Vector3 bottomRight = SphericalToCartesian(Spherical{ 1, theta + deltaTheta, phi + deltaPhi });

			rlVertex3f(bottomLeft.x, bottomLeft.y, bottomLeft.z);
			rlVertex3f(topRight.x, topRight.y, topRight.z);
			rlVertex3f(topLeft.x, topLeft.y, topLeft.z);

			rlVertex3f(bottomLeft.x, bottomLeft.y, bottomLeft.z);
			rlVertex3f(bottomRight.x, bottomRight.y, bottomRight.z);
			rlVertex3f(topRight.x, topRight.y, topRight.z);

			theta += deltaTheta;

			vertexBufferTheta[j] = tmpBottomLeft;
			tmpBottomLeft = bottomRight;
		}
		Vector3 topLeft = SphericalToCartesian(Spherical{ 1, endTheta, phi});;
		Vector3 topRight = vertexBufferTheta[nSegmentsTheta];
		Vector3 bottomRight = SphericalToCartesian(Spherical{ 1, endTheta, phi + deltaPhi});

		rlVertex3f(topRight.x, topRight.y, topRight.z);
		rlVertex3f(topLeft.x, topLeft.y, topLeft.z);
		rlVertex3f(bottomRight.x, bottomRight.y, bottomRight.z);
		
		vertexBufferTheta[vertexBufferTheta.size() - 1] = vertexBufferTheta[0];
		phi += deltaPhi;
	}
	
	rlEnd();
	rlPopMatrix();

}

void MyDrawSphereWiresPortion(Quaternion q, Sphere sph, float startTheta, float endTheta, float startPhi, float endPhi, int nSegmentsTheta, int nSegmentsPhi, Color color) {
	if (nSegmentsTheta < 3 || nSegmentsPhi < 2) return;

	float deltaPhi = (endPhi - startPhi) / nSegmentsPhi;
	float deltaTheta = (endTheta - startTheta) / nSegmentsTheta;
	
	std::vector<Vector3> vertexBufferTheta(nSegmentsTheta + 1);
	std::fill(vertexBufferTheta.begin(), vertexBufferTheta.end(), Vector3{ 0,1,0 });

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
	//

	rlScalef(sph.rayon, sph.rayon, sph.rayon);

	rlBegin(RL_LINES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	float phi = startPhi;
	for (int i = 0; i < nSegmentsPhi; i++)
	{
		float theta = startTheta;

		for (int j = 0; j < nSegmentsTheta; j++)
		{
			Vector3 topLeft = vertexBufferTheta[j];
			Vector3 bottomLeft = SphericalToCartesian(Spherical{ 1, theta, phi + deltaPhi });
			Vector3 topRight = vertexBufferTheta[j + 1];

			rlVertex3f(topLeft.x, topLeft.y, topLeft.z);
			rlVertex3f(topRight.x, topRight.y, topRight.z);

			rlVertex3f(topLeft.x, topLeft.y, topLeft.z);
			rlVertex3f(bottomLeft.x, bottomLeft.y, bottomLeft.z);

			theta += deltaTheta;

			vertexBufferTheta[j] = bottomLeft;
		}
		Vector3 topRight = vertexBufferTheta[nSegmentsTheta];
		Vector3 bottomRight = SphericalToCartesian(Spherical{ 1, endTheta, phi + deltaPhi});

		rlVertex3f(topRight.x, topRight.y, topRight.z);
		rlVertex3f(bottomRight.x, bottomRight.y, bottomRight.z);
		
		vertexBufferTheta[vertexBufferTheta.size() - 1] = vertexBufferTheta[0];
		phi += deltaPhi;
	}
	rlEnd();
	rlPopMatrix();
}

void MyDrawSphereEx2(Quaternion q, Sphere sph, int nSegmentsTheta, int nSegmentsPhi, Color color)
{
	MyDrawSpherePortion(q, sph, 0, 2 * PI, 0, PI, nSegmentsTheta, nSegmentsPhi, color);
}

void MyDrawSphereWiresEx2(Quaternion q, Sphere sph, int nSegmentsTheta, int nSegmentsPhi, Color color)
{
	MyDrawSphereWiresPortion(q, sph, 0, 2 * PI, 0, PI, nSegmentsTheta, nSegmentsPhi, color);
}

/* Use rlVertex3f method*/
void MyDrawQuad2(Quaternion q, Vector3 center, Vector2 size, Color color) {
	//Center - Hauteur / 2
	Vector3 point1 = Vector3SubtractValue(center, size.y / 2); // -z
	//Center + Hauteur / 2
	Vector3 point2 = Vector3AddValue(center, size.y / 2); //+z
	//Center - Largeur / 2
	Vector3 point3 = Vector3SubtractValue(center, size.x / 2); //-x
	//Center + Largeur / 2
	Vector3 point4 = Vector3AddValue(center, size.x / 2); //+x

	rlPushMatrix();
	// NOTE: Transformation is applied in inverse order (scale -> translate)
	rlTranslatef(center.x, center.y, center.z);
	
	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	//rlScalef(size.x, size.y, center.y);

	rlBegin(RL_TRIANGLES);
	//

	rlColor4ub(color.r, color.g, color.b, color.a);

	//Left
	rlVertex3f(point4.x, center.y, point1.z);
	rlVertex3f(point3.x, center.y, point1.z);
	rlVertex3f(point4.x, center.y, point2.z);

	//Right
	rlVertex3f(point3.x, center.y, point1.z);
	rlVertex3f(point3.x, center.y, point2.z);
	rlVertex3f(point4.x, center.y, point2.z);
	rlEnd();
	rlPopMatrix();
}

/* Use rlVertex3f method*/
void MyDrawQuadWire2(Quaternion q, Vector3 center, Vector2 size, Color color) {
	//Center - Hauteur / 2
	Vector3 point1 = Vector3SubtractValue(center, size.y / 2); // -z
	//Center + Hauteur / 2
	Vector3 point2 = Vector3AddValue(center, size.y / 2); //+z
	//Center - Largeur / 2
	Vector3 point3 = Vector3SubtractValue(center, size.x / 2); //-x
	//Center + Largeur / 2
	Vector3 point4 = Vector3AddValue(center, size.x / 2); //+x

	rlPushMatrix();
	// NOTE: Transformation is applied in inverse order (scale -> translate)
	rlTranslatef(center.x, center.y, center.z);
	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

	rlBegin(RL_LINES);


	rlColor4ub(color.r, color.g, color.b, color.a);

	//Left
	rlColor4ub(color.r, color.g, color.b, color.a);
	rlVertex3f(point3.x, center.y, point1.z);
	rlVertex3f(point3.x, center.y, point2.z);
	//Right
	rlVertex3f(point4.x, center.y, point1.z);
	rlVertex3f(point4.x, center.y, point2.z);
	//Up
	rlVertex3f(point3.x, center.y, point2.z);
	rlVertex3f(point4.x, center.y, point2.z);
	//Down
	rlVertex3f(point3.x, center.y, point1.z);
	rlVertex3f(point4.x, center.y, point1.z);

	//The Intersec Line
	rlVertex3f(point3.x, center.y, point1.z);
	rlVertex3f(point4.x, center.y, point2.z);
	rlEnd();
	rlPopMatrix();
}

void MyDrawCylinder(Quaternion q, Cylinder cyl, int nSegmentsTheta, bool drawCaps, Color color)
{
	MyDrawCylinderPortion(q, cyl, 0, 2 * PI, nSegmentsTheta, drawCaps, color);
}

void MyDrawCylinderWires(Quaternion q, Cylinder cyl, int nSegmentsTheta, bool drawCaps, Color color) {
	MyDrawCylinderWiresPortion(q, cyl, 0, 2 * PI, nSegmentsTheta, drawCaps, color);
}

void MyDrawDisk(Quaternion q, Vector3 center, float radius, int nSegmentsTheta, Color color) {
	MyDrawDiskPortion(q, center, radius, 0, 2 * PI, nSegmentsTheta, color);
}

void MyDrawDiskWires(Quaternion q, Vector3 center, float radius, int nSegmentsTheta, Color color) {
	MyDrawDiskWiresPortion(q, center, radius, 0, 2 * PI, nSegmentsTheta, color);
}

void MyDrawDiskPortion(Quaternion q, Vector3 center, float radius, float startTheta, float endTheta, int nSegmentsTheta, Color color)
{
	if (nSegmentsTheta < 3) return;
	int numVertex = nSegmentsTheta * 3;
	
	if (rlCheckBufferLimit(numVertex)) rlglDraw();
	
	rlPushMatrix();
	rlTranslatef(center.x, center.y, center.z);

	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	rlScalef(radius, 0, radius);

	rlBegin(RL_TRIANGLES);
	rlColor4ub(color.r, color.g, color.b, color.a);
	
	float deltaTheta = (endTheta - startTheta) / nSegmentsTheta;
	float theta = startTheta;
	Vector3 tmpBottomLeft = CylindricalToCartesian(Cylindrical{1, theta, 0});
	
	for (int i = 0; i < nSegmentsTheta; i++)
	{
		Vector3 bottomLeft = tmpBottomLeft;
		Vector3 bottomRight = CylindricalToCartesian({1, theta + deltaTheta, 0});
		
		rlVertex3f(0, 0, 0);
		rlVertex3f(bottomLeft.x, bottomLeft.y, bottomLeft.z);
		rlVertex3f(bottomRight.x, bottomRight.y, bottomRight.z);

		theta += deltaTheta;
		tmpBottomLeft = bottomRight;
	}
	
	rlEnd();
	rlPopMatrix();
}

void MyDrawDiskWiresPortion(Quaternion q, Vector3 center, float radius, float startTheta, float endTheta, int nSegmentsTheta, Color color)
{
	if (nSegmentsTheta < 3) return;
	
	rlPushMatrix();
	rlTranslatef(center.x, center.y, center.z);

	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	rlScalef(radius, radius, radius);

	rlBegin(RL_LINES);
	rlColor4ub(color.r, color.g, color.b, color.a);
	
	float deltaTheta = (endTheta - startTheta) / nSegmentsTheta;
	float theta = startTheta;
	Vector3 tmpBottomLeft = CylindricalToCartesian(Cylindrical{1, theta, 0});
	
	for (int i = 0; i < nSegmentsTheta; i++)
	{
		Vector3 bottomLeft = tmpBottomLeft;
		Vector3 bottomRight = CylindricalToCartesian({1,  theta + deltaTheta, 0});
		
		rlVertex3f(0, 0, 0);
		rlVertex3f(bottomLeft.x, bottomLeft.y, bottomLeft.z);
		rlVertex3f(bottomRight.x, bottomRight.y, bottomRight.z);
		rlVertex3f(bottomLeft.x, bottomLeft.y, bottomLeft.z);

		theta += deltaTheta;
		tmpBottomLeft = bottomRight;
	}
	
	Vector3 bottomRight = tmpBottomLeft;
	rlVertex3f(0, 0, 0);
	rlVertex3f(bottomRight.x, bottomRight.y, bottomRight.z);
	
	rlEnd();
	rlPopMatrix();
}

void MyDrawCapsule(Quaternion q, Capsule capsule, Color color) {
	rlPushMatrix();
	rlTranslatef(capsule.referential.origin.x, capsule.referential.origin.y, capsule.referential.origin.z);

	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

	rlScalef(capsule.radius, capsule.radius, capsule.radius);

	Vector3 up = { 0, 1, 0 };
	Vector3 down = { 0, 0, 0 };

	Quaternion qUp = QuaternionFromAxisAngle({ 0, 0, 1 }, 0.5 * PI);
	Quaternion qDown = QuaternionFromAxisAngle({ 0, 0, 1 }, 1.5 * PI);
	Quaternion qIdentity = QuaternionIdentity();

	Cylinder cylinder = { up, down, 1 };
	Sphere sphereUp = { up, 1 };
	Sphere sphereDown = { down, 1 };

	MyDrawCylinder(qIdentity, cylinder, 25, true, color);
	MyDrawSpherePortion(qUp, sphereUp, 0, PI, 0, PI, 30, 30, color);
	MyDrawSpherePortion(qDown, sphereDown, 0, PI, 0, PI, 25, 25, color);;

	rlPopMatrix();
}

void MyDrawCapsuleWires(Quaternion q, Capsule capsule, Color color) {
	rlPushMatrix();
	rlTranslatef(capsule.referential.origin.x, capsule.referential.origin.y, capsule.referential.origin.z);

	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

	rlScalef(capsule.radius, capsule.radius, capsule.radius);

	Vector3 up = { 0, 1, 0 };
	Vector3 down = { 0, 0, 0 };

	Quaternion qUp = QuaternionFromAxisAngle({ 0, 0, 1 }, 0.5 * PI);
	Quaternion qDown = QuaternionFromAxisAngle({ 0, 0, 1 }, 1.5 * PI);
	Quaternion qIdentity = QuaternionIdentity();

	Cylinder cylinder = { up, down, 1 };
	Sphere sphereUp = { up, 1 };
	Sphere sphereDown = { down, 1 };
	
	MyDrawCylinderWires(qIdentity, cylinder, 25, true, color);
	MyDrawSphereWiresPortion(qUp, sphereUp, 0, PI, 0, PI, 30, 30, color);
	MyDrawSphereWiresPortion(qDown, sphereDown, 0, PI, 0, PI, 25, 25, color);

	rlPopMatrix();
}

void MyDrawRoundBox(Quaternion q, Vector3 center, Vector3 size, Color color) {
	rlPushMatrix();
	rlTranslatef(center.x, center.y, center.z);

	Quaternion qUp = QuaternionIdentity();
	Quaternion qLeft = QuaternionFromAxisAngle({ 1, 0, 0 }, PI * 0.5f);
	Quaternion qRight = QuaternionFromAxisAngle({ 1, 0, 0 }, PI * -0.5f);
	Quaternion qFront = QuaternionFromAxisAngle({ 0, 0, 1 }, PI * 0.5f);
	Quaternion qBack = QuaternionFromAxisAngle({ 0, 0, 1 }, PI * -0.5f);

	Quaternion qX = QuaternionFromAxisAngle({ 1, 0, 0 }, PI);
	Quaternion qY = QuaternionFromAxisAngle({ 0, 1, 0 }, 0);
	Quaternion qDown = QuaternionMultiply(qX, qY);

	Referential referential = Referential({ 0,1,0 });

	/*Quad quad = { referential, size };
	Quad quadDown = { Referential({0, 0, 0}), size };*/

	/*MyDrawQuad2(qUp, quad.referential.origin, {quad.extension.x, quad.extension.z}, color);
	MyDrawQuad2(qLeft, quad.referential.origin, { quad.extension.x, quad.extension.z }, color);
	MyDrawQuad2(qRight, quad.referential.origin, {quad.extension.x, quad.extension.z}, color);
	MyDrawQuad2(qFront, quad.referential.origin, { quad.extension.x, quad.extension.z }, color);
	MyDrawQuad2(qBack, quad.referential.origin, { quad.extension.x, quad.extension.z }, color);
	MyDrawQuad2(qDown, quadDown.referential.origin, {quadDown.extension.x, quadDown.extension.z}, color);*/

	Capsule capsule = { referential, 1 };
	Capsule capsuleDown = { Referential({0, 0, 0}), 1 };

	//MyDrawCapsule(qUp, capsule, color);
	MyDrawCapsule(qLeft, capsule, color);
	//MyDrawCapsule(qRight, capsule, color);
	//MyDrawCapsule(qFront, capsule, color);
	//MyDrawCapsule(qBack, capsule, color);
	//MyDrawCapsule(qDown, capsuleDown, color);

	rlPopMatrix();
}

void MyDrawRoundBoxWires(Quaternion q, Vector3 center, Vector3 size, Color color) {
	rlPushMatrix();
	rlTranslatef(center.x, center.y, center.z);

	Quaternion qUp = QuaternionIdentity();
	Quaternion qLeft = QuaternionFromAxisAngle({ 1, 0, 0 }, PI * 0.5f);
	Quaternion qRight = QuaternionFromAxisAngle({ 1, 0, 0 }, PI * -0.5f);
	Quaternion qFront = QuaternionFromAxisAngle({ 0, 0, 1 }, PI * 0.5f);
	Quaternion qBack = QuaternionFromAxisAngle({ 0, 0, 1 }, PI * -0.5f);

	Quaternion qX = QuaternionFromAxisAngle({ 1, 0, 0 }, PI);
	Quaternion qY = QuaternionFromAxisAngle({ 0, 1, 0 }, 0);
	Quaternion qDown = QuaternionMultiply(qX, qY);

	Referential referential = Referential({ 0,1,0 });

	Quad quad = { referential, size };
	Quad quadDown = { Referential({0, 0, 0}), size };

	/*MyDrawQuadWire2(qUp, quad.referential.origin, {quad.extension.x, quad.extension.z}, color);
	MyDrawQuadWire2(qLeft, quad.referential.origin, { quad.extension.x, quad.extension.z }, color);
	MyDrawQuadWire2(qRight, quad.referential.origin, {quad.extension.x, quad.extension.z}, color);
	MyDrawQuadWire2(qFront, quad.referential.origin, { quad.extension.x, quad.extension.z }, color);
	MyDrawQuadWire2(qBack, quad.referential.origin, { quad.extension.x, quad.extension.z }, color);
	MyDrawQuadWire2(qDown, quadDown.referential.origin, {quadDown.extension.x, quadDown.extension.z}, color);*/

	Referential referentialUp = Referential({ 0,1,0 });
	Referential referentialDown = Referential({ 0,0,0 });
	Referential referentialLeft = Referential({ 1,1,0 });
	Referential referentialRight = Referential({ -1,1,0 });

	Capsule capsuleUp = { referential, 1 };
	Capsule capsuleLeft = { referential, 1 };
	Capsule capsuleRight = { referential, 1 };
	Capsule capsuleDown = { referential, 1 };

	/*MyDrawCapsuleWires(qUp, capsule, color);
	MyDrawCapsuleWires(qLeft, capsule, color);
	MyDrawCapsuleWires(qRight, capsule, color);
	MyDrawCapsuleWires(qFront, capsule, color);
	MyDrawCapsuleWires(qBack, capsule, color);
	MyDrawCapsuleWires(qDown, capsuleDown, color);*/

	rlPopMatrix();
}
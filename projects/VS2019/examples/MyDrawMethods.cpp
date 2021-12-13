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
	Vector3 PQ = Vector3Subtract(cyl.pt2, cyl.pt1);
	//Quaternion qVector = QuaternionFromVector3ToVector3({0, 1, 0}, Vector3Normalize(AB));
	//Quaternion qMult = QuaternionMultiply(q, qVector);

	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	rlScalef(cyl.radius, Vector3Length(PQ), cyl.radius); // norme

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
		MyDrawDiskPortion(QuaternionMultiply(qX, qY), { Referential({0, 0, 0}), 1 }, startTheta, endTheta, nSegmentsTheta, color);
		MyDrawDiskPortion(QuaternionIdentity(), { Referential({0, 1, 0}), 1 }, startTheta, endTheta, nSegmentsTheta, color);
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
	//Quaternion qVector = QuaternionFromVector3ToVector3({0, 1, 0}, Vector3Normalize(AB));
	//Quaternion qMult = QuaternionMultiply(q, qVector);

	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(q, &vect, &angle);
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
		MyDrawDiskWiresPortion(QuaternionMultiply(qX, qY), { Referential({0, 0, 0}), 1 }, startTheta, endTheta, nSegmentsTheta, color);
		MyDrawDiskWiresPortion(QuaternionIdentity(), { Referential({0, 1, 0}), 1 }, startTheta, endTheta, nSegmentsTheta, color);
	}

	rlEnd();
	rlPopMatrix();
}

void MyDrawSpherePortion(Quaternion q, Sphere sph, float startTheta, float endTheta, float startPhi, float endPhi, int nSegmentsTheta, int nSegmentsPhi, Color color) {
	if (nSegmentsTheta < 3 || nSegmentsPhi < 2) return;

	float deltaPhi = (endPhi - startPhi) / nSegmentsPhi;
	float deltaTheta = (endTheta - startTheta) / nSegmentsTheta;

	std::vector<Vector3> vertexBufferTheta(nSegmentsTheta + 1);
	for (int i = 0; i < nSegmentsTheta + 1; i++)
		vertexBufferTheta[i] = SphericalToCartesian({ 1,Lerp(startTheta,endTheta,i / (float)nSegmentsTheta),startPhi });

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

		vertexBufferTheta[vertexBufferTheta.size() - 1] = SphericalToCartesian(Spherical{ 1, theta + deltaTheta, phi + deltaPhi });

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
	for (size_t i = 0; i < nSegmentsTheta + 1; i++)
		vertexBufferTheta[i] = SphericalToCartesian({ 1,Lerp(startTheta,endTheta,i / (float)nSegmentsTheta),startPhi });
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
			Vector3 topLeft = SphericalToCartesian(Spherical{ 1, theta , phi });// vertexBufferTheta[j];
			Vector3 bottomLeft = SphericalToCartesian(Spherical{ 1, theta , phi + deltaPhi });//tmpBottomLeft;
			Vector3 topRight = SphericalToCartesian(Spherical{ 1, theta + deltaTheta, phi });//vertexBufferTheta[j + 1];
			/*Vector3 topLeft = vertexBufferTheta[j];
			Vector3 bottomLeft = SphericalToCartesian(Spherical{ 1, theta, phi + deltaPhi });
			Vector3 topRight = vertexBufferTheta[j + 1];*/

			rlVertex3f(topLeft.x, topLeft.y, topLeft.z);
			rlVertex3f(topRight.x, topRight.y, topRight.z);

			rlVertex3f(topLeft.x, topLeft.y, topLeft.z);
			rlVertex3f(bottomLeft.x, bottomLeft.y, bottomLeft.z);

			theta += deltaTheta;

			vertexBufferTheta[j] = bottomLeft;
		}
		Vector3 topRight = SphericalToCartesian(Spherical{ 1, theta + deltaTheta, phi });//vertexBufferTheta[j + 1];
		//Vector3 topRight = vertexBufferTheta[nSegmentsTheta];
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
void MyDrawQuad2(Quaternion q, Quad quad, Color color) {
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
	QuaternionToAxisAngle(q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	//rlScalef(size.x, size.y, center.y);

	rlBegin(RL_TRIANGLES);
	//

	rlColor4ub(color.r, color.g, color.b, color.a);

	//Left
	rlVertex3f(x, 0, zNeg);
	rlVertex3f(xNeg, 0, zNeg);
	rlVertex3f(x, 0, z);

	//Right
	rlVertex3f(xNeg, 0, zNeg);
	rlVertex3f(xNeg, 0, z);
	rlVertex3f(x, 0, z);
	rlEnd();
	rlPopMatrix();
}

/* Use rlVertex3f method*/
void MyDrawQuadWire2(Quaternion q, Quad quad, Color color) {
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
	QuaternionToAxisAngle(q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

	rlBegin(RL_LINES);


	rlColor4ub(color.r, color.g, color.b, color.a);

	//Left
	rlVertex3f(xNeg, 0, zNeg);
	rlVertex3f(xNeg, 0, z);
	//Right
	rlVertex3f(x, 0, zNeg);
	rlVertex3f(x, 0, z);
	//Up
	rlVertex3f(xNeg, 0, z);
	rlVertex3f(x, 0, z);
	//Down
	rlVertex3f(xNeg, 0, zNeg);
	rlVertex3f(x, 0, zNeg);

	//The Intersec Line
	rlVertex3f(xNeg, 0, zNeg);
	rlVertex3f(x, 0, z);

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

void MyDrawDisk(Quaternion q, Disk disk, int nSegmentsTheta, Color color) {
	MyDrawDiskPortion(q, disk, 0, 2 * PI, nSegmentsTheta, color);
}

void MyDrawDiskWires(Quaternion q, Disk disk, int nSegmentsTheta, Color color) {
	MyDrawDiskWiresPortion(q, disk, 0, 2 * PI, nSegmentsTheta, color);
}

void MyDrawDiskPortion(Quaternion q, Disk disk, float startTheta, float endTheta, int nSegmentsTheta, Color color)
{
	if (nSegmentsTheta < 3) return;
	int numVertex = nSegmentsTheta * 3;
	
	if (rlCheckBufferLimit(numVertex)) rlglDraw();
	
	rlPushMatrix();
	rlTranslatef(disk.referential.origin.x, disk.referential.origin.y, disk.referential.origin.z);

	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	rlScalef(disk.radius, 0, disk.radius);

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

void MyDrawDiskWiresPortion(Quaternion q, Disk disk, float startTheta, float endTheta, int nSegmentsTheta, Color color)
{
	if (nSegmentsTheta < 3) return;
	
	rlPushMatrix();
	rlTranslatef(disk.referential.origin.x, disk.referential.origin.y, disk.referential.origin.z);

	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);
	rlScalef(disk.radius, 0, disk.radius);

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

	//rlScalef(capsule.radius, 1, capsule.radius);

	Vector3 up = { 0, capsule.height, 0 };
	Vector3 down = { 0, 0, 0 };

	Quaternion qUp = QuaternionFromAxisAngle({ 0, 0, 1 }, 0.5 * PI);
	Quaternion qDown = QuaternionFromAxisAngle({ 0, 0, 1 }, -0.5 * PI);
	Quaternion qIdentity = QuaternionIdentity();

	Referential ref = Referential(down);
	ref.RotateByQuaternion(qIdentity);
	Cylinder cylinder = Cylinder(ref, capsule.radius, capsule.height);

	Sphere sphereUp = { up, capsule.radius };
	Sphere sphereDown = { down, capsule.radius };

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

	//rlScalef(capsule.radius, 1, capsule.radius);

	Vector3 up = {0, capsule.height, 0};
	Vector3 down = { 0, 0, 0 };

	Quaternion qUp = QuaternionFromAxisAngle({ 0, 0, 1 }, 0.5 * PI);
	Quaternion qDown = QuaternionFromAxisAngle({ 0, 0, 1 }, 1.5 * PI);
	Quaternion qIdentity = QuaternionIdentity();

	Referential ref = Referential(down);
	ref.RotateByQuaternion(qIdentity);
	Cylinder cylinder = Cylinder(ref, capsule.radius, capsule.height);

	Sphere sphereUp = { up, capsule.radius };
	Sphere sphereDown = { down, capsule.radius };
	
	MyDrawCylinderWires(qIdentity, cylinder, 25, true, color);
	MyDrawSphereWiresPortion(qUp, sphereUp, 0, PI, 0, PI, 30, 30, color);
	MyDrawSphereWiresPortion(qDown, sphereDown, 0, PI, 0, PI, 25, 25, color);

	rlPopMatrix();
}

void MyDrawRoundBox(Quaternion q, RoundedBox roundedBox, Color color) {
	rlPushMatrix();
	rlTranslatef(roundedBox.ref.origin.x, roundedBox.ref.origin.y, roundedBox.ref.origin.z);

	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

	rlScalef(roundedBox.extension.x, roundedBox.extension.y, roundedBox.extension.z);

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
	referentialQuadRight.RotateByQuaternion(qRight);
	referentialQuadDown.RotateByQuaternion(qDown);

	Quad quadUp = { referentialQuadUp, {1, 1, 1} };
	Quad quadFront = { referentialQuadFront, {1, 1, 1} };
	Quad quadBack = { referentialQuadBack, {1, 1, 1} };
	Quad quadLeft = { referentialQuadLeft, {1, 1, 1} };
	Quad quadRight = { referentialQuadRight, {1, 1, 1} };
	Quad quadDown = { referentialQuadDown, {1, 1, 1} };

	MyDrawQuad2(qUp, quadUp, color);
	MyDrawQuad2(qLeft, quadLeft, color);
	MyDrawQuad2(qRight, quadRight, color);
	MyDrawQuad2(qFront, quadFront, color);
	MyDrawQuad2(qBack, quadBack, color);
	MyDrawQuad2(qDown, quadDown, color);

	Referential referentiaCapsUp = Referential({ 0, 0.5f, -0.5f });
	Referential referentialCapsDown = Referential({ 0, -0.5f, -0.5f });
	Referential referentialCapsRight = Referential({ 0, -0.5f, -0.5f });
	Referential referentialCapsLeft = Referential({ 0, -0.5f, 0.5f });

	Capsule capsuleUp = { referentiaCapsUp, roundedBox.radius, 1 };
	Capsule capsuleDown = { referentialCapsDown, roundedBox.radius, 1 };
	Capsule capsuleRight = { referentialCapsRight, roundedBox.radius, 1 };
	Capsule capsuleLeft = { referentialCapsLeft, roundedBox.radius, 1 };

	MyDrawCapsule(qLeft, capsuleUp, color);
	MyDrawCapsule(qUp, capsuleRight, color);
	MyDrawCapsule(qUp, capsuleLeft, color);
	MyDrawCapsule(qLeft, capsuleDown, color);

	Referential referentiaCapsUpBack = Referential({ -1, 0.5f, -0.5f });
	Referential referentialCapsDownBack = Referential({ -1, -0.5f, -0.5f });
	Referential referentialCapsLeftBack = Referential({ -1, -0.5f, -0.5f });
	Referential referentialCapsRightBack = Referential({ -1, -0.5f, 0.5f });

	Capsule capsuleUpBack = { referentiaCapsUpBack, roundedBox.radius, 1 };
	Capsule capsuleDownBack = { referentialCapsDownBack, roundedBox.radius, 1 };
	Capsule capsuleLeftBack = { referentialCapsLeftBack, roundedBox.radius, 1 };
	Capsule capsuleRightBack = { referentialCapsRightBack, roundedBox.radius, 1 };

	MyDrawCapsule(qLeft, capsuleUpBack, color);
	MyDrawCapsule(qLeft, capsuleDownBack, color);
	MyDrawCapsule(qUp, capsuleLeftBack, color);
	MyDrawCapsule(qUp, capsuleRightBack, color);

	Referential referentiaCapsUpLeft = Referential({ -1, 0.5f, -0.5f });
	Referential referentialCapsDownLeft = Referential({ -1, -0.5f, -0.5f });

	Capsule capsuleUpRight = { referentiaCapsUpLeft, roundedBox.radius, 1 };
	Capsule capsuleDownRight = { referentialCapsDownLeft, roundedBox.radius, 1 };

	MyDrawCapsule(qFront, capsuleUpRight, color);
	MyDrawCapsule(qFront, capsuleDownRight, color);

	Referential referentiaCapsUpRight = Referential({ -1, 0.5f, 0.5f });
	Referential referentialCapsDownRight = Referential({ -1, -0.5f, 0.5f });

	Capsule capsuleUpLeft = { referentiaCapsUpRight, roundedBox.radius, 1 };
	Capsule capsuleDownLeft = { referentialCapsDownRight, roundedBox.radius, 1 };

	MyDrawCapsule(qFront, capsuleUpLeft, color);
	MyDrawCapsule(qFront, capsuleDownLeft, color);

	rlPopMatrix();
}

/*void MyDrawRoundBox(Quaternion q, RoundedBox roundexBox, Color color) {
	rlPushMatrix();
	//rlTranslatef(roundexBox.ref.origin.x, roundexBox.ref.origin.y, roundexBox.ref.origin.z);

	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(q, &vect, &angle);
	//rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

	//rlScalef(roundexBox.extension.x, roundexBox.extension.y, roundexBox.extension.z);

	Quaternion qUp = QuaternionIdentity();
	Quaternion qLeft = QuaternionFromAxisAngle(roundexBox.ref.i, PI * 0.5f);
	Quaternion qRight = QuaternionFromAxisAngle(roundexBox.ref.i, PI * -0.5f);
	Quaternion qFront = QuaternionFromAxisAngle(roundexBox.ref.k, PI * -0.5f);
	Quaternion qBack = QuaternionFromAxisAngle(roundexBox.ref.k, PI * 0.5f);
	//Quaternion qDown = QuaternionFromAxisAngle({ 1, 0, 0 }, PI * 1);

	Quaternion qX = QuaternionFromAxisAngle({ 1, 0, 0 }, PI);
	Quaternion qY = QuaternionFromAxisAngle({ 0, 1, 0 }, 0);
	Quaternion qDown = QuaternionMultiply(qX, qY);

	Referential referentialQuadUp = Referential(Vector3Add(roundexBox.ref.origin, Vector3Scale(roundexBox.ref.j, roundexBox.extension.y)));
	Referential referentialQuadFront = Referential(Vector3Add(roundexBox.ref.origin, Vector3Scale(roundexBox.ref.i, roundexBox.extension.x)));
	Referential referentialQuadBack = Referential(Vector3Add(roundexBox.ref.origin, Vector3Scale(Vector3Negate(roundexBox.ref.i), roundexBox.extension.x)));
	Referential referentialQuadLeft = Referential(Vector3Add(roundexBox.ref.origin, Vector3Scale(roundexBox.ref.k, roundexBox.extension.z)));
	Referential referentialQuadRight = Referential(Vector3Add(roundexBox.ref.origin, Vector3Scale(Vector3Negate(roundexBox.ref.k), roundexBox.extension.z)));
	Referential referentialQuadDown = Referential(Vector3Add(roundexBox.ref.origin, Vector3Scale(Vector3Negate(roundexBox.ref.j), roundexBox.extension.y)));

	referentialQuadUp.RotateByQuaternion(qUp);
	referentialQuadFront.RotateByQuaternion(qFront);
	referentialQuadBack.RotateByQuaternion(qBack);
	referentialQuadLeft.RotateByQuaternion(qLeft);
	referentialQuadRight.RotateByQuaternion(qRight);
	referentialQuadDown.RotateByQuaternion(qDown);

	Quad quadUp = { referentialQuadUp, {roundexBox.extension.x, 0, roundexBox.extension.z} };
	Quad quadFront = { referentialQuadFront, {roundexBox.extension.y, 0, roundexBox.extension.z} };
	Quad quadBack = { referentialQuadBack, {roundexBox.extension.y, 0, roundexBox.extension.z} };
	Quad quadLeft = { referentialQuadLeft, {roundexBox.extension.x, 0, roundexBox.extension.y} };
	Quad quadRight = { referentialQuadRight, {roundexBox.extension.x, 0, roundexBox.extension.y} };
	Quad quadDown = { referentialQuadDown, {roundexBox.extension.x, 0, roundexBox.extension.z} };

	MyDrawQuad2(qUp, quadUp, RED);
	MyDrawQuad2(qLeft, quadLeft, YELLOW);
	MyDrawQuad2(qRight, quadRight, GREEN);
	MyDrawQuad2(qFront, quadFront, BLUE);
	MyDrawQuad2(qBack, quadBack, PURPLE);
	MyDrawQuad2(qDown, quadDown, DARKGRAY);


	Vector3 iPrime = Vector3Add(roundexBox.ref.origin, Vector3Scale(roundexBox.ref.i, roundexBox.extension.x));
	Vector3 jPrime = Vector3Add(roundexBox.ref.origin, Vector3Scale(roundexBox.ref.j, roundexBox.extension.y));
	Vector3 zPrime = Vector3Add(roundexBox.ref.origin, Vector3Scale(roundexBox.ref.k, roundexBox.extension.z));
	Referential referentialCapsUp = Referential(Vector3Add(iPrime, jPrime));
	Referential referentialCapsDown = Referential(Vector3Add(iPrime, Vector3Negate(jPrime)));
	Referential referentialCapsRight = Referential(Vector3Add(iPrime, zPrime));
	Referential referentialCapsLeft = Referential(Vector3Add(iPrime, Vector3Negate(zPrime)));

	Capsule capsuleUp = { referentialCapsUp, roundexBox.radius };
	Capsule capsuleDown = { referentialCapsDown, roundexBox.radius };
	Capsule capsuleRight = { referentialCapsRight, roundexBox.radius };
	Capsule capsuleLeft = { referentialCapsLeft, roundexBox.radius };

	MyDrawCapsule(qLeft, capsuleUp, color);
	MyDrawCapsule(qUp, capsuleRight, color);
	MyDrawCapsule(qUp, capsuleLeft, color);
	MyDrawCapsule(qLeft, capsuleDown, color);
	/
	Referential referentiaCapsUpBack = Referential({-1, 0.5f, -0.5f});
	Referential referentialCapsDownBack = Referential({ -1, -0.5f, -0.5f });
	Referential referentialCapsLeftBack = Referential({ -1, -0.5f, -0.5f });
	Referential referentialCapsRightBack = Referential({ -1, -0.5f, 0.5f });

	Capsule capsuleUpBack = { referentiaCapsUpBack, roundexBox.radius };
	Capsule capsuleDownBack = { referentialCapsDownBack, roundexBox.radius };
	Capsule capsuleLeftBack = { referentialCapsLeftBack, roundexBox.radius };
	Capsule capsuleRightBack = { referentialCapsRightBack, roundexBox.radius };

	MyDrawCapsule(qLeft, capsuleUpBack, color);
	MyDrawCapsule(qLeft, capsuleDownBack, color);
	MyDrawCapsule(qUp, capsuleLeftBack, color);
	MyDrawCapsule(qUp, capsuleRightBack, color);

	Referential referentiaCapsUpLeft = Referential({ -1, 0.5f, -0.5f });
	Referential referentialCapsDownLeft = Referential({ -1, -0.5f, -0.5f });

	Capsule capsuleUpRight = { referentiaCapsUpLeft, roundexBox.radius };
	Capsule capsuleDownRight = { referentialCapsDownLeft, roundexBox.radius };

	MyDrawCapsule(qFront, capsuleUpRight, color);
	MyDrawCapsule(qFront, capsuleDownRight, color);

	Referential referentiaCapsUpRight = Referential({ -1, 0.5f, 0.5f });
	Referential referentialCapsDownRight = Referential({ -1, -0.5f, 0.5f });

	Capsule capsuleUpLeft = { referentiaCapsUpRight, roundexBox.radius };
	Capsule capsuleDownLeft = { referentialCapsDownRight, roundexBox.radius };

	MyDrawCapsule(qFront, capsuleUpLeft, color);
	MyDrawCapsule(qFront, capsuleDownLeft, color);

	rlPopMatrix();
}*/

void MyDrawRoundBoxWires(Quaternion q, RoundedBox roundexBox, Color color) {
	rlPushMatrix();
	rlTranslatef(roundexBox.ref.origin.x, roundexBox.ref.origin.y, roundexBox.ref.origin.z);

	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

	rlScalef(roundexBox.extension.x, roundexBox.extension.y, roundexBox.extension.z);

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

	Quad quadUp = { referentialQuadUp, {1, 1, 1} };
	Quad quadFront = { referentialQuadFront, {1, 1, 1} };
	Quad quadBack = { referentialQuadBack, {1, 1, 1} };
	Quad quadLeft = { referentialQuadLeft, {1, 1, 1} };
	Quad quadRight = { referentialQuadRight, {1, 1, 1} };
	Quad quadDown = { referentialQuadDown, {1, 1, 1} };

	MyDrawQuadWire2(qUp, quadUp, color);
	MyDrawQuadWire2(qLeft, quadLeft, color);
	MyDrawQuadWire2(qRight, quadRight, color);
	MyDrawQuadWire2(qFront, quadFront, color);
	MyDrawQuadWire2(qBack, quadBack, color);
	MyDrawQuadWire2(qDown, quadDown, color);

	Referential referentiaCapsUp = Referential({ 0, 0.5f, -0.5f });
	Referential referentialCapsDown = Referential({ 0, -0.5f, -0.5f });
	Referential referentialCapsRight = Referential({ 0, -0.5f, -0.5f });
	Referential referentialCapsLeft = Referential({ 0, -0.5f, 0.5f });

	referentiaCapsUp.RotateByQuaternion(qLeft);
	referentialCapsRight.RotateByQuaternion(qUp);
	referentialCapsLeft.RotateByQuaternion(qUp);
	referentialCapsDown.RotateByQuaternion(qLeft);

	Capsule capsuleUp = { referentiaCapsUp, roundexBox.radius, 1 };
	Capsule capsuleDown = { referentialCapsDown, roundexBox.radius, 1 };
	Capsule capsuleRight = { referentialCapsRight, roundexBox.radius, 1 };
	Capsule capsuleLeft = { referentialCapsLeft, roundexBox.radius, 1 };

	MyDrawCapsuleWires(qLeft, capsuleUp, color);
	MyDrawCapsuleWires(qUp, capsuleRight, color);
	MyDrawCapsuleWires(qUp, capsuleLeft, color);
	MyDrawCapsuleWires(qLeft, capsuleDown, color);

	Referential referentiaCapsUpBack = Referential({ -1, 0.5f, -0.5f });
	Referential referentialCapsDownBack = Referential({ -1, -0.5f, -0.5f });
	Referential referentialCapsLeftBack = Referential({ -1, -0.5f, -0.5f });
	Referential referentialCapsRightBack = Referential({ -1, -0.5f, 0.5f });

	referentiaCapsUpBack.RotateByQuaternion(qLeft);
	referentialCapsDownBack.RotateByQuaternion(qLeft);
	referentialCapsLeftBack.RotateByQuaternion(qUp);
	referentialCapsRightBack.RotateByQuaternion(qUp);

	Capsule capsuleUpBack = { referentiaCapsUpBack, roundexBox.radius, 1 };
	Capsule capsuleDownBack = { referentialCapsDownBack, roundexBox.radius, 1 };
	Capsule capsuleLeftBack = { referentialCapsLeftBack, roundexBox.radius, 1 };
	Capsule capsuleRightBack = { referentialCapsRightBack, roundexBox.radius, 1 };

	MyDrawCapsuleWires(qLeft, capsuleUpBack, color);
	MyDrawCapsuleWires(qLeft, capsuleDownBack, color);
	MyDrawCapsuleWires(qUp, capsuleLeftBack, color);
	MyDrawCapsuleWires(qUp, capsuleRightBack, color);

	Referential referentiaCapsUpLeft = Referential({ -1, 0.5f, -0.5f });
	Referential referentialCapsDownLeft = Referential({ -1, -0.5f, -0.5f });

	referentiaCapsUpLeft.RotateByQuaternion(qFront);
	referentialCapsDownLeft.RotateByQuaternion(qFront);

	Capsule capsuleUpRight = { referentiaCapsUpLeft, roundexBox.radius, 1 };
	Capsule capsuleDownRight = { referentialCapsDownLeft, roundexBox.radius, 1 };

	MyDrawCapsuleWires(qFront, capsuleUpRight, color);
	MyDrawCapsuleWires(qFront, capsuleDownRight, color);

	Referential referentiaCapsUpRight = Referential({ -1, 0.5f, 0.5f });
	Referential referentialCapsDownRight = Referential({ -1, -0.5f, 0.5f });

	referentiaCapsUpRight.RotateByQuaternion(qFront);
	referentialCapsDownRight.RotateByQuaternion(qFront);

	Capsule capsuleUpLeft = { referentiaCapsUpRight, roundexBox.radius, 1 };
	Capsule capsuleDownLeft = { referentialCapsDownRight, roundexBox.radius, 1 };

	MyDrawCapsuleWires(qFront, capsuleUpLeft, color);
	MyDrawCapsuleWires(qFront, capsuleDownLeft, color);

	rlPopMatrix();
}

void MyDrawBox(Quaternion q, Box box, Color color) {
	rlPushMatrix();
	rlTranslatef(box.ref.origin.x, box.ref.origin.y, box.ref.origin.z);

	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

	rlScalef(box.extension.x, box.extension.y, box.extension.z);

	Quaternion qUp = QuaternionIdentity();
	Quaternion qLeft = QuaternionFromAxisAngle({ 1, 0, 0 }, PI * 0.5f);
	Quaternion qRight = QuaternionFromAxisAngle({ 1, 0, 0 }, PI * -0.5f);
	Quaternion qFront = QuaternionFromAxisAngle({ 0, 0, 1 }, PI * -0.5f);
	Quaternion qBack = QuaternionFromAxisAngle({ 0, 0, 1 }, PI * 0.5f);
	//Quaternion qDown = QuaternionFromAxisAngle({ 1, 0, 0 }, PI * 1);

	Quaternion qX = QuaternionFromAxisAngle({ 1, 0, 0 }, PI);
	Quaternion qY = QuaternionFromAxisAngle({ 0, 1, 0 }, 0);
	Quaternion qDown = QuaternionMultiply(qX, qY);

	Referential referentialQuadUp = Referential({ -0.5f, 0.5f, 0 });
	Referential referentialQuadFront = Referential({ 0, 0, 0 });
	Referential referentialQuadBack = Referential({ -1, 0, 0 });
	Referential referentialQuadLeft = Referential({ -0.5f, 0, 0.5f });
	Referential referentialQuadRight = Referential({ -0.5f, 0, -0.5f });
	Referential referentialQuadDown = Referential({ -0.5f, -0.5f, 0 });

	referentialQuadUp.RotateByQuaternion(qUp);
	referentialQuadFront.RotateByQuaternion(qFront);
	referentialQuadBack.RotateByQuaternion(qBack);
	referentialQuadLeft.RotateByQuaternion(qLeft);
	referentialQuadDown.RotateByQuaternion(qDown);

	Quad quadUp = { referentialQuadUp, {1, 1, 1} };
	Quad quadFront = { referentialQuadFront, {1, 1, 1} };
	Quad quadBack = { referentialQuadBack, {1, 1, 1} };
	Quad quadLeft = { referentialQuadLeft, {1, 1, 1} };
	Quad quadRight = { referentialQuadRight, {1, 1, 1} };
	Quad quadDown = { referentialQuadDown, {1, 1, 1} };

	MyDrawQuad2(qUp, quadUp, color);
	MyDrawQuad2(qLeft, quadLeft, color);
	MyDrawQuad2(qRight, quadRight, color);
	MyDrawQuad2(qFront, quadFront, color);
	MyDrawQuad2(qBack, quadBack, color);
	MyDrawQuad2(qDown, quadDown, color);

	rlPopMatrix();
}

void MyDrawBoxWires(Quaternion q, Box box, Color color) {
	rlPushMatrix();
	rlTranslatef(box.ref.origin.x, box.ref.origin.y, box.ref.origin.z);

	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

	rlScalef(box.extension.x, box.extension.y, box.extension.z);

	Quaternion qUp = QuaternionIdentity();
	Quaternion qLeft = QuaternionFromAxisAngle({ 1, 0, 0 }, PI * 0.5f);
	Quaternion qRight = QuaternionFromAxisAngle({ 1, 0, 0 }, PI * -0.5f);
	Quaternion qFront = QuaternionFromAxisAngle({ 0, 0, 1 }, PI * -0.5f);
	Quaternion qBack = QuaternionFromAxisAngle({ 0, 0, 1 }, PI * 0.5f);
	//Quaternion qDown = QuaternionFromAxisAngle({ 1, 0, 0 }, PI * 1);

	Quaternion qX = QuaternionFromAxisAngle({ 1, 0, 0 }, PI);
	Quaternion qY = QuaternionFromAxisAngle({ 0, 1, 0 }, 0);
	Quaternion qDown = QuaternionMultiply(qX, qY);

	Referential referentialQuadUp = Referential({ -0.5f, 0.5f, 0 });
	Referential referentialQuadFront = Referential({ 0, 0, 0 });
	Referential referentialQuadBack = Referential({ -1, 0, 0 });
	Referential referentialQuadLeft = Referential({ -0.5f, 0, 0.5f });
	Referential referentialQuadRight = Referential({ -0.5f, 0, -0.5f });
	Referential referentialQuadDown = Referential({ -0.5f, -0.5f, 0 });

	referentialQuadUp.RotateByQuaternion(qUp);
	referentialQuadFront.RotateByQuaternion(qFront);
	referentialQuadBack.RotateByQuaternion(qBack);
	referentialQuadLeft.RotateByQuaternion(qLeft);
	referentialQuadDown.RotateByQuaternion(qDown);

	Quad quadUp = { referentialQuadUp, {1, 1, 1} };
	Quad quadFront = { referentialQuadFront, {1, 1, 1} };
	Quad quadBack = { referentialQuadBack, {1, 1, 1} };
	Quad quadLeft = { referentialQuadLeft, {1, 1, 1} };
	Quad quadRight = { referentialQuadRight, {1, 1, 1} };
	Quad quadDown = { referentialQuadDown, {1, 1, 1} };

	MyDrawQuadWire2(qUp, quadUp, color);
	MyDrawQuadWire2(qLeft, quadLeft, color);
	MyDrawQuadWire2(qRight, quadRight, color);
	MyDrawQuadWire2(qFront, quadFront, color);
	MyDrawQuadWire2(qBack, quadBack, color);
	MyDrawQuadWire2(qDown, quadDown, color);

	rlPopMatrix();
}
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

void MyDrawCylinderPortion(Cylinder cyl, float startTheta, float endTheta, int nSegmentsTheta, bool drawCaps, Color color) {
	if (nSegmentsTheta < 3) return;

	int numVertex = nSegmentsTheta * 6;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	rlPushMatrix();
	rlTranslatef(cyl.pt1.x, cyl.pt1.y, cyl.pt1.z);

	//ROTATION
	Vector3 PQ = Vector3Subtract(cyl.pt2, cyl.pt1);

	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(cyl.ref.q, &vect, &angle);
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
		MyDrawDiskPortion({ Referential({0, 0, 0}, QuaternionMultiply(qX, qY)), 1 }, startTheta, endTheta, nSegmentsTheta, color);
		MyDrawDiskPortion({ Referential({0, 1, 0}, QuaternionIdentity()), 1 }, startTheta, endTheta, nSegmentsTheta, color);
	}
	
	rlEnd();
	rlPopMatrix();
}

void MyDrawCylinderWiresPortion(Cylinder cyl, float startTheta, float endTheta, int nSegmentsTheta, bool drawCaps, Color color) {
	if (nSegmentsTheta < 3) return;

	int numVertex = nSegmentsTheta * (drawCaps ? 0 : 4) * 2;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	rlPushMatrix();
	rlTranslatef(cyl.pt1.x, cyl.pt1.y, cyl.pt1.z);

	//ROTATION
	Vector3 AB = Vector3Subtract(cyl.pt2, cyl.pt1);

	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(cyl.ref.q, &vect, &angle);
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
		MyDrawDiskWiresPortion({ Referential({0, 0, 0}, QuaternionMultiply(qX, qY)), 1 }, startTheta, endTheta, nSegmentsTheta, color);
		MyDrawDiskWiresPortion({ Referential({0, 1, 0}, QuaternionIdentity()), 1 }, startTheta, endTheta, nSegmentsTheta, color);
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
void MyDrawQuad2(Quad quad, Color color) {
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
void MyDrawQuadWire2(Quad quad, Color color) {
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

void MyDrawCylinder(Cylinder cyl, int nSegmentsTheta, bool drawCaps, Color color)
{
	MyDrawCylinderPortion(cyl, 0, 2 * PI, nSegmentsTheta, drawCaps, color);
}

void MyDrawCylinderWires(Cylinder cyl, int nSegmentsTheta, bool drawCaps, Color color) {
	MyDrawCylinderWiresPortion(cyl, 0, 2 * PI, nSegmentsTheta, drawCaps, color);
}

void MyDrawDisk(Disk disk, int nSegmentsTheta, Color color) {
	MyDrawDiskPortion(disk, 0, 2 * PI, nSegmentsTheta, color);
}

void MyDrawDiskWires(Disk disk, int nSegmentsTheta, Color color) {
	MyDrawDiskWiresPortion(disk, 0, 2 * PI, nSegmentsTheta, color);
}

void MyDrawDiskPortion(Disk disk, float startTheta, float endTheta, int nSegmentsTheta, Color color)
{
	if (nSegmentsTheta < 3) return;
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

void MyDrawCapsule(Capsule capsule, Color color) {
	rlPushMatrix();
	rlTranslatef(capsule.referential.origin.x, capsule.referential.origin.y, capsule.referential.origin.z);

	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(capsule.referential.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

	//rlScalef(capsule.radius, 1, capsule.radius);

	Vector3 up = { 0, capsule.height, 0 };
	Vector3 down = { 0, 0, 0 };

	Quaternion qUp = QuaternionFromAxisAngle({ 0, 0, 1 }, 0.5 * PI);
	Quaternion qDown = QuaternionFromAxisAngle({ 0, 0, 1 }, -0.5 * PI);
	Quaternion qIdentity = QuaternionIdentity();

	Referential ref = Referential(down, qIdentity);
	Cylinder cylinder = Cylinder(ref, capsule.radius, capsule.height);

	Sphere sphereUp = { up, capsule.radius };
	Sphere sphereDown = { down, capsule.radius };

	MyDrawCylinder(cylinder, 25, true, color);
	MyDrawSpherePortion(qUp, sphereUp, 0, PI, 0, PI, 30, 30, color);
	MyDrawSpherePortion(qDown, sphereDown, 0, PI, 0, PI, 25, 25, color);;

	rlPopMatrix();
}

void MyDrawCapsuleWires(Capsule capsule, Color color) {
	rlPushMatrix();
	rlTranslatef(capsule.referential.origin.x, capsule.referential.origin.y, capsule.referential.origin.z);

	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(capsule.referential.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

	//rlScalef(capsule.radius, 1, capsule.radius);

	Vector3 up = {0, capsule.height, 0};
	Vector3 down = { 0, 0, 0 };

	Quaternion qUp = QuaternionFromAxisAngle({ 0, 0, 1 }, 0.5 * PI);
	Quaternion qDown = QuaternionFromAxisAngle({ 0, 0, 1 }, 1.5 * PI);
	Quaternion qIdentity = QuaternionIdentity();

	Referential ref = Referential(down, qIdentity);
	Cylinder cylinder = Cylinder(ref, capsule.radius, capsule.height);

	Sphere sphereUp = { up, capsule.radius };
	Sphere sphereDown = { down, capsule.radius };
	
	MyDrawCylinderWires(cylinder, 25, true, color);
	MyDrawSphereWiresPortion(qUp, sphereUp, 0, PI, 0, PI, 30, 30, color);
	MyDrawSphereWiresPortion(qDown, sphereDown, 0, PI, 0, PI, 25, 25, color);

	rlPopMatrix();
}

void MyDrawRoundedBoxV2(RoundedBox roundedBox, Color color) {
	rlPushMatrix();

	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(roundedBox.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

	Vector3 posRef = Vector3Add(roundedBox.ref.origin, { -roundedBox.extension.x / 2, -roundedBox.extension.y / 2, roundedBox.extension.z / 2 });

	// référentiel de capsules
	Quaternion qLeft = QuaternionFromAxisAngle({ 1, 0, 0 }, -PI * 0.5f);
	Capsule capsLeftBottom = { Referential(posRef, qLeft), roundedBox.radius, roundedBox.extension.z };
	MyDrawCapsule(capsLeftBottom, color);

	Quaternion qFront = QuaternionFromAxisAngle({ 0, 0, 1 }, -PI * 0.5f);
	Capsule capsFrontBottom = { Referential(posRef, qFront), roundedBox.radius, roundedBox.extension.x };
	MyDrawCapsule(capsFrontBottom, color);

	Quaternion qUp = QuaternionIdentity();
	Capsule capsFrontLeft = { Referential(posRef, qUp), roundedBox.radius, roundedBox.extension.y };
	MyDrawCapsule(capsFrontLeft, color);
	// fin référentiel de capsules

	// toutes les autres capsules
	Capsule capsFrontTop = { Referential(Vector3Add(posRef, { 0, roundedBox.extension.y, 0 }), qFront), roundedBox.radius, roundedBox.extension.x };
	MyDrawCapsule(capsFrontTop, color);

	Capsule capsFrontRight = {Referential(Vector3Add(posRef, {roundedBox.extension.x, 0, 0}), qUp), roundedBox.radius, roundedBox.extension.y};
	MyDrawCapsule(capsFrontRight, color);

	Capsule capsRightBottom = { Referential(Vector3Add(posRef, { roundedBox.extension.x, 0, 0 }), qLeft), roundedBox.radius, roundedBox.extension.z };
	MyDrawCapsule(capsRightBottom, color);

	Capsule capsRightTop = { Referential(Vector3Add(posRef, { roundedBox.extension.x, roundedBox.extension.y, 0 }), qLeft), roundedBox.radius, roundedBox.extension.z };
	MyDrawCapsule(capsRightTop, color);

	Capsule capsLeftTop = { Referential(Vector3Add(posRef, { 0, roundedBox.extension.y, 0 }), qLeft), roundedBox.radius, roundedBox.extension.z };
	MyDrawCapsule(capsLeftTop, color);

	Capsule capsBackBottom = {Referential(Vector3Add(posRef, {0, 0, -roundedBox.extension.z}), qFront), roundedBox.radius, roundedBox.extension.x};
	MyDrawCapsule(capsBackBottom, color);

	Capsule capsBackTop = { Referential(Vector3Add(posRef, { 0, roundedBox.extension.y, -roundedBox.extension.z }), qFront), roundedBox.radius, roundedBox.extension.x };
	MyDrawCapsule(capsBackTop, color);

	Capsule capsBackLeft = { Referential(Vector3Add(posRef, { 0, 0, -roundedBox.extension.z }), qUp), roundedBox.radius, roundedBox.extension.y };
	MyDrawCapsule(capsBackLeft, color);

	Capsule capsBackRight = { Referential(Vector3Add(posRef, { roundedBox.extension.x, 0, -roundedBox.extension.z }), qUp), roundedBox.radius, roundedBox.extension.y };
	MyDrawCapsule(capsBackRight, color);

	Quaternion qFrontQuad = QuaternionFromAxisAngle({ 1, 0, 0 }, PI * 0.5f);
	Quad quadFront = { Referential(Vector3Add(posRef, {roundedBox.extension.x / 2, roundedBox.extension.y / 2, roundedBox.radius}), qFrontQuad), {roundedBox.extension.x, roundedBox.extension.z, roundedBox.extension.y } };
	MyDrawQuad2(quadFront, color);

	Quaternion qBackQuad = QuaternionFromAxisAngle({ 1, 0, 0 }, -PI * 0.5f);
	Quad quadBack = { Referential(Vector3Add(posRef, {roundedBox.extension.x / 2,  roundedBox.extension.y / 2, -(roundedBox.extension.z + roundedBox.radius)}), qBackQuad), {roundedBox.extension.x, roundedBox.extension.z, roundedBox.extension.y} };
	MyDrawQuad2(quadBack, color);

	Quaternion qRightQuad = QuaternionFromAxisAngle({0, 0, 1}, -PI * 0.5f);
	Quad quadRight = { Referential(Vector3Add(posRef, {roundedBox.extension.x + roundedBox.radius, roundedBox.extension.y / 2, -roundedBox.extension.z / 2}), qRightQuad), {roundedBox.extension.y, roundedBox.extension.x, roundedBox.extension.z} };
	MyDrawQuad2(quadRight, color);

	Quaternion qLeftQuad = QuaternionFromAxisAngle({ 0, 0, 1 }, PI * 0.5f);
	Quad quadLeft = { Referential(Vector3Add(posRef, {-roundedBox.radius, roundedBox.extension.y / 2, -roundedBox.extension.z / 2}), qLeftQuad), {roundedBox.extension.y, roundedBox.extension.x, roundedBox.extension.z } };
	MyDrawQuad2(quadLeft, color);

	Quaternion qTopQuad = QuaternionIdentity();
	Quad quadTop = { Referential(Vector3Add(posRef, {roundedBox.extension.x / 2, roundedBox.extension.y + roundedBox.radius, -roundedBox.extension.z / 2}), qTopQuad), {roundedBox.extension.x, roundedBox.extension.y, roundedBox.extension.z } };
	MyDrawQuad2(quadTop, color);

	Quaternion qBottomQuad = QuaternionFromAxisAngle({ 0, 0, 1 }, PI);
	Quad quadBottom = { Referential(Vector3Add(posRef, {roundedBox.extension.x / 2, -roundedBox.radius, -roundedBox.extension.z / 2}), qBottomQuad), {roundedBox.extension.x, roundedBox.extension.y, roundedBox.extension.z } };
	MyDrawQuad2(quadBottom, color);

	rlPopMatrix();
}

void MyDrawRoundBoxWiresV2(RoundedBox roundedBox, Color color) {
	rlPushMatrix();

	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(roundedBox.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

	Vector3 posRef = Vector3Add(roundedBox.ref.origin, { -roundedBox.extension.x / 2, -roundedBox.extension.y / 2, roundedBox.extension.z / 2 });

	// référentiel de capsules
	Quaternion qLeft = QuaternionFromAxisAngle({ 1, 0, 0 }, -PI * 0.5f);
	Capsule capsLeftBottom = { Referential(posRef, qLeft), roundedBox.radius, roundedBox.extension.z };
	MyDrawCapsuleWires(capsLeftBottom, color);

	Quaternion qFront = QuaternionFromAxisAngle({ 0, 0, 1 }, -PI * 0.5f);
	Capsule capsFrontBottom = { Referential(posRef, qFront), roundedBox.radius, roundedBox.extension.x };
	MyDrawCapsuleWires(capsFrontBottom, color);

	Quaternion qUp = QuaternionIdentity();
	Capsule capsFrontLeft = { Referential(posRef, qUp), roundedBox.radius, roundedBox.extension.y };
	MyDrawCapsuleWires(capsFrontLeft, color);
	// fin référentiel de capsules

	// toutes les autres capsules
	Capsule capsFrontTop = { Referential(Vector3Add(posRef, { 0, roundedBox.extension.y, 0 }), qFront), roundedBox.radius, roundedBox.extension.x };
	MyDrawCapsuleWires(capsFrontTop, color);

	Capsule capsFrontRight = { Referential(Vector3Add(posRef, {roundedBox.extension.x, 0, 0}), qUp), roundedBox.radius, roundedBox.extension.y };
	MyDrawCapsuleWires(capsFrontRight, color);

	Capsule capsRightBottom = { Referential(Vector3Add(posRef, { roundedBox.extension.x, 0, 0 }), qLeft), roundedBox.radius, roundedBox.extension.z };
	MyDrawCapsuleWires(capsRightBottom, color);

	Capsule capsRightTop = { Referential(Vector3Add(posRef, { roundedBox.extension.x, roundedBox.extension.y, 0 }), qLeft), roundedBox.radius, roundedBox.extension.z };
	MyDrawCapsuleWires(capsRightTop, color);

	Capsule capsLeftTop = { Referential(Vector3Add(posRef, { 0, roundedBox.extension.y, 0 }), qLeft), roundedBox.radius, roundedBox.extension.z };
	MyDrawCapsuleWires(capsLeftTop, color);

	Capsule capsBackBottom = { Referential(Vector3Add(posRef, {0, 0, -roundedBox.extension.z}), qFront), roundedBox.radius, roundedBox.extension.x };
	MyDrawCapsuleWires(capsBackBottom, color);

	Capsule capsBackTop = { Referential(Vector3Add(posRef, { 0, roundedBox.extension.y, -roundedBox.extension.z }), qFront), roundedBox.radius, roundedBox.extension.x };
	MyDrawCapsuleWires(capsBackTop, color);

	Capsule capsBackLeft = { Referential(Vector3Add(posRef, { 0, 0, -roundedBox.extension.z }), qUp), roundedBox.radius, roundedBox.extension.y };
	MyDrawCapsuleWires(capsBackLeft, color);

	Capsule capsBackRight = { Referential(Vector3Add(posRef, { roundedBox.extension.x, 0, -roundedBox.extension.z }), qUp), roundedBox.radius, roundedBox.extension.y };
	MyDrawCapsuleWires(capsBackRight, color);

	Quaternion qFrontQuad = QuaternionFromAxisAngle({ 1, 0, 0 }, PI * 0.5f);
	Quad quadFront = { Referential(Vector3Add(posRef, {roundedBox.extension.x / 2, roundedBox.extension.y / 2, roundedBox.radius}), qFrontQuad), {roundedBox.extension.x, roundedBox.extension.z, roundedBox.extension.y } };
	MyDrawQuadWire2(quadFront, color);

	Quaternion qBackQuad = QuaternionFromAxisAngle({ 1, 0, 0 }, -PI * 0.5f);
	Quad quadBack = { Referential(Vector3Add(posRef, {roundedBox.extension.x / 2,  roundedBox.extension.y / 2, -(roundedBox.extension.z + roundedBox.radius)}), qBackQuad), {roundedBox.extension.x, roundedBox.extension.z, roundedBox.extension.y} };
	MyDrawQuadWire2(quadBack, color);

	Quaternion qRightQuad = QuaternionFromAxisAngle({ 0, 0, 1 }, -PI * 0.5f);
	Quad quadRight = { Referential(Vector3Add(posRef, {roundedBox.extension.x + roundedBox.radius, roundedBox.extension.y / 2, -roundedBox.extension.z / 2}), qRightQuad), {roundedBox.extension.y, roundedBox.extension.x, roundedBox.extension.z} };
	MyDrawQuadWire2(quadRight, color);

	Quaternion qLeftQuad = QuaternionFromAxisAngle({ 0, 0, 1 }, PI * 0.5f);
	Quad quadLeft = { Referential(Vector3Add(posRef, {-roundedBox.radius, roundedBox.extension.y / 2, -roundedBox.extension.z / 2}), qLeftQuad), {roundedBox.extension.y, roundedBox.extension.x, roundedBox.extension.z } };
	MyDrawQuadWire2(quadLeft, color);

	Quaternion qTopQuad = QuaternionIdentity();
	Quad quadTop = { Referential(Vector3Add(posRef, {roundedBox.extension.x / 2, roundedBox.extension.y + roundedBox.radius, -roundedBox.extension.z / 2}), qTopQuad), {roundedBox.extension.x, roundedBox.extension.y, roundedBox.extension.z } };
	MyDrawQuadWire2(quadTop, color);

	Quaternion qBottomQuad = QuaternionFromAxisAngle({ 0, 0, 1 }, PI);
	Quad quadBottom = { Referential(Vector3Add(posRef, {roundedBox.extension.x / 2, -roundedBox.radius, -roundedBox.extension.z / 2}), qBottomQuad), {roundedBox.extension.x, roundedBox.extension.y, roundedBox.extension.z } };
	MyDrawQuadWire2(quadBottom, color);

	rlPopMatrix();
}

void MyDrawBox(Box box, Color color) {
	rlPushMatrix();
	rlTranslatef(box.ref.origin.x, box.ref.origin.y, box.ref.origin.z);

	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(box.ref.q, &vect, &angle);
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

	Referential referentialQuadUp = Referential({ -0.5f, 0.5f, 0 }, qUp);
	Referential referentialQuadFront = Referential({ 0, 0, 0 }, qFront);
	Referential referentialQuadBack = Referential({ -1, 0, 0 }, qBack);
	Referential referentialQuadLeft = Referential({ -0.5f, 0, 0.5f });
	Referential referentialQuadRight = Referential({ -0.5f, 0, -0.5f }, qLeft);
	Referential referentialQuadDown = Referential({ -0.5f, -0.5f, 0 }, qDown);

	Quad quadUp = { referentialQuadUp, {1, 1, 1} };
	Quad quadFront = { referentialQuadFront, {1, 1, 1} };
	Quad quadBack = { referentialQuadBack, {1, 1, 1} };
	Quad quadLeft = { referentialQuadLeft, {1, 1, 1} };
	Quad quadRight = { referentialQuadRight, {1, 1, 1} };
	Quad quadDown = { referentialQuadDown, {1, 1, 1} };

	MyDrawQuad2(quadUp, color);
	MyDrawQuad2(quadLeft, color);
	MyDrawQuad2(quadRight, color);
	MyDrawQuad2(quadFront, color);
	MyDrawQuad2(quadBack, color);
	MyDrawQuad2(quadDown, color);

	rlPopMatrix();
}

void MyDrawBoxWires(Box box, Color color) {
	rlPushMatrix();
	rlTranslatef(box.ref.origin.x, box.ref.origin.y, box.ref.origin.z);

	//ROTATION
	Vector3 vect;
	float angle;
	QuaternionToAxisAngle(box.ref.q, &vect, &angle);
	rlRotatef(angle * RAD2DEG, vect.x, vect.y, vect.z);

	rlScalef(box.extension.x, box.extension.y, box.extension.z);

	Quaternion qUp = QuaternionIdentity();
	Quaternion qLeft = QuaternionFromAxisAngle({ 1, 0, 0 }, PI * 0.5f);
	Quaternion qRight = QuaternionFromAxisAngle({ 1, 0, 0 }, PI * -0.5f);
	Quaternion qFront = QuaternionFromAxisAngle({ 0, 0, 1 }, PI * -0.5f);
	Quaternion qBack = QuaternionFromAxisAngle({ 0, 0, 1 }, PI * 0.5f);

	Quaternion qX = QuaternionFromAxisAngle({ 1, 0, 0 }, PI);
	Quaternion qY = QuaternionFromAxisAngle({ 0, 1, 0 }, 0);
	Quaternion qDown = QuaternionMultiply(qX, qY);

	Referential referentialQuadUp = Referential({ -0.5f, 0.5f, 0 }, qUp);
	Referential referentialQuadFront = Referential({ 0, 0, 0 }, qFront);
	Referential referentialQuadBack = Referential({ -1, 0, 0 }, qBack);
	Referential referentialQuadLeft = Referential({ -0.5f, 0, 0.5f }, qLeft);
	Referential referentialQuadRight = Referential({ -0.5f, 0, -0.5f }, qRight);
	Referential referentialQuadDown = Referential({ -0.5f, -0.5f, 0 }, qDown);

	Quad quadUp = { referentialQuadUp, {1, 1, 1} };
	Quad quadFront = { referentialQuadFront, {1, 1, 1} };
	Quad quadBack = { referentialQuadBack, {1, 1, 1} };
	Quad quadLeft = { referentialQuadLeft, {1, 1, 1} };
	Quad quadRight = { referentialQuadRight, {1, 1, 1} };
	Quad quadDown = { referentialQuadDown, {1, 1, 1} };

	MyDrawQuadWire2(quadUp, color);
	MyDrawQuadWire2(quadLeft, color);
	MyDrawQuadWire2(quadRight, color);
	MyDrawQuadWire2(quadFront, color);
	MyDrawQuadWire2(quadBack, color);
	MyDrawQuadWire2(quadDown, color);

	rlPopMatrix();
}
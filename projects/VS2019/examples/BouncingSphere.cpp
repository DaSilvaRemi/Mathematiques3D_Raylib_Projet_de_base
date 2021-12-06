/*******************************************************************************************
*
*   raylib [core] example - Basic window
*
*   Welcome to raylib!
*
*   To test examples, just press F6 and execute raylib_compile_execute script
*   Note that compiled executable is placed in the same folder as .c file
*
*   You can find all basic examples on C:\raylib\raylib\examples folder or
*   raylib official webpage: www.raylib.com
*
*   Enjoy using raylib. :)
*
*   This example has been created using raylib 1.0 (www.raylib.com)
*   raylib is licensed under an unmodified zlib/libpng license (View raylib.h for details)
*
*   Copyright (c) 2014 Ramon Santamaria (@raysan5)
*
********************************************************************************************/

#include "raylib.h"
#include <raymath.h>
#include "rlgl.h"
#include <math.h>
#include <float.h>
#include <vector>
#include <iostream>
#include "Tools.h"
#include "Test.h"

#if defined(PLATFORM_DESKTOP)
#define GLSL_VERSION            330
#else   // PLATFORM_RPI, PLATFORM_ANDROID, PLATFORM_WEB
#define GLSL_VERSION            100
#endif

#define EPSILON 1.e-6f


template <typename T> int sgn(T val) {
	return (T(0) < val) - (val < T(0));
}

void MyUpdateOrbitalCamera(Camera* camera, float deltaTime)
{
	static Spherical sphPos = { 20,PI / 4.f,PI / 4.f };
	static Spherical sphSpeed = { 10,.4f,.4f };
	float rhoMin = 4;
	float rhoMax = 40;

	static Vector2 prevMousePos = { 0,0 };
	Vector2 mousePos = GetMousePosition();
	Vector2 mouseVect = Vector2Subtract(mousePos, prevMousePos);
	prevMousePos = mousePos;

	Spherical sphDelta = { -GetMouseWheelMove() * sphSpeed.rho * deltaTime * 10,
	IsMouseButtonDown(MOUSE_RIGHT_BUTTON) ? mouseVect.x * sphSpeed.theta * deltaTime  * - 1 : 0,
	IsMouseButtonDown(MOUSE_RIGHT_BUTTON) ? mouseVect.y * sphSpeed.phi * deltaTime * -1 : 0 };

	Spherical newSphPos = sphPos + sphDelta;
	newSphPos = { Clamp(newSphPos.rho,rhoMin,rhoMax),
	newSphPos.theta,
	Clamp(newSphPos.phi,PI / 100.f,.99f * PI) };

	sphPos = newSphPos;

	camera->position = SphericalToCartesian(sphPos);

}

int main(int argc, char* argv[])
{
	// Initialization
	//--------------------------------------------------------------------------------------
	float screenSizeCoef = 0.7f;
	const int screenWidth = 1920 * screenSizeCoef;
	const int screenHeight = 1080 * screenSizeCoef;

	InitWindow(screenWidth, screenHeight, "Bouncy Sphere");

	SetTargetFPS(60);

	//CAMERA
	Vector3 cameraPos = { 8.0f, 15.0f, 14.0f };
	Camera camera = { 0 };
	camera.position = cameraPos;
	camera.target = { 0.0f, 0.0f, 0.0f };
	camera.up = { 0.0f, 1.0f, 0.0f };
	camera.fovy = 45.0f;
	camera.type = CAMERA_PERSPECTIVE;
	SetCameraMode(camera, CAMERA_CUSTOM);  // Set an orbital camera mode

	//TEST CONVERSION CARTESIAN->CYLINDRICAL
	/*Vector3 pos = {1,1,1};
	Cylindrical cyl = CartesianToCylindrical(pos);
	cyl = cyl + cyl;*/

	time_t t = 0;
	time(&t);
	srand(t);

	std::vector<RoundedBox> roundedBoxes;
	std::vector<Quaternion> quaternions;

	for (int i = 0; i < 5; i++) {
		Vector3 randomVect = { i + 10, random_float(1, 3), i + 10 };
		float randomAngle = random_float(0.1f, 0.5f);
		Quaternion qRandom = QuaternionFromAxisAngle(randomVect, PI * randomAngle);

		randomVect = { random_float(-20, 20) , 3 , random_float(-20, 20) };
		Referential randomRef = Referential(randomVect);
		randomRef.RotateByQuaternion(qRandom);


		float randomSize = random_float(1, 3);
		RoundedBox randomRoundedBox = { randomRef, {randomSize, randomSize, randomSize}, 0.25f };

		roundedBoxes.push_back(randomRoundedBox);
		quaternions.push_back(qRandom);
	}

	Vector3 omega = { 2.95f, 1, 0 };
	Vector3 vitesse = { 3.5f, 0, 0 };

	// Main game loop
	while (!WindowShouldClose())    // Detect window close button or ESC key
	{
		// Update
		//----------------------------------------------------------------------------------
		// TODO: Update your variables here
		//----------------------------------------------------------------------------------

		float deltaTime = GetFrameTime();
		float time = (float)GetTime();

		MyUpdateOrbitalCamera(&camera, deltaTime);

		// Draw
		//----------------------------------------------------------------------------------
		BeginDrawing();

		ClearBackground(RAYWHITE);

		BeginMode3D(camera);
		{
			Vector3 size = { 5, 5, 5 };
			RoundedBox boxBack = {Referential({-4, 0, 0}), size,  0.25f };
			RoundedBox boxFront = { Referential({ 4, 0, 0 }), size, 0.25f };
			Quaternion qFront = QuaternionFromAxisAngle({ 0, 1, 0 }, PI * 0.5f);
			Quaternion qBack = QuaternionFromAxisAngle({ 0, 1, 0 }, PI * -0.5f);
			boxFront.ref.RotateByQuaternion(qFront);
			boxBack.ref.RotateByQuaternion(qBack);
			boxFront.ref.q = qFront;
			boxBack.ref.q = qBack;

			//MyDrawRoundBox(qFront, boxBack, BLUE);
			//MyDrawRoundBoxWires(qFront, boxBack, WHITE);

			//MyDrawRoundBox(qFront, boxFront, BLUE);
			MyDrawRoundBoxWires(qFront, boxFront, WHITE);

			Vector3 nextOmega = { 3.05f, 1, 0 };// Vector3Add(omega, Vector3Scale(vitesse, deltaTime));

			Quaternion qTime = QuaternionFromAxisAngle({ 0, 0, 1 }, PI * .2f * time);

			Sphere sphere = { omega, 1 };

			//MyDrawSphereEx2(qTime, sphere, 25, 25, BLUE);
			MyDrawSphereWiresEx2(qTime, sphere, 25, 25, RED);

			Segment seg = {  omega, nextOmega };
			MyDrawSegment(QuaternionIdentity(), seg, DARKPURPLE);

			Vector3 interPt;
			Vector3 interNormal;
			bool isIntersec = IntersecRoundedBox(seg, boxBack, interPt, interNormal);
			

			if (isIntersec) {
				Vector3Reflect(vitesse, interNormal);
				DrawLine3D(interPt, interNormal, BLUE);
				DrawSphere(interPt, .2f, DARKBROWN);
			}

			isIntersec = IntersecRoundedBox(seg, boxBack, interPt, interNormal);
			

			if (isIntersec) {
				Vector3Reflect(vitesse, interNormal);
				DrawLine3D(interPt, interNormal, BLUE);
				DrawSphere(interPt, .2f, DARKBROWN);
			}

			//nextOmega = Vector3Add(omega, Vector3Scale(vitesse, deltaTime));
			//omega = nextOmega;


			/*Quaternion qTime = QuaternionFromAxisAngle({1, 0, 0}, PI * .2f * time);
			Quaternion q = QuaternionIdentity();

			Sphere sphere = { omega, 1 };
			
			MyDrawSphereEx2(qTime, sphere, 25, 25, BLUE);
			MyDrawSphereWiresEx2(qTime, sphere, 25, 25, WHITE);


			Segment seg = { omega, Vector3Add(omega, vitesse) };
			MyDrawSegment(q, seg, RED);

			Quad quad = { Referential({0, 0, 0}), {50, 1, 50} };
			MyDrawQuad2(q, quad, DARKGREEN);
			MyDrawQuadWire2(q, quad, WHITE);

			Vector3 size = { 50, 15, 1.5f };
			Box boxBack = { Referential({ 25, 7, 25}), size };
			Box boxLeft = { Referential({ 25, 7, -25 }), size };
			Box boxRight = { Referential({ -25, 7, -25 }), size };
			Box boxFront = { Referential({ 25, 7, -25 }), size };		

			Quaternion qFront = QuaternionIdentity();
			Quaternion qLeft = QuaternionFromAxisAngle({ 0, 1, 0 }, PI * 0.5f);
			Quaternion qRight = QuaternionFromAxisAngle({ 0, 1, 0 }, PI * -0.5f);

			MyDrawBox(qFront, boxBack, GRAY);
			MyDrawBoxWires(qFront, boxBack, WHITE);

			MyDrawBox(qLeft, boxLeft, GRAY);
			MyDrawBoxWires(qLeft, boxLeft, WHITE);

			MyDrawBox(qLeft, boxRight, GRAY);
			MyDrawBoxWires(qLeft, boxRight, WHITE);

			MyDrawBox(qFront, boxFront, GRAY);
			MyDrawBoxWires(qFront, boxFront, WHITE);


			Vector3 nextOmega = Vector3Add(omega, Vector3Scale(vitesse, deltaTime));

			seg = { nextOmega, Vector3Add(nextOmega, vitesse) };

			for (int i = 0; i < roundedBoxes.size(); i++) {
				MyDrawRoundBox(quaternions.at(i), roundedBoxes.at(i), GREEN);
				MyDrawRoundBoxWires(quaternions.at(i), roundedBoxes.at(i), WHITE);

				Vector3 interPt;
				Vector3 interNormal;
				bool isIntersec = IntersecRoundedBox(seg, roundedBoxes.at(i), interPt, interNormal);
				
				if (isIntersec) {
					vitesse = Vector3Reflect(vitesse, interNormal);
				}
			}

			Vector3 interPt;
			Vector3 interNormal;
			bool isIntersec = InterSegQuad(seg, quad, interPt, interNormal);

			if (isIntersec) {
				vitesse = Vector3Reflect(vitesse, interNormal);
			}

			isIntersec = IntersecBox(seg, boxBack, interPt, interNormal);

			if (isIntersec) {
				vitesse = Vector3Reflect(vitesse, interNormal);
			}

			isIntersec = IntersecBox(seg, boxLeft, interPt, interNormal);

			if (isIntersec) {
				vitesse = Vector3Reflect(vitesse, interNormal);
			}

			isIntersec = IntersecBox(seg, boxRight, interPt, interNormal);

			if (isIntersec) {
				vitesse = Vector3Reflect(vitesse, interNormal);
			}

			isIntersec = IntersecBox(seg, boxFront, interPt, interNormal);

			if (isIntersec) {
				vitesse = Vector3Reflect(vitesse, interNormal);
			}

			vitesse.y -= 1 * deltaTime;

			nextOmega = Vector3Add(omega, Vector3Scale(vitesse, deltaTime));
			omega = nextOmega;*/

			/*RoundedBox roundedBox = {Referential({0, 0, 0}), {2, 2, 2}, 0.25f};
			Segment segment = { { -4, 0, -5 } , { 4, 0, 5 } };
			MyDrawRoundBox(q, roundedBox, BLUE);
			MyDrawRoundBoxWires(q, roundedBox, WHITE);
			MyDrawSegment(q, segment, RED);

			Vector3 interPt;
			Vector3 interNormal;
			bool isIntersecRoundBox = IntersecRoundedBox(segment, roundedBox, interPt, interNormal);

			if (isIntersecRoundBox) {
				DrawLine3D(interPt, interNormal, BLUE);
				DrawSphere(interPt, .2f, DARKBROWN);
			}*/
			//
			//3D REFERENTIAL
			DrawGrid(20, 1.0f);        // Draw a grid
			DrawLine3D({ 0 }, { 0,10,0 }, DARKGRAY);
			DrawSphere({ 10,0,0 }, .2f, RED);
			DrawSphere({ 0,10,0 }, .2f, GREEN);
			DrawSphere({ 0,0,10 }, .2f, BLUE);

			//TestDisplayCylinder();
		}
		EndMode3D();

		EndDrawing();
		//----------------------------------------------------------------------------------
	}

	// De-Initialization
	//--------------------------------------------------------------------------------------  
	CloseWindow();        // Close window and OpenGL context
	//--------------------------------------------------------------------------------------

	return 0;
}
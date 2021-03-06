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

#if defined(PLATFORM_DESKTOP)
#define GLSL_VERSION            330
#else   // PLATFORM_RPI, PLATFORM_ANDROID, PLATFORM_WEB
#define GLSL_VERSION            100
#endif

#define EPSILON 1.e-6f


template <typename T>
int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

/*void MyUpdateOrbitalCamera(Camera* camera, float deltaTime)
{
	static Spherical sphPos = { 20,PI / 4.f,PI / 4.f };
	static Spherical sphSpeed = { 10,.4f,.4f };
	float rhoMin = 4;
	float rhoMax = 40;

	static Vector2 prevMousePos = { 0,0 };
	Vector2 mousePos = GetMousePosition();
	Vector2 mouseVect = Vector2Subtract(mousePos, prevMousePos);
	prevMousePos = mousePos;

	Spherical sphDelta = { -GetMouseWheelMove() * sphSpeed.rho * deltaTime,
	IsMouseButtonDown(MOUSE_RIGHT_BUTTON) ? mouseVect.x * sphSpeed.theta * deltaTime  * - 1 : 0,
	IsMouseButtonDown(MOUSE_RIGHT_BUTTON) ? mouseVect.y * sphSpeed.phi * deltaTime * -1 : 0 };

	Spherical newSphPos = sphPos + sphDelta;
	newSphPos = { Clamp(newSphPos.rho,rhoMin,rhoMax),
	newSphPos.theta,
	Clamp(newSphPos.phi,PI / 100.f,.99f * PI) };

	sphPos = newSphPos;

	camera->position = SphericalToCartesian(sphPos);

}*/

/*int main(int argc, char* argv[])
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
	Vector3 pos = { 1,1,1 };
	Cylindrical cyl = CartesianToCylindrical(pos);
	cyl = cyl + cyl;

	/*	TEST INTERSECTIONS	*/

/*Quaternion qOrient = QuaternionFromAxisAngle({ 0,0,1 }, PI * .2f);

//TEST INTERSECTION SEGMENT PLANE
Segment segment = { { -4, 0, -5 } , { 4, 0, 5 } };
Plane plane = { Vector3RotateByQuaternion({ 0,1,0 }, qOrient), 2 };

Vector3 interSectPt = { 0, 0, 0 };
Vector3 interSecNormal = { 0, 0, 0 };
bool planeHaveIntersec = InterSegPlane(segment, plane, interSectPt, interSecNormal);

//TEST INTERSECTION SEGMENT PLANE
Sphere sphere = { {0, 0, 0}, 2 };
bool sphereHaveIntersec = InterSegSphere(segment, sphere, interSectPt, interSecNormal);

//TEST INTERSECTION SEGMENT CYLYNDRE
Cylinder cylinder = { {0, 0, 0}, {0, 4, 0}, 2};

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
        //
        //3D REFERENTIAL
        DrawGrid(20, 1.0f);        // Draw a grid
        DrawLine3D({ 0 }, { 0,10,0 }, DARKGRAY);
        DrawSphere({ 10,0,0 }, .2f, RED);
        DrawSphere({ 0,10,0 }, .2f, GREEN);
        DrawSphere({ 0,0,10 }, .2f, BLUE);

        //INTERSEC BETWEEN SEGMENT AND PLANE
        /*MyDrawQuadWire2(qOrient, Vector3Scale(plane.normal, plane.d), {2, 2}, WHITE);
        MyDrawQuad2(qOrient, Vector3Scale(plane.normal, plane.d), { 2, 2 }, BLUE);

        if (planeHaveIntersec) {
            DrawSphere(interSectPt, .2f, DARKBROWN);
        }
        DrawLine3D(segment.pt1, segment.pt2, DARKGREEN);*/

// INTERSEC BETWEEN SEGMENT AND SPHERE
/*Quaternion qOrient = QuaternionFromAxisAngle({1,0,0}, PI * .5f);
MyDrawSphereEx2(qOrient, sphere, 40, 20, BLUE);
MyDrawSphereWiresEx2(qOrient, sphere, 40, 20, WHITE);

if (sphereHaveIntersec) {
    DrawLine3D(interSectPt, interSecNormal, DARKPURPLE);
    DrawSphere(interSectPt, .2f, DARKBROWN);
}
DrawLine3D(segment.pt1, segment.pt2, DARKGREEN);

/*MyDrawCylinder(qOrient, cylinder, 25, false, BLUE);
MyDrawCylinderWires(qOrient, cylinder, 25, false, WHITE);*/
/*}
EndMode3D();

EndDrawing();
//----------------------------------------------------------------------------------
}

// De-Initialization
//--------------------------------------------------------------------------------------  
CloseWindow();        // Close window and OpenGL context
//--------------------------------------------------------------------------------------

return 0;
}*/
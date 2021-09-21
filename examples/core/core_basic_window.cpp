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

#if defined(PLATFORM_DESKTOP)
#define GLSL_VERSION            330
#else   // PLATFORM_RPI, PLATFORM_ANDROID, PLATFORM_WEB
#define GLSL_VERSION            100
#endif

#define EPSILON 1.e-6f


template <typename T> int sgn(T val) {
return (T(0) < val) - (val < T(0));
}

struct Cylindrical {
	float rho;
	float theta;
	float y;

	inline Cylindrical operator+(Cylindrical a) {
		return { a.rho + rho,a.theta + theta,a.y+y };
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

struct Line {
	Vector3 pt;
	Vector3 dir;
};

struct Segment {
	Vector3 pt1;
	Vector2 pt2;
};

struct Plane {
	Vector3 normal;
	float d;
};

Cylindrical CartesianToCylindrical(Vector3 cart)
{
Cylindrical cyl;
cyl.rho = sqrtf(cart.x * cart.x + cart.z * cart.z);
cyl.y = cart.y;

if (cyl.rho < EPSILON)cyl.theta = 0;
else
{
cyl.theta = atan2f(cart.x, cart.z);
if (cyl.theta < 0)cyl.theta += PI*2;
}
return cyl;
}

Vector3 CylindricalToCartesian(Cylindrical cyl)
{
return Vector3{ cyl.rho*sinf(cyl.theta),cyl.y,cyl.rho * cosf(cyl.theta) };
}

Vector3 SphericalToCartesian(Spherical sph)
{
return Vector3{ sph.rho * sinf(sph.phi)*sinf(sph.theta),
sph.rho*cosf(sph.phi),
sph.rho * sinf(sph.phi)*cosf(sph.theta) };
}

bool InterSegPlane(Segment seg, Plane plane, Vector3& interPt, Vector3& interNormal) {
	Vector3 AB = Vector3Substrac();
}


void MyUpdateOrbitalCamera(Camera* camera, float deltaTime)
{
static Spherical sphPos = {10,PI/4.f,PI / 4.f };
static Spherical sphSpeed = { 10,.4f,.4f };
float rhoMin = 4;
float rhoMax = 40;

static Vector2 prevMousePos = { 0,0 };
Vector2 mousePos = GetMousePosition();
Vector2 mouseVect = Vector2Subtract(mousePos ,prevMousePos);
prevMousePos = mousePos;

Spherical sphDelta = { -GetMouseWheelMove()*sphSpeed.rho*deltaTime,
IsMouseButtonDown(MOUSE_RIGHT_BUTTON)? mouseVect.x * sphSpeed.theta * deltaTime:0,
IsMouseButtonDown(MOUSE_RIGHT_BUTTON) ? mouseVect.y * sphSpeed.phi * deltaTime:0 };

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
float screenSizeCoef = 1.4f;
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
Vector3 pos={ 1,1,1 };
Cylindrical cyl = CartesianToCylindrical(pos);
cyl = cyl + cyl;

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

void MyDrawSphereEx2(Vector3 centerPos, float radius, int nSegmentsTheta, int nSegmentsPhi, Color color)
{
	if (nSegmentsTheta < 3 || nSegmentsPhi < 2) return;

	std::vector<Vector3> vertexBufferTheta(nSegmentsTheta + 1);
	std::fill(vertexBufferTheta.begin(), vertexBufferTheta.end(), Vector3{ 0,radius,0 });

	int numVertex = nSegmentsTheta * nSegmentsPhi * 6;
	if (rlCheckBufferLimit(numVertex)) rlglDraw();

	rlPushMatrix();
	// NOTE: Transformation is applied in inverse order (scale -> translate)
	rlTranslatef(centerPos.x, centerPos.y, centerPos.z);
	rlScalef(radius, radius, radius);

	rlBegin(RL_TRIANGLES);
	rlColor4ub(color.r, color.g, color.b, color.a);

	float deltaPhi = PI / nSegmentsPhi;
	float deltaTheta = 2 * PI / nSegmentsTheta;

	float phi = 0;
	for (int i = 0; i < nSegmentsPhi; i++)
	{
		float theta = 0;
		Vector3 tmpBottomLeft = SphericalToCartesian(Spherical{ radius,theta,phi + deltaPhi });

		for (int j = 0; j < nSegmentsTheta; j++)
		{
			Vector3 topLeft = vertexBufferTheta[j];
			Vector3 bottomLeft = tmpBottomLeft;
			Vector3 topRight = vertexBufferTheta[j + 1];
			Vector3 bottomRight = SphericalToCartesian(Spherical{ radius,theta + deltaTheta,phi + deltaPhi });


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
		vertexBufferTheta[vertexBufferTheta.size() - 1] = vertexBufferTheta[0];
		phi += deltaPhi;
	}
	rlEnd();
	rlPopMatrix();
}
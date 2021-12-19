/*******************************************************************************************
*
*   Projet Maths 3D - 21_E3FI_3I_SI5
*
*	DA SILVA Rémi  &&  BAILLEUL Valentin
*	E3FI - 1l
*
*	Décembre 2021
*
*	Dirigé et encadré par David BILEMDJIAN
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
#include "MyDrawMethods.h"

#if defined(PLATFORM_DESKTOP)
#define GLSL_VERSION            330
#else   // PLATFORM_RPI, PLATFORM_ANDROID, PLATFORM_WEB
#define GLSL_VERSION            100
#endif

#define EPSILON 1.e-6f


template <typename T> int sgn(T val) {
	return (T(0) < val) - (val < T(0));
}

/**
*	@brief Création de la scène panoramique
*	@param camera Camera dont on set les arguments dans la fonction 'main'
*	@param deltaTime Durée que l'on set à GetFrameTime()
*/
void MyUpdateOrbitalCamera(Camera* camera, float deltaTime)
{
	static Spherical sphPos = { 20 , PI / 4 , PI / 4 };
	static Spherical sphSpeed = { 10 , 0.4f , 0.4f };
	float rhoMin = 4;
	float rhoMax = 40;

	static Vector2 prevMousePos = { 0,0 };
	Vector2 mousePos = GetMousePosition();
	Vector2 mouseVect = Vector2Subtract(mousePos, prevMousePos);
	prevMousePos = mousePos;

	Spherical sphDelta = { -GetMouseWheelMove() * sphSpeed.rho * deltaTime * 10,
		IsMouseButtonDown(MOUSE_RIGHT_BUTTON) ? mouseVect.x * sphSpeed.theta * deltaTime * -1 : 0,
		IsMouseButtonDown(MOUSE_RIGHT_BUTTON) ? mouseVect.y * sphSpeed.phi * deltaTime * -1 : 0
	};

	Spherical newSphPos = sphPos + sphDelta;

	newSphPos = { Clamp(newSphPos.rho , rhoMin , rhoMax),
		newSphPos.theta,
		Clamp(newSphPos.phi , PI / 100 , PI * 0.99f)
	};

	sphPos = newSphPos;

	camera->position = SphericalToCartesian(sphPos);
}

/**
*	@brief Méthode principale du projet
*/
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
	camera.target = { 0, 0, 0 };
	camera.up = { 0, 1, 0 };
	camera.fovy = 45.0f;
	camera.type = CAMERA_PERSPECTIVE;

	SetCameraMode(camera, CAMERA_CUSTOM);  // Set an orbital camera mode

	//		Création des RoundedBox de la scène de jeu
	// 
	// Liste des RoundedBox
	time_t t = 0;
	time(&t);
	srand(t);

	std::vector<RoundedBox> roundedBoxes;
	Quaternion q = QuaternionIdentity();
	float x = random_float(-11.0f, -9.0f);
	float y = random_float(1.0f, 5.0f);
	float z = random_float(-11.0f, -9.0f);
	for (int i = 0; i < 9; i++) {
		y = random_float(1.0f, 5.0f);
		Vector3 size = { random_float(1.7f, 3.5f), y, random_float(1.7f, 3.5f) };
		if (i != 0 && i % 3 == 0) {
			x = -10.0f;
			z += 10.0f;
		}
		Vector3 pos = { x , y , z };
		RoundedBox roundedBox = { Referential(pos, q), size, 0.5f };
		roundedBoxes.push_back(roundedBox);
		x += 10.0f;
	}

	// RoundedBox décentrées :
	y = random_float(1.0f, 5.0f);
	roundedBoxes.push_back(RoundedBox{ Referential(Vector3{ random_float(-6.0f, -4.0f), y, random_float(-6.0f, -4.0f) }, q), { random_float(1.7f, 3.5f), y,random_float(1.7f, 3.5f) }, random_float(0.3f, 0.7f) });
	y = random_float(1.0f, 5.0f);
	roundedBoxes.push_back(RoundedBox{ Referential(Vector3{ random_float(6.0f, 4.0f), y, random_float(-6.0f, -4.0f) }, q), { random_float(1.7f, 3.5f), y,random_float(1.7f, 3.5f) }, random_float(0.3f, 0.7f) });
	y = random_float(1.0f, 5.0f);
	roundedBoxes.push_back(RoundedBox{ Referential(Vector3{ random_float(-6.0f, -4.0f), y, random_float(6.0f, 4.0f) }, q), { random_float(1.7f, 3.5f), y,random_float(1.7f, 3.5f) }, random_float(0.3f, 0.7f) });
	y = random_float(1.0f, 5.0f);
	roundedBoxes.push_back(RoundedBox{ Referential(Vector3{ random_float(6.0f, 4.0f), y, random_float(6.0f, 4.0f) }, q), { random_float(1.7f, 3.5f), y,random_float(1.7f, 3.5f) }, random_float(0.3f, 0.7f) });
	//		*fin* Création des RoundedBox de la scène de jeu


	// Position initiale de la balle dans la scène de jeu
	Vector3 omega = { 1, 5, 5 };

	// Vecteur vitesse initial
	Vector3 vitesse = { 1, 0, 3 };


	// Main game loop
	while (!WindowShouldClose())    // Detect window close button or ESC key
	{
		// Update
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
			Quaternion qUp = QuaternionIdentity();
			Quaternion qTime = QuaternionFromAxisAngle({ 1, 0, 0 }, PI * .2f * time);

			// Liste des points d'intersection dans la frame actuelle :
			std::vector<Vector3> intersecPoints;

			// Balle du jeu
			Sphere sphere = { omega, 1 };
			MyDrawSphereEx2(qTime, sphere, 25, 25, BLUE);
			MyDrawSphereWiresEx2(qTime, sphere, 25, 25, WHITE);

			// Segment de détection de collision de la balle de jeu
			Vector3 nextOmega = Vector3Add(omega, Vector3Scale(vitesse, deltaTime));
			Segment seg = { omega, Vector3Add(nextOmega, Vector3Scale(vitesse, deltaTime)) };
			MyDrawSegment(qUp, seg, DARKGREEN);

			// Quaternions de rotation de l'extérieur de la scène de jeu
			Quaternion qFront = QuaternionFromAxisAngle({ 0, 0, 1 }, PI * 0.5f);
			Quaternion qBack = QuaternionFromAxisAngle({ 0, 0, 1 }, PI * -0.5f);
			Quaternion qLeft = QuaternionFromAxisAngle({ 1, 0, 0 }, PI * 0.5f);
			Quaternion qRight = QuaternionFromAxisAngle({ 1, 0, 0 }, PI * -0.5f);

			// Création mathématique de la scène de jeu
			Quad quad = { Referential({0, 0, 0}), {25, 1, 25} };
			Quad quadLeft = { Referential({ 0, 7, -12.5f }, qLeft), { 25, 1, 15 } };
			Quad quadRight = { Referential({ 0, 7, 12.5f }, qRight), { 25, 1, 15 } };
			Quad quadBack = { Referential({ -12.5f, 7, 0 }, qBack), { 15, 1, 25 } };
			Quad quadFront = { Referential({ 12.5f, 7, 0 }, qFront), { 15, 1, 25 } };

			// Dessin de la scène de jeu		// Wires de la scène
			MyDrawQuad2(quad, DARKGREEN);		MyDrawQuadWire2(quad, WHITE);
			MyDrawQuad2(quadBack, GRAY);		MyDrawQuadWire2(quadBack, WHITE);
			MyDrawQuad2(quadLeft, GRAY);		MyDrawQuadWire2(quadLeft, WHITE);
			MyDrawQuad2(quadRight, GRAY);		MyDrawQuadWire2(quadRight, WHITE);
			MyDrawQuad2(quadFront, GRAY);		MyDrawQuadWire2(quadFront, WHITE);

			//		Dessin des RoundedBox de jeu qui se trouvent dans la liste 'roundedBoxes'
			for (int i = 0; i < roundedBoxes.size(); i++) {
				MyDrawRoundedBoxV2(roundedBoxes.at(i), DARKBLUE);
				MyDrawRoundBoxWiresV2(roundedBoxes.at(i), WHITE);

				Vector3 interPt;
				Vector3 interNormal;

				bool isIntersec = IntersecSegRoundedBox(seg, roundedBoxes.at(i), interPt, interNormal);
				// Y a-t-il une intersection ?
				if (isIntersec) {
					// On tourne le vecteur vitesse en fonction de la manière dont la sphère a intersecté l'obstacle
					//////////////vitesse = Vector3Reflect(vitesse, interNormal);
					//////////////qTime = QuaternionInvert(qTime);
					intersecPoints.push_back(interPt);
					intersecPoints.push_back(interNormal);
				}
			}
			//		*fin* Dessin des RoundedBox de jeu qui se trouvent dans la liste 'roundedBoxes'

			//		Recherche de toutes les intersections avec les murs de la scène
			Vector3 interPt;
			Vector3 interNormal;


			bool isIntersec = InterSegQuad(seg, quad, interPt, interNormal);
			if (isIntersec) {
				//DrawSphere(interPt, 0.25f, DARKBROWN);
				intersecPoints.push_back(interPt);
				intersecPoints.push_back(interNormal);
			}

			isIntersec = InterSegQuad(seg, quadFront, interPt, interNormal);
			if (isIntersec) {
				//DrawSphere(interPt, 0.25f, DARKBROWN);
				intersecPoints.push_back(interPt);
				intersecPoints.push_back(interNormal);
			}

			isIntersec = InterSegQuad(seg, quadBack, interPt, interNormal);
			if (isIntersec) {
				//DrawSphere(interPt, 0.25f, DARKBROWN);
				intersecPoints.push_back(interPt);
				intersecPoints.push_back(interNormal);
			}

			isIntersec = InterSegQuad(seg, quadLeft, interPt, interNormal);
			if (isIntersec) {
				//DrawSphere(interPt, 0.25f, DARKBROWN);
				intersecPoints.push_back(interPt);
				intersecPoints.push_back(interNormal);
			}

			isIntersec = InterSegQuad(seg, quadRight, interPt, interNormal);
			if (isIntersec) {
				//DrawSphere(interPt, 0.25f, DARKBROWN);
				intersecPoints.push_back(interPt);
				intersecPoints.push_back(interNormal);
			}
			//		*fin* Recherche de toutes les intersections avec les murs de la scène


			//		Recherche de l'interPt le plus proche de la position actuelle de la sphère
			if (intersecPoints.size() > 0)
			{ // On a trouvé des points d'intersection dans la frame actuelle
				Vector3 interPtNearest = intersecPoints.at(0);
				Vector3 interNormalNearest = intersecPoints.at(1);

				for (int i = 2; i < intersecPoints.size(); i += 2)
				{
					if (Vector3Distance(omega, interPtNearest) > Vector3Distance(omega, intersecPoints.at(i))) {
						interPtNearest = intersecPoints.at(i);
						interNormalNearest = intersecPoints.at(i + 1);
					}
				}

				vitesse = Vector3Reflect(vitesse, interNormalNearest);
				qTime = QuaternionInvert(qTime);
			}
			//		*fin* Recherche de l'interPt le plus proche de la position actuelle de la sphère

			// Vider la liste quand on a checké tous les points de la frame actuelle
			intersecPoints.clear();

			vitesse.y -= 1 * deltaTime;

			// Modifications du segment de détection de collision appliqué à la sphère de jeu
			nextOmega = Vector3Add(omega, Vector3Scale(vitesse, deltaTime));
			omega = nextOmega;


			// Dessin du repère "unitaire"
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
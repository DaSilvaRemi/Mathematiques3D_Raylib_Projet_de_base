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
#include "MyDrawMethods.h"

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

void MyUpdateOrbitalCamera(Camera* camera, float deltaTime)
{
    static Spherical sphPos = {20,PI / 4.f,PI / 4.f};
    static Spherical sphSpeed = {10, .4f, .4f};
    float rhoMin = 4;
    float rhoMax = 40;

    static Vector2 prevMousePos = {0, 0};
    Vector2 mousePos = GetMousePosition();
    Vector2 mouseVect = Vector2Subtract(mousePos, prevMousePos);
    prevMousePos = mousePos;

    Spherical sphDelta = {
        -GetMouseWheelMove() * sphSpeed.rho * deltaTime * 10,
        IsMouseButtonDown(MOUSE_RIGHT_BUTTON) ? mouseVect.x * sphSpeed.theta * deltaTime * - 1 : 0,
        IsMouseButtonDown(MOUSE_RIGHT_BUTTON) ? mouseVect.y * sphSpeed.phi * deltaTime * -1 : 0
    };

    Spherical newSphPos = sphPos + sphDelta;
    newSphPos = {
        Clamp(newSphPos.rho, rhoMin, rhoMax),
        newSphPos.theta,
        Clamp(newSphPos.phi,PI / 100.f, .99f * PI)
    };

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
    Vector3 cameraPos = {8.0f, 15.0f, 14.0f};
    Camera camera = {0};
    camera.position = cameraPos;
    camera.target = {0.0f, 0.0f, 0.0f};
    camera.up = {0.0f, 1.0f, 0.0f};
    camera.fovy = 45.0f;
    camera.type = CAMERA_PERSPECTIVE;
    SetCameraMode(camera, CAMERA_CUSTOM); // Set an orbital camera mode

    std::vector<RoundedBox> roundedBoxes;
    Quaternion q = QuaternionIdentity();
    Vector3 size = {3, 3, 3};

    float x = -10;
    float z = -10;

    for (int i = 0; i < 9; i++)
    {
        if (i != 0 && i % 3 == 0)
        {
            x = -10;
            z += 10;
        }

        Vector3 pos = {x, 2, z};
        RoundedBox roundedBox = {Referential(pos, q), size, 0.5f};
        roundedBoxes.push_back(roundedBox);

        x += 10;
    }

    Vector3 omega = {1, 5, 5};
    Vector3 vitesse = {1, 0, 3};

    // Main game loop
    while (!WindowShouldClose()) // Detect window close button or ESC key
    {
        // Update
        //----------------------------------------------------------------------------------
        // TODO: Update your variables here
        //----------------------------------------------------------------------------------

        float deltaTime = GetFrameTime();
        float time = static_cast<float>(GetTime());

        MyUpdateOrbitalCamera(&camera, deltaTime);

        // Draw
        //----------------------------------------------------------------------------------
        BeginDrawing();

        ClearBackground(RAYWHITE);

        BeginMode3D(camera);
        {
            Quaternion qUp = QuaternionIdentity();
            Quaternion qTime = QuaternionFromAxisAngle({1, 0, 0}, PI * .2f * time);

            Sphere sphere = {omega, 1};

            MyDrawSphereEx2(qTime, sphere, 25, 25, BLUE);
            MyDrawSphereWiresEx2(qTime, sphere, 25, 25, WHITE);

            Vector3 nextOmega = Vector3Add(omega, Vector3Scale(vitesse, deltaTime));
            Segment seg = {omega, Vector3Add(nextOmega, Vector3Scale(vitesse, deltaTime))};
            MyDrawSegment(qUp, seg, DARKGREEN);

            Quaternion qFront = QuaternionFromAxisAngle({0, 0, 1}, PI * 0.5f);
            Quaternion qBack = QuaternionFromAxisAngle({0, 0, 1}, PI * -0.5f);
            Quaternion qLeft = QuaternionFromAxisAngle({1, 0, 0}, PI * 0.5f);
            Quaternion qRight = QuaternionFromAxisAngle({1, 0, 0}, PI * -0.5f);

            Quad quad = {Referential({0, 0, 0}), {25, 1, 25}};
            Quad quadLeft = {Referential({0, 7, -12.5f}, qLeft), {25, 1, 15}};
            Quad quadRight = {Referential({0, 7, 12.5f}, qRight), {25, 1, 15}};
            Quad quadBack = {Referential({-12.5f, 7, 0}, qBack), {15, 1, 25}};
            Quad quadFront = {Referential({12.5f, 7, 0}, qFront), {15, 1, 25}};

            MyDrawQuad2(quad, DARKGREEN);
            MyDrawQuadWire2(quad, WHITE);

            MyDrawQuad2(quadBack, GRAY);
            MyDrawQuadWire2(quadBack, WHITE);

            MyDrawQuad2(quadLeft, GRAY);
            MyDrawQuadWire2(quadLeft, WHITE);

            MyDrawQuad2(quadRight, GRAY);
            MyDrawQuadWire2(quadRight, WHITE);

            MyDrawQuad2(quadFront, GRAY);
            MyDrawQuadWire2(quadFront, WHITE);

            for (int i = 0; i < roundedBoxes.size(); i++)
            {
                MyDrawRoundedBoxV2(roundedBoxes.at(i), GREEN);
                MyDrawRoundBoxWiresV2(roundedBoxes.at(i), WHITE);

                Vector3 interPt;
                Vector3 interNormal;
                bool isIntersec = IntersecSegRoundedBox(seg, roundedBoxes.at(i), interPt, interNormal);

                if (isIntersec)
                {
                    vitesse = Vector3Reflect(vitesse, interNormal);
                    qTime = QuaternionInvert(qTime);
                }
            }

            Vector3 interPt;
            Vector3 interNormal;
            bool isIntersec = InterSegQuad(seg, quad, interPt, interNormal);

            if (isIntersec)
            {
                DrawSphere(interPt, 0.25f, DARKBROWN);
                vitesse = Vector3Reflect(vitesse, interNormal);
                qTime = QuaternionInvert(qTime);
            }

            isIntersec = InterSegQuad(seg, quadFront, interPt, interNormal);

            if (isIntersec)
            {
                DrawSphere(interPt, 0.25f, DARKBROWN);
                vitesse = Vector3Reflect(vitesse, interNormal);
                qTime = QuaternionInvert(qTime);
            }

            isIntersec = InterSegQuad(seg, quadBack, interPt, interNormal);

            if (isIntersec)
            {
                DrawSphere(interPt, 0.25f, DARKBROWN);
                vitesse = Vector3Reflect(vitesse, interNormal);
                qTime = QuaternionInvert(qTime);
            }

            isIntersec = InterSegQuad(seg, quadLeft, interPt, interNormal);

            if (isIntersec)
            {
                DrawSphere(interPt, 0.25f, DARKBROWN);
                vitesse = Vector3Reflect(vitesse, interNormal);
                qTime = QuaternionInvert(qTime);
            }

            isIntersec = InterSegQuad(seg, quadRight, interPt, interNormal);

            if (isIntersec)
            {
                DrawSphere(interPt, 0.25f, DARKBROWN);
                vitesse = Vector3Reflect(vitesse, interNormal);
                qTime = QuaternionInvert(qTime);
            }

            vitesse.y -= 1 * deltaTime;

            nextOmega = Vector3Add(omega, Vector3Scale(vitesse, deltaTime));
            omega = nextOmega;

            //
            //3D REFERENTIAL
            DrawGrid(20, 1.0f); // Draw a grid
            DrawLine3D({0}, {0, 10, 0}, DARKGRAY);
            DrawSphere({10, 0, 0}, .2f, RED);
            DrawSphere({0, 10, 0}, .2f, GREEN);
            DrawSphere({0, 0, 10}, .2f, BLUE);
        }
        EndMode3D();

        EndDrawing();
        //----------------------------------------------------------------------------------
    }

    // De-Initialization
    //--------------------------------------------------------------------------------------  
    CloseWindow(); // Close window and OpenGL context
    //--------------------------------------------------------------------------------------

    return 0;
}

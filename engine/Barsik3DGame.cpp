#include "Barsik3DGame.h"
#include <algorithm>

void Barsik3DGame::create() {

    debug = false;

    pDepthBuffer = new float[Width * Height];
qDebug() << " =============== try load tex" ;
    Texture = new BarsikSprite(QString(":/images/ttt.png"));
qDebug() << " =============== tex loaded" ;
qDebug() << " =============== try load mesh" ;
    meshCube.LoadFromObjectFile(":/obj/spyro.obj", true);
qDebug() << " =============== mesh loaded" ;
    meshCube.save();
/*
            // SOUTH
    meshCube.tris.push_back(triangle(vec3d(0.0f, 0.0f, 0.0f, 1.0f), vec3d(0.0f, 1.0f, 0.0f, 1.0f), vec3d(1.0f, 1.0f, 0.0f, 1.0f), vec2d(0.0f, 1.0f, 1.0f), vec2d(0.0f, 0.0f, 1.0f), vec2d(1.0f, 0.0f, 1.0f), Qt::white, 1.0f));
    meshCube.tris.push_back(triangle(vec3d(0.0f, 0.0f, 0.0f, 1.0f), vec3d(1.0f, 1.0f, 0.0f, 1.0f), vec3d(1.0f, 0.0f, 0.0f, 1.0f), vec2d(0.0f, 1.0f, 1.0f), vec2d(1.0f, 0.0f, 1.0f), vec2d(1.0f, 1.0f, 1.0f), Qt::white, 1.0f));

            // EAST
    meshCube.tris.push_back(triangle(vec3d(1.0f, 0.0f, 0.0f, 1.0f), vec3d(1.0f, 1.0f, 0.0f, 1.0f), vec3d(1.0f, 1.0f, 1.0f, 1.0f), vec2d(0.0f, 1.0f, 1.0f), vec2d(0.0f, 0.0f, 1.0f), vec2d(1.0f, 0.0f, 1.0f), Qt::white, 1.0f));
    meshCube.tris.push_back(triangle(vec3d(1.0f, 0.0f, 0.0f, 1.0f), vec3d(1.0f, 1.0f, 1.0f, 1.0f), vec3d(1.0f, 0.0f, 1.0f, 1.0f), vec2d(0.0f, 1.0f, 1.0f), vec2d(1.0f, 0.0f, 1.0f), vec2d(1.0f, 1.0f, 1.0f), Qt::white, 1.0f));

            // NORTH
    meshCube.tris.push_back(triangle(vec3d(1.0f, 0.0f, 1.0f, 1.0f), vec3d(1.0f, 1.0f, 1.0f, 1.0f), vec3d(0.0f, 1.0f, 1.0f, 1.0f), vec2d(0.0f, 1.0f, 1.0f), vec2d(0.0f, 0.0f, 1.0f), vec2d(1.0f, 0.0f, 1.0f), Qt::white, 1.0f));
    meshCube.tris.push_back(triangle(vec3d(1.0f, 0.0f, 1.0f, 1.0f), vec3d(0.0f, 1.0f, 1.0f, 1.0f), vec3d(0.0f, 0.0f, 1.0f, 1.0f), vec2d(0.0f, 1.0f, 1.0f), vec2d(1.0f, 0.0f, 1.0f), vec2d(1.0f, 1.0f, 1.0f), Qt::white, 1.0f));

            // WEST
    meshCube.tris.push_back(triangle(vec3d(0.0f, 0.0f, 1.0f, 1.0f), vec3d(0.0f, 1.0f, 1.0f, 1.0f), vec3d(0.0f, 1.0f, 0.0f, 1.0f), vec2d(0.0f, 1.0f, 1.0f), vec2d(0.0f, 0.0f, 1.0f), vec2d(1.0f, 0.0f, 1.0f), Qt::white, 1.0f));
    meshCube.tris.push_back(triangle(vec3d(0.0f, 0.0f, 1.0f, 1.0f), vec3d(0.0f, 1.0f, 0.0f, 1.0f), vec3d(0.0f, 0.0f, 0.0f, 1.0f), vec2d(0.0f, 1.0f, 1.0f), vec2d(1.0f, 0.0f, 1.0f), vec2d(1.0f, 1.0f, 1.0f), Qt::white, 1.0f));

            // TOP
    meshCube.tris.push_back(triangle(vec3d(0.0f, 1.0f, 0.0f, 1.0f), vec3d(0.0f, 1.0f, 1.0f, 1.0f), vec3d(1.0f, 1.0f, 1.0f, 1.0f), vec2d(0.0f, 1.0f, 1.0f), vec2d(0.0f, 0.0f, 1.0f), vec2d(1.0f, 0.0f, 1.0f), Qt::white, 1.0f));
    meshCube.tris.push_back(triangle(vec3d(0.0f, 1.0f, 0.0f, 1.0f), vec3d(1.0f, 1.0f, 1.0f, 1.0f), vec3d(1.0f, 1.0f, 0.0f, 1.0f), vec2d(0.0f, 1.0f, 1.0f), vec2d(1.0f, 0.0f, 1.0f), vec2d(1.0f, 1.0f, 1.0f), Qt::white, 1.0f));

            // BOTTOM
    meshCube.tris.push_back(triangle(vec3d(1.0f, 0.0f, 1.0f, 1.0f), vec3d(0.0f, 0.0f, 1.0f, 1.0f), vec3d(0.0f, 0.0f, 0.0f, 1.0f), vec2d(0.0f, 1.0f, 1.0f), vec2d(0.0f, 0.0f, 1.0f), vec2d(1.0f, 0.0f, 1.0f), Qt::white, 1.0f));
    meshCube.tris.push_back(triangle(vec3d(1.0f, 0.0f, 1.0f, 1.0f), vec3d(0.0f, 0.0f, 0.0f, 1.0f), vec3d(1.0f, 0.0f, 0.0f, 1.0f), vec2d(0.0f, 1.0f, 1.0f), vec2d(1.0f, 0.0f, 1.0f), vec2d(1.0f, 1.0f, 1.0f), Qt::white, 1.0f));
*/

    fTheta = 0.0f;
    fYaw = 0.0f;
    move_x = 0.0f;
    move_y = 0.0f;
    move_z = 0.0f;
    rotate_x = 0.0f;
    rotate_y = 0.0f;
    //rotate_z = 3.14f*180.0f/180.0f;
    rotate_z = 0;

    vCamera.x = 0;
    vCamera.y = 0;
    vCamera.z = 0;

    for(auto b : Key)   { b = false; }
    matProj = Matrix_MakeProjection(60.0f,
                                    (float)Height / (float)Width,
                                    0.1f, 1000.0f);
}

void Barsik3DGame::update() {

    float speed = 0.05f;

    vec3d vForward = Vector_Mul(vLookDir, speed);
/*
    if(Key[BUTTON_UP]) {
        vCamera.z += 0.05f;
    }
    if(Key[BUTTON_DOWN]) {
        vCamera.z -= 0.05f;
    }
*/

    if(Key[BUTTON_UP]) {
        vCamera = Vector_Add(vCamera,vForward);
    }
    if(Key[BUTTON_DOWN]) {
        vCamera = Vector_Sub(vCamera,vForward);
    }
    if(Key[BUTTON_LEFT]) {
        fYaw += speed;
    }
    if(Key[BUTTON_RIGHT]) {
        fYaw -= speed;
    }
    if(Key[BUTTON_ACT1]) {
        vCamera.y -= speed;
    }
    if(Key[BUTTON_ACT2]) {
        vCamera.y += speed;
    }
}

void Barsik3DGame::render(QPainter *p) {

    // Set up rotation matrices
    mat4x4 matRotZ, matRotX, matRotY, matRotXYZ;
    //fTheta += 0.015f;

    matRotX = Matrix_MakeRotationX(rotate_x);
    matRotY = Matrix_MakeRotationY(rotate_y + fTheta);
    matRotZ = Matrix_MakeRotationZ(rotate_z + fTheta);
    matRotXYZ = Matrix_MakeIdentity();
    matRotXYZ = Matrix_MultiplyMatrix(matRotZ, matRotY);
    matRotXYZ = Matrix_MultiplyMatrix(matRotXYZ, matRotX);

    mat4x4 matTrans;
    matTrans = Matrix_MakeTranslation(move_x, move_y, move_z);

    mat4x4 matWorld;
    matWorld = Matrix_MakeIdentity();
    matWorld = Matrix_MultiplyMatrix(matWorld, matRotXYZ);
    matWorld = Matrix_MultiplyMatrix(matWorld, matTrans);

    vec3d vUp = { 0.0f, 1.0f, 0.0f };
    vec3d vTarget = { 0.0f, 0.0f, 1.0f };
    mat4x4 matCameraRot = Matrix_MakeRotationY(fYaw);
    vLookDir = Matrix_MultiplyVector(matCameraRot,vTarget);
    vTarget = Vector_Add(vCamera,vLookDir);

    mat4x4 matCamera = Matrix_PointAt(vCamera, vTarget, vUp);
    mat4x4 matView   = Matrix_QuickInverse(matCamera);

    QVector<triangle> vecTrianglesToRaster;

    QPen pen = p->pen();
    QBrush brush = p->brush();

    for(auto tri : meshCube.tris)   {
        triangle triProjected, triTransformed, triViewed;

        triTransformed.p[0] = Matrix_MultiplyVector(matWorld,tri.p[0]);
        triTransformed.p[1] = Matrix_MultiplyVector(matWorld,tri.p[1]);
        triTransformed.p[2] = Matrix_MultiplyVector(matWorld,tri.p[2]);
        triTransformed.t[0] = tri.t[0];
        triTransformed.t[1] = tri.t[1];
        triTransformed.t[2] = tri.t[2];

        vec3d normal, line1, line2;
        line1  = Vector_Sub(triTransformed.p[1], triTransformed.p[0]);
        line2  = Vector_Sub(triTransformed.p[2], triTransformed.p[0]);
        normal = Vector_CrossProduct(line1, line2);
        normal = Vector_Normalise(normal);

        vec3d vCameraRay = Vector_Sub(triTransformed.p[0],vCamera);

        if(Vector_DotProduct(normal,vCameraRay) < 0.0f) {

                // lighting
            vec3d light_direction = { 0.3f, 1.0f, -0.5f };
            light_direction = Vector_Normalise(light_direction);

            float l = sqrtf(light_direction.x*light_direction.x + light_direction.y*light_direction.y + light_direction.z*light_direction.z);
            light_direction.x /= l; light_direction.y /= l; light_direction.z /= l;

                // How similar is normal to light direction
            float dp = normal.x * light_direction.x + normal.y * light_direction.y + normal.z * light_direction.z;
            triTransformed.dp = dp;
            triTransformed.color = Qt::white;

            triViewed.p[0] = Matrix_MultiplyVector(matView, triTransformed.p[0]);
            triViewed.p[1] = Matrix_MultiplyVector(matView, triTransformed.p[1]);
            triViewed.p[2] = Matrix_MultiplyVector(matView, triTransformed.p[2]);
            triViewed.color = triTransformed.color;
            triViewed.t[0] = triTransformed.t[0];
            triViewed.t[1] = triTransformed.t[1];
            triViewed.t[2] = triTransformed.t[2];

            int nClippedTriangles = 0;
            triangle clipped[2];
            nClippedTriangles = Triangle_ClipAgainstPlane(
                                { 0.0f, 0.0f, 0.1f },
                                { 0.0f, 0.0f, 1.0f },
                        triViewed, clipped[0], clipped[1]);


            for(int n = 0; n < nClippedTriangles; n++)  {

                    // from 3D to 2D
                triProjected.p[0] = Matrix_MultiplyVector(matProj,clipped[n].p[0]);
                triProjected.p[1] = Matrix_MultiplyVector(matProj,clipped[n].p[1]);
                triProjected.p[2] = Matrix_MultiplyVector(matProj,clipped[n].p[2]);
                triProjected.dp = clipped[n].dp;
                triProjected.t[0] = clipped[n].t[0];
                triProjected.t[1] = clipped[n].t[1];
                triProjected.t[2] = clipped[n].t[2];
                triProjected.t[0].u = triProjected.t[0].u / triProjected.p[0].w;
                triProjected.t[1].u = triProjected.t[1].u / triProjected.p[1].w;
                triProjected.t[2].u = triProjected.t[2].u / triProjected.p[2].w;
                triProjected.t[0].v = triProjected.t[0].v / triProjected.p[0].w;
                triProjected.t[1].v = triProjected.t[1].v / triProjected.p[1].w;
                triProjected.t[2].v = triProjected.t[2].v / triProjected.p[2].w;
                triProjected.t[0].w = 1.0f / triProjected.p[0].w;
                triProjected.t[1].w = 1.0f / triProjected.p[1].w;
                triProjected.t[2].w = 1.0f / triProjected.p[2].w;

                triProjected.p[0] = Vector_Div(triProjected.p[0],triProjected.p[0].w);
                triProjected.p[1] = Vector_Div(triProjected.p[1],triProjected.p[1].w);
                triProjected.p[2] = Vector_Div(triProjected.p[2],triProjected.p[2].w);

                    // inver X and Y
                //triProjected.p[0].x *= -1.0f;
                //triProjected.p[1].x *= -1.0f;
                //triProjected.p[2].x *= -1.0f;
                triProjected.p[0].y *= -1.0f;
                triProjected.p[1].y *= -1.0f;
                triProjected.p[2].y *= -1.0f;

                    // scale
                vec3d vOffsetView = { 1.0f, 1.0f, 0.0f };
                triProjected.p[0] = Vector_Add(triProjected.p[0], vOffsetView);
                triProjected.p[1] = Vector_Add(triProjected.p[1], vOffsetView);
                triProjected.p[2] = Vector_Add(triProjected.p[2], vOffsetView);
                triProjected.p[0].x *= 0.5f * Width;
                triProjected.p[0].y *= 0.5f * Height;
                triProjected.p[1].x *= 0.5f * Width;
                triProjected.p[1].y *= 0.5f * Height;
                triProjected.p[2].x *= 0.5f * Width;
                triProjected.p[2].y *= 0.5f * Height;

                triProjected.dp = dp;
                triProjected.color = clipped[n].color;
                vecTrianglesToRaster.push_back(triProjected);
            }
        }
    }
/*
    std::sort(vecTrianglesToRaster.begin(),vecTrianglesToRaster.end(),
         [](triangle &t1, triangle t2)
    {
        float z1 = (t1.p[0].z + t1.p[1].z + t1.p[2].z) / 3.0f;
        float z2 = (t2.p[0].z + t2.p[1].z + t2.p[2].z) / 3.0f;
        return z1 > z2;
    });
*/
    for (int i = 0; i < Width*Height; i++)
        pDepthBuffer[i] = 0.0f;

    for(auto &triToRaster : vecTrianglesToRaster)   {
        //DrawTriangle(triProjected, p, triProjected.getRasterColor());

        triangle clipped[2];
        QList<triangle> listTriangles;

        // Add initial triangle
        listTriangles.push_back(triToRaster);
        int nNewTriangles = 1;

        for (int p = 0; p < 4; p++)
        {
            int nTrisToAdd = 0;
            while (nNewTriangles > 0)
            {
                // Take triangle from front of queue
                triangle test = listTriangles.front();
                listTriangles.pop_front();
                nNewTriangles--;

                // Clip it against a plane. We only need to test each
                // subsequent plane, against subsequent new triangles
                // as all triangles after a plane clip are guaranteed
                // to lie on the inside of the plane. I like how this
                // comment is almost completely and utterly justified
                switch (p)
                {
                    case 0:	nTrisToAdd = Triangle_ClipAgainstPlane({ 0.0f, 0.0f, 0.0f }, { 0.0f, 1.0f, 0.0f }, test, clipped[0], clipped[1]); break;
                    case 1:	nTrisToAdd = Triangle_ClipAgainstPlane({ 0.0f, (float)Height - 1, 0.0f }, { 0.0f, -1.0f, 0.0f }, test, clipped[0], clipped[1]); break;
                    case 2:	nTrisToAdd = Triangle_ClipAgainstPlane({ 0.0f, 0.0f, 0.0f }, { 1.0f, 0.0f, 0.0f }, test, clipped[0], clipped[1]); break;
                    case 3:	nTrisToAdd = Triangle_ClipAgainstPlane({ (float)Width - 1, 0.0f, 0.0f }, { -1.0f, 0.0f, 0.0f }, test, clipped[0], clipped[1]); break;
                }

                // Clipping may yield a variable number of triangles, so
                // add these new ones to the back of the queue for subsequent
                // clipping against next planes
                for (int w = 0; w < nTrisToAdd; w++)
                    listTriangles.push_back(clipped[w]);
            }
            nNewTriangles = listTriangles.size();
        }


        // Draw the transformed, viewed, clipped, projected, sorted, clipped triangles
        for (auto &t : listTriangles)
        {

            TexturedTriangle(t.p[0].x, t.p[0].y, t.t[0].u, t.t[0].v, t.t[0].w,
                             t.p[1].x, t.p[1].y, t.t[1].u, t.t[1].v, t.t[1].w,
                             t.p[2].x, t.p[2].y, t.t[2].u, t.t[2].v, t.t[2].w,
                    Texture, p);

            //DrawTriangle(t, p, t.getRasterColor());

        }



    }

    p->setPen(pen);
    p->setBrush(brush);

}

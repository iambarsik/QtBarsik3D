#ifndef BARSIK3DGAME_H
#define BARSIK3DGAME_H

#include <QtWidgets>
#include <QDebug>

#include "BarsikSprite.h"

// buttons for controll
#define BUTTON_COUNT    13

#define BUTTON_UP       0
#define BUTTON_RIGHT    1
#define BUTTON_DOWN     2
#define BUTTON_LEFT     3
#define BUTTON_ACT1     4
#define BUTTON_ACT2     5
#define BUTTON_ACT3     6
#define BUTTON_ACT4     7
#define BUTTON_ACT5     8
#define BUTTON_ACT6     9
#define BUTTON_ACT7     10
#define BUTTON_ACT8     11
#define BUTTON_MENU     12

struct vec2d
{
    float u;
    float v;
    float w;
    vec2d() {
        u = 0;
        v = 0;
        w = 1;
    }
    vec2d(float u, float v, float w = 1.0f) {
        this->u = u;
        this->v = v;
        this->w = w;
    }
};

struct vec3d    {
    float x;
    float y;
    float z;
    float w;
    vec3d() {
        x = 0;
        y = 0;
        z = 0;
        w = 1;
    }
    vec3d(float x, float y, float z, float w = 1) {
        this->x = x;
        this->y = y;
        this->z = z;
        this->w = w;
    }
};

struct triangle {
    vec3d p[3];
    vec2d t[3];
    QColor color;
    float dp;
    triangle()  {}
    triangle(vec3d p1, vec3d p2, vec3d p3)  {
        this->p[0] = p1;
        this->p[1] = p2;
        this->p[2] = p3;
    }
    triangle(vec3d p1, vec3d p2, vec3d p3,
             vec2d t1, vec2d t2, vec2d t3,
             QColor color, float dp)    {
        this->p[0] = p1;
        this->p[1] = p2;
        this->p[2] = p3;
        this->t[0] = t1;
        this->t[1] = t2;
        this->t[2] = t3;
        this->color = color;
        this->dp = dp;
    }
    triangle(vec3d p1, vec3d p2, vec3d p3,
             vec2d t1, vec2d t2, vec2d t3)    {
        this->p[0] = p1;
        this->p[1] = p2;
        this->p[2] = p3;
        this->t[0] = t1;
        this->t[1] = t2;
        this->t[2] = t3;
    }
    QColor getRasterColor() {
        int r = (float)color.red()*dp;
        int g = (float)color.green()*dp;
        int b = (float)color.blue()*dp;
        if(r > 255) r = 255;
        if(g > 255) g = 255;
        if(b > 255) b = 255;
        if(r < 20) r = 20;
        if(g < 20) g = 20;
        if(b < 20) b = 20;
        QColor c(r, g, b);
        return c;
    }
};

struct mesh {
    QVector<triangle> tris;
    bool LoadFromObjectFile(QString sFilename, bool bHasTexture = false)   {
        QFile f(sFilename);
        if (!f.open(QFile::ReadOnly))
            return false;

        QVector<vec3d> verts;
        QVector<vec2d> texs;
        QTextStream in(&f);
        while (!in.atEnd()) {
            QString sLine = in.readLine();
            if(sLine[0] == 'v') {
                if(sLine[1] == ' ') {
                    vec3d v;
                    QStringList list = sLine.split(" ");
                    v.x = list[1].toFloat();
                    v.y = list[2].toFloat();
                    v.z = list[3].toFloat();
                    verts.push_back(v);
                } else if(sLine[1] == 't') {
                    vec2d v;
                    QStringList list = sLine.split(" ");
                    v.u = list[1].toFloat();
                    v.v = list[2].toFloat();
                    texs.push_back(v);
                } else if(sLine[1] == 'n') {

                }
            }
            if (!bHasTexture)   {
                if(sLine[0] == 'f') {
                    int f[3];
                    QStringList list = sLine.split(" ");
                    f[0] = list[1].toInt();
                    f[1] = list[2].toInt();
                    f[2] = list[3].toInt();
                    if(f[0] - 1 >=0 && f[0] - 1 < verts.size() &&
                            f[1] - 1 >=0 && f[1] - 1 < verts.size() &&
                            f[2] - 1 >=0 && f[2] - 1 < verts.size()) {
                        tris.push_back({ verts[f[0] - 1], verts[f[1] - 1], verts[f[2] - 1] });
                    }
                }
            } else {
                if (sLine[0] == 'f') {

                    QString tokens[8];
                    int nTokenCount = -1;

                    for(int i = 1; i < sLine.length(); i++) {
                        QChar ch = sLine[i];
                        char c = ch.toLatin1();
                        if (c == ' ' || c == '/' || c =='\n')   {
                            nTokenCount++;
                        } else {
                            if(nTokenCount < 8)
                                //tokens[nTokenCount].append(1, c);
                                tokens[nTokenCount] += ch;
                        }
                    }

                    //tokens[nTokenCount].pop_back();

                    if(tokens[0].toInt() - 1 >= 0 && tokens[0].toInt() - 1 < verts.size() &&
                       tokens[1].toInt() - 1 >= 0 && tokens[1].toInt() - 1 < verts.size() &&
                       tokens[2].toInt() - 1 >= 0 && tokens[2].toInt() - 1 < verts.size() &&
                       tokens[3].toInt() - 1 >= 0 && tokens[3].toInt() - 1 < verts.size() &&
                       tokens[4].toInt() - 1 >= 0 && tokens[4].toInt() - 1 < verts.size() &&
                            tokens[5].toInt() - 1 >= 0 && tokens[5].toInt() - 1 < verts.size())
                    {

                        tris.push_back({ verts[tokens[0].toInt() - 1], verts[tokens[2].toInt() - 1], verts[tokens[4].toInt() - 1],
                            texs[tokens[1].toInt() - 1], texs[tokens[3].toInt() - 1], texs[tokens[5].toInt() - 1] });

                    }
                    if(nTokenCount >= 7) {
                        tris.push_back({ verts[tokens[0].toInt() - 1], verts[tokens[2].toInt() - 1], verts[tokens[6].toInt() - 1],
                            texs[tokens[1].toInt() - 1], texs[tokens[3].toInt() - 1], texs[tokens[7].toInt() - 1] });
                        tris.push_back({ verts[tokens[0].toInt() - 1], verts[tokens[6].toInt() - 1], verts[tokens[4].toInt() - 1],
                            texs[tokens[1].toInt() - 1], texs[tokens[7].toInt() - 1], texs[tokens[5].toInt() - 1] });
                        tris.push_back({ verts[tokens[2].toInt() - 1], verts[tokens[4].toInt() - 1], verts[tokens[6].toInt() - 1],
                            texs[tokens[3].toInt() - 1], texs[tokens[5].toInt() - 1], texs[tokens[7].toInt() - 1] });
                    }
                            /*
                    qDebug() << "v " << tris.end()->p[0].x << " " << tris.end()->p[0].y << " " << tris.end()->p[0].z;
                    qDebug() << "v " << tris.end()->p[1].x << " " << tris.end()->p[1].y << " " << tris.end()->p[1].z;
                    qDebug() << "v " << tris.end()->p[2].x << " " << tris.end()->p[2].y << " " << tris.end()->p[2].z;

                    qDebug() << "token 0 " << tokens[0];
*/
                }
            }
        }
        return true;
    }
    void save() {
        QFile file("game_data.txt");
        file.open(QIODevice::ReadWrite);
        QTextStream out(&file);
        for(int i = 0; i < tris.size(); i++)    {
            out << "v " << tris[i].p[0].x << " " << tris[i].p[0].y << " " << tris[i].p[0].z << "\n";
            out << "v " << tris[i].p[1].x << " " << tris[i].p[1].y << " " << tris[i].p[1].z << "\n";
            out << "v " << tris[i].p[2].x << " " << tris[i].p[2].y << " " << tris[i].p[2].z << "\n";
            out << "vt " << tris[i].t[0].u << " " << tris[i].t[0].v << "\n";
            out << "vt " << tris[i].t[1].u << " " << tris[i].t[1].v << "\n";
        }
        file.close();
    }
};
struct mat4x4   {
    float m[4][4] = { 0.0f };
};



//========================================================================

class Barsik3DGame  {
public:
    float move_x;
    float move_y;
    float move_z;
    float rotate_x;
    float rotate_y;
    float rotate_z;
    vec3d vCamera;

    bool Key[BUTTON_COUNT];

public:
    void init(int w, int h) { Width = w; Height = h; create();}
    void create();
    void update();
    void render(QPainter *p);

private:
    int Width, Height;
    mesh meshCube;
    mat4x4 matProj;

    BarsikSprite *Texture;

    vec3d vLookDir;

    float fYaw;
    float fTheta;

    bool debug;

    float *pDepthBuffer = nullptr;

    vec3d Matrix_MultiplyVector(mat4x4 &m, vec3d &i)    {
        vec3d v;
        v.x = i.x * m.m[0][0] + i.y * m.m[1][0] + i.z * m.m[2][0] + i.w * m.m[3][0];
        v.y = i.x * m.m[0][1] + i.y * m.m[1][1] + i.z * m.m[2][1] + i.w * m.m[3][1];
        v.z = i.x * m.m[0][2] + i.y * m.m[1][2] + i.z * m.m[2][2] + i.w * m.m[3][2];
        v.w = i.x * m.m[0][3] + i.y * m.m[1][3] + i.z * m.m[2][3] + i.w * m.m[3][3];
        return v;
    }
    mat4x4 Matrix_MakeIdentity()    {
        mat4x4 matrix;
        matrix.m[0][0] = 1.0f;
        matrix.m[1][1] = 1.0f;
        matrix.m[2][2] = 1.0f;
        matrix.m[3][3] = 1.0f;
        return matrix;
    }
    mat4x4 Matrix_MakeRotationX(float fAngleRad)    {
        mat4x4 matrix;
        matrix.m[0][0] = 1.0f;
        matrix.m[1][1] = cosf(fAngleRad);
        matrix.m[1][2] = sinf(fAngleRad);
        matrix.m[2][1] = -sinf(fAngleRad);
        matrix.m[2][2] = cosf(fAngleRad);
        matrix.m[3][3] = 1.0f;
        return matrix;
    }
    mat4x4 Matrix_MakeRotationY(float fAngleRad)    {
        mat4x4 matrix;
        matrix.m[0][0] = cosf(fAngleRad);
        matrix.m[0][2] = sinf(fAngleRad);
        matrix.m[2][0] = -sinf(fAngleRad);
        matrix.m[1][1] = 1.0f;
        matrix.m[2][2] = cosf(fAngleRad);
        matrix.m[3][3] = 1.0f;
        return matrix;
    }
    mat4x4 Matrix_MakeRotationZ(float fAngleRad)    {
        mat4x4 matrix;
        matrix.m[0][0] = cosf(fAngleRad);
        matrix.m[0][1] = sinf(fAngleRad);
        matrix.m[1][0] = -sinf(fAngleRad);
        matrix.m[1][1] = cosf(fAngleRad);
        matrix.m[2][2] = 1.0f;
        matrix.m[3][3] = 1.0f;
        return matrix;
    }
    mat4x4 Matrix_MakeTranslation(float x, float y, float z)    {
        mat4x4 matrix;
        matrix.m[0][0] = 1.0f;
        matrix.m[1][1] = 1.0f;
        matrix.m[2][2] = 1.0f;
        matrix.m[3][3] = 1.0f;
        matrix.m[3][0] = x;
        matrix.m[3][1] = y;
        matrix.m[3][2] = z;
        return matrix;
    }
    mat4x4 Matrix_MakeProjection(float fFovDegrees, float fAspectRatio, float fNear, float fFar)    {
        float fFovRad = 1.0f / tanf(fFovDegrees * 0.5f / 180.0f * 3.14159f);
        mat4x4 matrix;
        matrix.m[0][0] = fAspectRatio * fFovRad;
        matrix.m[1][1] = fFovRad;
        matrix.m[2][2] = fFar / (fFar - fNear);
        matrix.m[3][2] = (-fFar * fNear) / (fFar - fNear);
        matrix.m[2][3] = 1.0f;
        matrix.m[3][3] = 0.0f;
        return matrix;
    }
    mat4x4 Matrix_MultiplyMatrix(mat4x4 &m1, mat4x4 &m2)    {
        mat4x4 matrix;
        for (int c = 0; c < 4; c++)
            for (int r = 0; r < 4; r++)
                matrix.m[r][c] = m1.m[r][0] * m2.m[0][c] +
                                 m1.m[r][1] * m2.m[1][c] +
                                 m1.m[r][2] * m2.m[2][c] +
                                 m1.m[r][3] * m2.m[3][c];
        return matrix;
    }
    mat4x4 Matrix_PointAt(vec3d &pos, vec3d &target, vec3d &up) {
        // Calculate new forward direction
        vec3d newForward = Vector_Sub(target, pos);
        newForward = Vector_Normalise(newForward);

        // Calculate new Up direction
        vec3d a = Vector_Mul(newForward, Vector_DotProduct(up, newForward));
        vec3d newUp = Vector_Sub(up, a);
        newUp = Vector_Normalise(newUp);

        // New Right direction is easy, its just cross product
        vec3d newRight = Vector_CrossProduct(newUp, newForward);

        // Construct Dimensioning and Translation Matrix
        mat4x4 matrix;
        matrix.m[0][0] = newRight.x;	matrix.m[0][1] = newRight.y;	matrix.m[0][2] = newRight.z;	matrix.m[0][3] = 0.0f;
        matrix.m[1][0] = newUp.x;		matrix.m[1][1] = newUp.y;		matrix.m[1][2] = newUp.z;		matrix.m[1][3] = 0.0f;
        matrix.m[2][0] = newForward.x;	matrix.m[2][1] = newForward.y;	matrix.m[2][2] = newForward.z;	matrix.m[2][3] = 0.0f;
        matrix.m[3][0] = pos.x;			matrix.m[3][1] = pos.y;			matrix.m[3][2] = pos.z;			matrix.m[3][3] = 1.0f;
        return matrix;
    }
    mat4x4 Matrix_QuickInverse(mat4x4 &m)   { // Only for Rotation/Translation Matrices
        mat4x4 matrix;
        matrix.m[0][0] = m.m[0][0]; matrix.m[0][1] = m.m[1][0]; matrix.m[0][2] = m.m[2][0]; matrix.m[0][3] = 0.0f;
        matrix.m[1][0] = m.m[0][1]; matrix.m[1][1] = m.m[1][1]; matrix.m[1][2] = m.m[2][1]; matrix.m[1][3] = 0.0f;
        matrix.m[2][0] = m.m[0][2]; matrix.m[2][1] = m.m[1][2]; matrix.m[2][2] = m.m[2][2]; matrix.m[2][3] = 0.0f;
        matrix.m[3][0] = -(m.m[3][0] * matrix.m[0][0] + m.m[3][1] * matrix.m[1][0] + m.m[3][2] * matrix.m[2][0]);
        matrix.m[3][1] = -(m.m[3][0] * matrix.m[0][1] + m.m[3][1] * matrix.m[1][1] + m.m[3][2] * matrix.m[2][1]);
        matrix.m[3][2] = -(m.m[3][0] * matrix.m[0][2] + m.m[3][1] * matrix.m[1][2] + m.m[3][2] * matrix.m[2][2]);
        matrix.m[3][3] = 1.0f;
        return matrix;
    }
    vec3d Vector_Add(vec3d &v1, vec3d &v2)  {
        return { v1.x + v2.x, v1.y + v2.y, v1.z + v2.z };
    }
    vec3d Vector_Sub(vec3d &v1, vec3d &v2)  {
        return { v1.x - v2.x, v1.y - v2.y, v1.z - v2.z };
    }
    vec3d Vector_Mul(vec3d &v1, float k)    {
        return { v1.x * k, v1.y * k, v1.z * k };
    }
    vec3d Vector_Div(vec3d &v1, float k)    {
        return { v1.x / k, v1.y / k, v1.z / k };
    }
    float Vector_DotProduct(vec3d &v1, vec3d &v2)   {
        return v1.x*v2.x + v1.y*v2.y + v1.z * v2.z;
    }
    float Vector_Length(vec3d &v)   {
        return sqrtf(Vector_DotProduct(v, v));
    }
    vec3d Vector_Normalise(vec3d &v)    {
        float l = Vector_Length(v);
        return { v.x / l, v.y / l, v.z / l };
    }
    vec3d Vector_CrossProduct(vec3d &v1, vec3d &v2) {
        vec3d v;
        v.x = v1.y * v2.z - v1.z * v2.y;
        v.y = v1.z * v2.x - v1.x * v2.z;
        v.z = v1.x * v2.y - v1.y * v2.x;
        return v;
    }
    vec3d Vector_IntersectPlane(vec3d &plane_p, vec3d &plane_n, vec3d &lineStart, vec3d &lineEnd, float &t)   {
        plane_n = Vector_Normalise(plane_n);
        float plane_d = -Vector_DotProduct(plane_n, plane_p);
        float ad = Vector_DotProduct(lineStart, plane_n);
        float bd = Vector_DotProduct(lineEnd, plane_n);
        t = (-plane_d - ad) / (bd - ad);
        vec3d lineStartToEnd = Vector_Sub(lineEnd, lineStart);
        vec3d lineToIntersect = Vector_Mul(lineStartToEnd, t);
        return Vector_Add(lineStart, lineToIntersect);
    }

    int Triangle_ClipAgainstPlane(vec3d plane_p, vec3d plane_n, triangle &in_tri, triangle &out_tri1, triangle &out_tri2)   {
        // Make sure plane normal is indeed normal
        plane_n = Vector_Normalise(plane_n);

        // Return signed shortest distance from point to plane, plane normal must be normalised
        auto dist = [&](vec3d &p)   {
            vec3d n = Vector_Normalise(p);
            return (plane_n.x * p.x + plane_n.y * p.y + plane_n.z * p.z - Vector_DotProduct(plane_n, plane_p));
        };

        // Create two temporary storage arrays to classify points either side of plane
        // If distance sign is positive, point lies on "inside" of plane
        vec3d* inside_points[3];  int nInsidePointCount = 0;
        vec3d* outside_points[3]; int nOutsidePointCount = 0;
        vec2d* inside_tex[3]; int nInsideTexCount = 0;
        vec2d* outside_tex[3]; int nOutsideTexCount = 0;

        // Get signed distance of each point in triangle to plane
        float d0 = dist(in_tri.p[0]);
        float d1 = dist(in_tri.p[1]);
        float d2 = dist(in_tri.p[2]);

        if (d0 >= 0) {
            inside_points[nInsidePointCount++] = &in_tri.p[0];
            inside_tex[nInsideTexCount++] = &in_tri.t[0];
        } else {
            outside_points[nOutsidePointCount++] = &in_tri.p[0];
            outside_tex[nOutsideTexCount++] = &in_tri.t[0];
        }
        if (d1 >= 0) {
            inside_points[nInsidePointCount++] = &in_tri.p[1];
            inside_tex[nInsideTexCount++] = &in_tri.t[1];
        } else {
            outside_points[nOutsidePointCount++] = &in_tri.p[1];
            outside_tex[nOutsideTexCount++] = &in_tri.t[1];
        }
        if (d2 >= 0) {
            inside_points[nInsidePointCount++] = &in_tri.p[2];
            inside_tex[nInsideTexCount++] = &in_tri.t[2];
        } else {
            outside_points[nOutsidePointCount++] = &in_tri.p[2];
            outside_tex[nOutsideTexCount++] = &in_tri.t[2];
        }

        // Now classify triangle points, and break the input triangle into
        // smaller output triangles if required. There are four possible
        // outcomes...

        if (nInsidePointCount == 0) {
            // All points lie on the outside of plane, so clip whole triangle
            // It ceases to exist

            return 0; // No returned triangles are valid
        }

        if (nInsidePointCount == 3)     {
            // All points lie on the inside of plane, so do nothing
            // and allow the triangle to simply pass through
            out_tri1 = in_tri;

            return 1; // Just the one returned original triangle is valid
        }

        if (nInsidePointCount == 1 && nOutsidePointCount == 2)  {
            // Triangle should be clipped. As two points lie outside
            // the plane, the triangle simply becomes a smaller triangle

            // Copy appearance info to new triangle
            //out_tri1.col =  in_tri.col;
            //out_tri1.sym = in_tri.sym;
            if(debug)
                out_tri1.color = Qt::blue;
            else
                out_tri1.color =  in_tri.color;
            out_tri1.dp =  in_tri.dp;

            // The inside point is valid, so keep that...
            out_tri1.p[0] = *inside_points[0];
            out_tri1.t[0] = *inside_tex[0];

            // but the two new points are at the locations where the
            // original sides of the triangle (lines) intersect with the plane
            float t;
            out_tri1.p[1] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[0], t);
            out_tri1.t[1].u = t * (outside_tex[0]->u - inside_tex[0]->u) + inside_tex[0]->u;
            out_tri1.t[1].v = t * (outside_tex[0]->v - inside_tex[0]->v) + inside_tex[0]->v;
            out_tri1.t[1].w = t * (outside_tex[0]->w - inside_tex[0]->w) + inside_tex[0]->w;
            out_tri1.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[1], t);
            out_tri1.t[2].u = t * (outside_tex[1]->u - inside_tex[0]->u) + inside_tex[0]->u;
            out_tri1.t[2].v = t * (outside_tex[1]->v - inside_tex[0]->v) + inside_tex[0]->v;
            out_tri1.t[2].w = t * (outside_tex[1]->w - inside_tex[0]->w) + inside_tex[0]->w;
            return 1; // Return the newly formed single triangle
        }

        if (nInsidePointCount == 2 && nOutsidePointCount == 1)  {
            // Triangle should be clipped. As two points lie inside the plane,
            // the clipped triangle becomes a "quad". Fortunately, we can
            // represent a quad with two new triangles

            // Copy appearance info to new triangles
            //out_tri1.col =  in_tri.col;
            //out_tri1.sym = in_tri.sym;
            if(debug)
                out_tri1.color =  Qt::green;
            else
                out_tri1.color =  in_tri.color;
            out_tri1.dp =  in_tri.dp;

            //out_tri2.col =  in_tri.col;
            //out_tri2.sym = in_tri.sym;
            if(debug)
                out_tri2.color =  Qt::red;
            else
                out_tri2.color =  in_tri.color;
            out_tri2.dp =  in_tri.dp;

            // The first triangle consists of the two inside points and a new
            // point determined by the location where one side of the triangle
            // intersects with the plane
            out_tri1.p[0] = *inside_points[0];
            out_tri1.p[1] = *inside_points[1];
            out_tri1.t[0] = *inside_tex[0];
            out_tri1.t[1] = *inside_tex[1];
            float t;
            out_tri1.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[0], t);
            out_tri1.t[2].u = t * (outside_tex[0]->u - inside_tex[0]->u) + inside_tex[0]->u;
            out_tri1.t[2].v = t * (outside_tex[0]->v - inside_tex[0]->v) + inside_tex[0]->v;
            out_tri1.t[2].w = t * (outside_tex[0]->w - inside_tex[0]->w) + inside_tex[0]->w;
            // The second triangle is composed of one of he inside points, a
            // new point determined by the intersection of the other side of the
            // triangle and the plane, and the newly created point above
            out_tri2.p[0] = *inside_points[1];
            out_tri2.t[0] = *inside_tex[1];
            out_tri2.p[1] = out_tri1.p[2];
            out_tri2.t[1] = out_tri1.t[2];
            out_tri2.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[1], *outside_points[0], t);
            out_tri2.t[2].u = t * (outside_tex[0]->u - inside_tex[1]->u) + inside_tex[1]->u;
            out_tri2.t[2].v = t * (outside_tex[0]->v - inside_tex[1]->v) + inside_tex[1]->v;
            out_tri2.t[2].w = t * (outside_tex[0]->w - inside_tex[1]->w) + inside_tex[1]->w;

            return 2; // Return two newly formed triangles which form a quad
        }
        return 0;
    }

    void DrawTriangle(triangle tri, QPainter *p, QColor Color)  {
        QPen pen = p->pen();
        QBrush brush = p->brush();

        p->setPen(QPen(Color,1,Qt::SolidLine));
        p->setBrush(QBrush(Color,Qt::SolidPattern));
        QVector<QPointF> polygon;
        polygon.push_back(QPointF(tri.p[0].x,tri.p[0].y));
        polygon.push_back(QPointF(tri.p[1].x,tri.p[1].y));
        polygon.push_back(QPointF(tri.p[2].x,tri.p[2].y));
        p->drawConvexPolygon(polygon);

        p->setPen(pen);
        p->setBrush(brush);
    }

    void TexturedTriangle(int x1, int y1, float u1, float v1, float w1,
                          int x2, int y2, float u2, float v2, float w2,
                          int x3, int y3, float u3, float v3, float w3,
                          BarsikSprite *tex, QPainter *p)
        {
            if (y2 < y1)
            {
                std::swap(y1, y2);
                std::swap(x1, x2);
                std::swap(u1, u2);
                std::swap(v1, v2);
                std::swap(w1, w2);
            }

            if (y3 < y1)
            {
                std::swap(y1, y3);
                std::swap(x1, x3);
                std::swap(u1, u3);
                std::swap(v1, v3);
                std::swap(w1, w3);
            }

            if (y3 < y2)
            {
                std::swap(y2, y3);
                std::swap(x2, x3);
                std::swap(u2, u3);
                std::swap(v2, v3);
                std::swap(w2, w3);
            }

            int dy1 = y2 - y1;
            int dx1 = x2 - x1;
            float dv1 = v2 - v1;
            float du1 = u2 - u1;
            float dw1 = w2 - w1;

            int dy2 = y3 - y1;
            int dx2 = x3 - x1;
            float dv2 = v3 - v1;
            float du2 = u3 - u1;
            float dw2 = w3 - w1;

            float tex_u, tex_v, tex_w;

            float dax_step = 0, dbx_step = 0,
                du1_step = 0, dv1_step = 0,
                du2_step = 0, dv2_step = 0,
                dw1_step=0, dw2_step=0;

            if (dy1) dax_step = dx1 / (float)abs(dy1);
            if (dy2) dbx_step = dx2 / (float)abs(dy2);

            if (dy1) du1_step = du1 / (float)abs(dy1);
            if (dy1) dv1_step = dv1 / (float)abs(dy1);
            if (dy1) dw1_step = dw1 / (float)abs(dy1);

            if (dy2) du2_step = du2 / (float)abs(dy2);
            if (dy2) dv2_step = dv2 / (float)abs(dy2);
            if (dy2) dw2_step = dw2 / (float)abs(dy2);

            if (dy1)
            {
                for (int i = y1; i <= y2; i++)
                {
                    int ax = x1 + (float)(i - y1) * dax_step;
                    int bx = x1 + (float)(i - y1) * dbx_step;

                    float tex_su = u1 + (float)(i - y1) * du1_step;
                    float tex_sv = v1 + (float)(i - y1) * dv1_step;
                    float tex_sw = w1 + (float)(i - y1) * dw1_step;

                    float tex_eu = u1 + (float)(i - y1) * du2_step;
                    float tex_ev = v1 + (float)(i - y1) * dv2_step;
                    float tex_ew = w1 + (float)(i - y1) * dw2_step;

                    if (ax > bx)
                    {
                        std::swap(ax, bx);
                        std::swap(tex_su, tex_eu);
                        std::swap(tex_sv, tex_ev);
                        std::swap(tex_sw, tex_ew);
                    }

                    tex_u = tex_su;
                    tex_v = tex_sv;
                    tex_w = tex_sw;

                    float tstep = 1.0f / ((float)(bx - ax));
                    float t = 0.0f;

                    for (int j = ax; j < bx; j++)
                    {
                        tex_u = (1.0f - t) * tex_su + t * tex_eu;
                        tex_v = (1.0f - t) * tex_sv + t * tex_ev;
                        tex_w = (1.0f - t) * tex_sw + t * tex_ew;
                        if (tex_w > pDepthBuffer[i*Width + j])
                        {
                            //Draw(j, i, tex->SampleGlyph(tex_u / tex_w, tex_v / tex_w), tex->SampleColour(tex_u / tex_w, tex_v / tex_w));

                            p->setPen(QPen(tex->SampleColour(tex_u / tex_w, tex_v / tex_w),1));
                            p->drawPoint(j,i);

                            pDepthBuffer[i*Width + j] = tex_w;
                        }
                        t += tstep;
                    }

                }
            }

            dy1 = y3 - y2;
            dx1 = x3 - x2;
            dv1 = v3 - v2;
            du1 = u3 - u2;
            dw1 = w3 - w2;

            if (dy1) dax_step = dx1 / (float)abs(dy1);
            if (dy2) dbx_step = dx2 / (float)abs(dy2);

            du1_step = 0, dv1_step = 0;
            if (dy1) du1_step = du1 / (float)abs(dy1);
            if (dy1) dv1_step = dv1 / (float)abs(dy1);
            if (dy1) dw1_step = dw1 / (float)abs(dy1);

            if (dy1)
            {
                for (int i = y2; i <= y3; i++)
                {
                    int ax = x2 + (float)(i - y2) * dax_step;
                    int bx = x1 + (float)(i - y1) * dbx_step;

                    float tex_su = u2 + (float)(i - y2) * du1_step;
                    float tex_sv = v2 + (float)(i - y2) * dv1_step;
                    float tex_sw = w2 + (float)(i - y2) * dw1_step;

                    float tex_eu = u1 + (float)(i - y1) * du2_step;
                    float tex_ev = v1 + (float)(i - y1) * dv2_step;
                    float tex_ew = w1 + (float)(i - y1) * dw2_step;

                    if (ax > bx)
                    {
                        std::swap(ax, bx);
                        std::swap(tex_su, tex_eu);
                        std::swap(tex_sv, tex_ev);
                        std::swap(tex_sw, tex_ew);
                    }

                    tex_u = tex_su;
                    tex_v = tex_sv;
                    tex_w = tex_sw;

                    float tstep = 1.0f / ((float)(bx - ax));
                    float t = 0.0f;

                    for (int j = ax; j < bx; j++)
                    {
                        tex_u = (1.0f - t) * tex_su + t * tex_eu;
                        tex_v = (1.0f - t) * tex_sv + t * tex_ev;
                        tex_w = (1.0f - t) * tex_sw + t * tex_ew;

                        if (tex_w > pDepthBuffer[i*Width + j])
                        {
                            //Draw(j, i, tex->SampleGlyph(tex_u / tex_w, tex_v / tex_w), tex->SampleColour(tex_u / tex_w, tex_v / tex_w));
                            p->setPen(QPen(tex->SampleColour(tex_u / tex_w, tex_v / tex_w),1));
                            p->drawPoint(j,i);

                            pDepthBuffer[i*Width + j] = tex_w;
                        }
                        t += tstep;
                    }
                }
            }
        }
};


#endif // BARSIK3DGAME_H

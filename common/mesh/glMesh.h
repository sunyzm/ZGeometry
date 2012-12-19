#ifndef _GLMESH_H_
#define _GLMESH_H_

#include <string>
#include <gl/glut.h>
#include <util/color.h>
#include "Mesh.h"

void setupObject(const CQrot& qrot, const Vector3D& trans);
void setupEye(double eyeZ);

void glDrawText(GLint x, GLint y, std::string s, GLfloat r, GLfloat g, GLfloat b);
void glDrawMesh(const CMesh* tmesh, const CQrot& rot, const Vector3D& trans, const GLfloat* clr);
void glDrawMeshColor(const CMesh* tmesh, const CQrot& rot, const Vector3D& trans, const std::vector<double>& vSigColor);
void glDrawAxes(double length);
void glDrawNumber(int nb, Vector3D& center);
void glDrawNumbers(int vid, Vector3D& center, double offset);
void glDrawSingleFeature( const CMesh* tmesh, const CQrot& rot, const Vector3D& trans, int vp);
void glDrawSignature( const CMesh* tmesh, const CQrot& rot, const Vector3D& trans, const std::vector<double>& vSigColor);
void glDrawFeatures( const CMesh* tmesh, const CQrot& rot, const Vector3D& trans, const std::vector<int>& vFeatures, int selected = -1 );

void glColorCoded(float v, float pf);
void glBlueCoded(float v, float p);
void glGreenCoded(float v, float p);
void glFalseColor(float v, float p);

#endif
#include "glMesh.h"
#include <cassert>
#include <iostream>

using namespace std;

FalseColorMap g_fcm;

void glDrawText(GLint x, GLint y, string s, GLfloat r, GLfloat g, GLfloat b)
{
	int ss = (int) s.size();
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glOrtho(0.0, glutGet(GLUT_WINDOW_WIDTH), 0.0, glutGet(GLUT_WINDOW_HEIGHT), -1.0, 1.0);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glColor3f(r,g,b);
	glRasterPos2i(x, y);
	for(int i = 0; i < ss; i++)
	{
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, s[i]);
	}
	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
}

void glDrawAxes(double length)
{
	glBegin(GL_LINES);
	glColor3f(1.0, 0.0, 0.0);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(length, 0.0, 0.0);
	glColor3f(0.0, 1.0, 0.0);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(0.0, length, 0.0);
	glColor3f(0.0, 0.0, 1.0);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(0.0, 0.0, length); 
	glEnd(); 
}

void glColorCoded(float v, float pf)
{
	int ic = v;
	float f = v - ic;
	switch(ic)
	{
	case 0: glColor4f(1,f,0,pf); break;     // red -> yellow
	case 1: glColor4f(1-f,1,0,pf); break;   // yellow -> green
	case 2: glColor4f(0,1,f,pf); break;		// green -> cyan
	case 3: glColor4f(0,1-f,1,pf); break;	// cyan -> blue
	case 4: glColor4f(f,0,1,pf); break;     // blue -> purple 
	case 5: glColor4f(1,0,1-f,pf); break;   // purple -> red
	}
}

void glFalseColor(float v, float p)
{
	int floor = v * 255.0;
	glColor4f(g_fcm.RedMap[floor], g_fcm.GreenMap[floor], g_fcm.BlueMap[floor], p);
}

void glBlueCoded(float v, float p)
{
	if(v<0.2) glColor4f(0.1, 0.1, v*3.5+0.2, p);
	else if(v<0.4) glColor4f(0.1, (v-0.2)+0.1, 0.9, p);
	else glColor4f((v-0.4)*1.3+0.1, (v-0.2)+0.1, 0.9, p);
}

void glGreenCoded(float v, float p)
{
	//(0.1,0.2,0.1) -> (0.8,1.0,0.8)
	if(v<0.4) glColor4f(0.1, v*2.0+0.2, 0.1, p);
	else if(v<0.8) glColor4f(0.1, 1.0, (v-0.4)+0.1, p);
	else glColor4f((v-0.8)*3.0+0.1, 1.0, (v-0.4)+0.1, p);
	//glColor4f(0.3,v,0.3,p);
}

void glDrawNumber(int nb, Vector3D& center)
{
	//	glColor3f(0.5,0.5,0.5);
	glRasterPos3d(center.x, center.y, center.z);
	switch (nb)
	{
	case 0: glutBitmapCharacter(GLUT_BITMAP_8_BY_13, '0'); break;
	case 1: glutBitmapCharacter(GLUT_BITMAP_8_BY_13, '1'); break;
	case 2: glutBitmapCharacter(GLUT_BITMAP_8_BY_13, '2'); break;
	case 3: glutBitmapCharacter(GLUT_BITMAP_8_BY_13, '3'); break;
	case 4: glutBitmapCharacter(GLUT_BITMAP_8_BY_13, '4'); break;
	case 5: glutBitmapCharacter(GLUT_BITMAP_8_BY_13, '5'); break;
	case 6: glutBitmapCharacter(GLUT_BITMAP_8_BY_13, '6'); break;
	case 7: glutBitmapCharacter(GLUT_BITMAP_8_BY_13, '7'); break;
	case 8: glutBitmapCharacter(GLUT_BITMAP_8_BY_13, '8'); break;
	case 9: glutBitmapCharacter(GLUT_BITMAP_8_BY_13, '9'); break;
	case 10: glutBitmapCharacter(GLUT_BITMAP_8_BY_13, 'A'); break;
	case 11: glutBitmapCharacter(GLUT_BITMAP_8_BY_13, 'B'); break;
	case 12: glutBitmapCharacter(GLUT_BITMAP_8_BY_13, 'C'); break;
	case 13: glutBitmapCharacter(GLUT_BITMAP_8_BY_13, 'D'); break;
	case 14: glutBitmapCharacter(GLUT_BITMAP_8_BY_13, 'E'); break;
	case 15: glutBitmapCharacter(GLUT_BITMAP_8_BY_13, 'F'); break;
	}
}

void glDrawNumbers(int vid, Vector3D& center, double offset)
{
	int id = vid;
	Vector3D vec = Vector3D(center.x-offset,center.y,center.z);
	while(id>0)
	{
		glDrawNumber(id % 10, vec);
		id /= 10;
		vec += Vector3D(-offset,0,0); 
	}
}

void glDrawMesh(const CMesh* tmesh, const CQrot& rot, const Vector3D& trans, const GLfloat* clr)
{
	if(!tmesh) return;

	float specReflection[] = { 0.0f, 0.0f, 0.0f, 1.0f };
	glMaterialfv(GL_FRONT, GL_SPECULAR, specReflection);
	//glMateriali(GL_FRONT, GL_SHININESS, 96);
	glPushMatrix();

	setupObject(rot, trans);

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1.0, 1.0);


	{	// just display mesh in single color
		glBegin(GL_TRIANGLES);
		for (int i = 0; i < tmesh->m_nFace; i++)
		{
			if(!tmesh->m_pFace[i].m_piEdge) continue;
			for (int j = 0; j < 3; j++)
			{
				glColor4f(clr[0], clr[1], clr[2], 1.0); 
				int pi = tmesh->m_pFace[i].m_piVertex[j];
				Vector3D norm = tmesh->m_pVertex[pi].getNormal();
				glNormal3f(norm.x, norm.y, norm.z);
				Vector3D vt = tmesh->m_pVertex[pi].m_vPosition;
				//vt -= tmesh->m_Center;
				glVertex3f(vt.x, vt.y, vt.z);
			}
		}
		glEnd();
	}

	glDisable(GL_POLYGON_OFFSET_FILL);
	glDisable(GL_LIGHTING);

	//display boundary edge in black
	//	glLineWidth(2.0);
	glBegin(GL_LINES);	
	for(int i = 0; i < tmesh->m_nEdge; i++)
	{
		if(tmesh->m_pEdge[i].m_iTwinEdge < 0) 
		{
			int p1 = tmesh->m_pEdge[i].m_iVertex[0];
			int p2 = tmesh->m_pEdge[i].m_iVertex[1];
			glLineWidth(2.0);
			glColor4f(0.0, 0.0, 0.0, 1.0);
			if(tmesh->m_pVertex[p1].m_bIsHole) 
			{
				glColor4f(0.0, 0.0, 1.0, 1.0);		//show blue edge on holes
				//glLineWidth(1.0);
			}
			Vector3D v1 = tmesh->m_pVertex[p1].m_vPosition;
			//v1 -= tmesh->m_Center;
			Vector3D v2 = tmesh->m_pVertex[p2].m_vPosition;
			//v2 -= tmesh->m_Center;
			glVertex3d(v1.x, v1.y, v1.z);
			glVertex3d(v2.x, v2.y, v2.z);
		}
	}
	glEnd();

	glEnable(GL_LIGHTING);
	glPopMatrix();
}

void glDrawMeshColor(const CMesh* tmesh, const CQrot& rot, const Vector3D& trans, const std::vector<double>& vSigColor)
{
	if(!tmesh) return;

	const int shapeSize = tmesh->getVerticesNum();
	assert( shapeSize == vSigColor.size());

	float specReflection[] = { 0.0f, 0.0f, 0.0f, 1.0f };
	glMaterialfv(GL_FRONT, GL_SPECULAR, specReflection);
	//glMateriali(GL_FRONT, GL_SHININESS, 96);
	glPushMatrix();

	setupObject(rot, trans);

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1.0, 1.0);

	{	//display signature value in false color
	 	glBegin(GL_TRIANGLES);
	 	for (int i = 0; i < tmesh->m_nFace; i++)
	 	{
	 		if(!tmesh->m_pFace[i].m_piEdge) continue;
	 		for (int j = 0; j < 3; j++)
	 		{
	 			int pi = tmesh->m_pFace[i].m_piVertex[j];
	 			float scaleVal = vSigColor[pi];
	 			glFalseColor(scaleVal, 1.0f);				
	 				
	 			Vector3D norm = tmesh->m_pVertex[pi].getNormal();
	 			glNormal3f(norm.x, norm.y, norm.z);
	 				
	 			Vector3D vt = tmesh->m_pVertex[pi].m_vPosition;
	 			//vt -= tmesh->m_Center;
	 			glVertex3f(vt.x, vt.y, vt.z);
	 		}
	 	}
	 	glEnd();
	}

	glDisable(GL_POLYGON_OFFSET_FILL);
	glDisable(GL_LIGHTING);

	//display boundary edge in black
	//	glLineWidth(2.0);
	glBegin(GL_LINES);	
	for(int i = 0; i < tmesh->m_nEdge; i++)
	{
		if(tmesh->m_pEdge[i].m_iTwinEdge < 0) 
		{
			int p1 = tmesh->m_pEdge[i].m_iVertex[0];
			int p2 = tmesh->m_pEdge[i].m_iVertex[1];
			glLineWidth(2.0);
			glColor4f(0.0, 0.0, 0.0, 1.0);
			if(tmesh->m_pVertex[p1].m_bIsHole) 
			{
				glColor4f(0.0, 0.0, 1.0, 1.0);		//show blue edge on holes
				//glLineWidth(1.0);
			}
			Vector3D v1 = tmesh->m_pVertex[p1].m_vPosition;
			//v1 -= tmesh->m_Center;
			Vector3D v2 = tmesh->m_pVertex[p2].m_vPosition;
			//v2 -= tmesh->m_Center;
			glVertex3d(v1.x, v1.y, v1.z);
			glVertex3d(v2.x, v2.y, v2.z);
		}
	}
	glEnd();

	glEnable(GL_LIGHTING);
	glPopMatrix();
}

void glDrawSignature( const CMesh* tmesh, const CQrot& rot, const Vector3D& trans, const vector<double>& vSigColor)
{
	glPushMatrix();
	setupObject(rot, trans);

	int shapeSize = tmesh->getVerticesNum();
	for (int vi = 0; vi < shapeSize; ++vi)
	{
		double col = vSigColor[vi];
		if (col < 0 || col > 1) { cerr << "invalid color!" << endl; continue; }
		Vector3D vt = tmesh->m_pVertex[vi].m_vPosition;
		glFalseColor((float)col, 1.0);
		glPushMatrix();
		glTranslated(vt.x, vt.y, vt.z);
		glutSolidSphere(tmesh->m_edge*0.6, 4, 4);
		glPopMatrix();
	}
	glPopMatrix();
}

void setupObject(const CQrot& qrot, const Vector3D& trans)
{
	glTranslated(trans.x, trans.y, trans.z);
	double rot[16];
	qrot.convert( rot );
	glMultMatrixd(( GLdouble*)  rot );
}

void setupEye(double eyeZ)	//eyeZ should be g_EyeZ
{
	glLoadIdentity();
	gluLookAt(0, 0, eyeZ, 0, 0, 0, 0, 1, 0);
}

void glDrawFeatures( const CMesh* tmesh, const CQrot& rot, const Vector3D& trans, const std::vector<int>& vFeatures, int selected /*= -1 */ )
{
	if (!tmesh) return;

	glPushMatrix();
	setupObject(rot, trans);

	int featureSize = (int)vFeatures.size();
	for (int i = 0; i < featureSize; ++i)
	{
		Vector3D vt = tmesh->m_pVertex[vFeatures[i]].m_vPosition;
		if (selected == i)
			glColor4f(1.0f, 0.5f, 0.0f, 1.0f);
		else 
			glColor4f(1.0f, 1.0f, 0.0f, 1.0f);
		glPushMatrix();
		glTranslated(vt.x, vt.y, vt.z);
		glutSolidSphere(tmesh->m_edge*2, 8, 8);
		glPopMatrix();
	}
	glPopMatrix();
}

void glDrawSingleFeature( const CMesh* tmesh, const CQrot& rot, const Vector3D& trans, int vp) 
{
	if (!tmesh) return;

	glPushMatrix();
	setupObject(rot, trans);

	Vector3D vt = tmesh->m_pVertex[vp].m_vPosition;
	glColor4f(1.0f, 0.5f, 0.0f, 1.0f);
	glPushMatrix();
	glTranslated(vt.x, vt.y, vt.z);
	glutSolidSphere(tmesh->m_edge*2, 8, 8);
	glPopMatrix();
	glPopMatrix();
}
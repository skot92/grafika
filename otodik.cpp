#include <GLFW\glfw3.h>                           // (or others, depending on the system in use)
#include <cmath>
#include <iostream>
#include <stdio.h>

#define PI 3.141592653589793

GLdouble updateFrequency = 0.01, lastUpdate;

/*======================================*/

/**
* 3D vektor
*/
typedef struct { GLdouble x, y, z; } VECTOR3;

/**
* 4D vektor
*/
typedef struct { GLdouble x, y, z, w; } VECTOR4;

/**
* 4x4 m�trix, sorfolytonosan t�rolva
*/
typedef GLdouble MATRIX4[4][4];

/**
* egy laphoz tartoz� cs�csok indexei
*/
typedef struct { GLint v[4]; } FACE;

typedef GLdouble MATRIX7[7][7];

/*======================================*/

/**
* visszadja az (x,y,z) vektort
*/
VECTOR3 initVector3(GLdouble x, GLdouble y, GLdouble z) {
	VECTOR3 P;

	P.x = x;
	P.y = y;
	P.z = z;

	return P;
}

/**
* visszadja az (x,y,z,w) vektort
*/
VECTOR4 initVector4(GLdouble x, GLdouble y, GLdouble z, GLdouble w) {
	VECTOR4 P;

	P.x = x;
	P.y = y;
	P.z = z;
	P.w = w;

	return P;
}

VECTOR3 convertToInhomogen(VECTOR4 vector) {
	return initVector3(vector.x / vector.w, vector.y / vector.w, vector.z / vector.w);
}

VECTOR4 convertToHomogen(VECTOR3 vector) {
	return initVector4(vector.x, vector.y, vector.z, 1.0);
}

GLdouble getVectorLength(VECTOR3 vec) {
	return sqrt(pow(vec.x, 2) + pow(vec.y, 2) + pow(vec.z, 2));
}

VECTOR3 normalizeVector(VECTOR3 vector) {
	GLdouble length = getVectorLength(vector);
	return initVector3(vector.x / length, vector.y / length, vector.z / length);
}

GLdouble dotProduct(VECTOR3 a, VECTOR3 b) {
	return (a.x * b.x + a.y * b.y + a.z * b.z);
}

VECTOR3 crossProduct(VECTOR3 a, VECTOR3 b) {

	return initVector3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}


VECTOR3 vecSub(VECTOR3 a, VECTOR3 b) {
	return initVector3(
		b.x - a.x,
		b.y - a.y,
		b.z - a.z);
}

/**
* visszaadja az 'a' vektor hossz�t
*/
float length(VECTOR3 a) {
	return sqrt(dotProduct(a, a));
}

/**
* visszaadja az 'a' vektor normaliz�ltj�t
*/
VECTOR3 normalize(VECTOR3 a) {
	float len = length(a);

	return initVector3(a.x / len, a.y / len, a.z / len);
}



/**
* visszadja a (v0,v1,v2,v3) indexekhez tartoz� lapot
*/
FACE initFace(GLint v0, GLint v1, GLint v2, GLint v3) {
	FACE f;

	f.v[0] = v0;
	f.v[1] = v1;
	f.v[2] = v2;
	f.v[3] = v3;

	return f;
}






/**
* felt�lti az A m�trixot az egys�gm�trixszal
*/
void initIdentityMatrix(MATRIX4 A)
{
	int i, j;

	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			A[i][j] = 0.0f;

	for (i = 0; i < 4; i++)
		A[i][i] = 1.0f;
}

void initRotationMatrixX(MATRIX4 A, float alpha)
{
	float c = cos(alpha);
	float s = sin(alpha);

	initIdentityMatrix(A);

	A[1][1] = c;
	A[1][2] = -s;
	A[2][1] = s;
	A[2][2] = c;
}

void initRotationMatrixY(MATRIX4 A, float alpha)
{
	float c = cos(alpha);
	float s = sin(alpha);

	initIdentityMatrix(A);

	A[0][0] = c;
	A[0][2] = s;
	A[2][0] = -s;
	A[2][2] = c;
}

void initRotationMatrixZ(MATRIX4 A, float alpha)
{
	float c = cos(alpha);
	float s = sin(alpha);

	initIdentityMatrix(A);

	A[0][0] = c;
	A[0][1] = -s;
	A[1][0] = s;
	A[1][1] = c;
}

/**
* felt�lti az A m�trixot a (0,0,s) k�z�ppont� centr�lis vet�t�s m�trixsz�val
*/
void initPersProjMatrix(MATRIX4 A, GLdouble s)
{
	initIdentityMatrix(A);

	A[2][2] = 0.0f;
	A[3][2] = -1.0f / s;
}

/**
* felt�lti az A m�trixot a (wx,wy):(wx+ww,wy+wh) ablakb�l (vx,vy):(vx+vw,vy+vh) n�zetbe
* t�rt�n� transzform�ci� m�trixsz�val
*/
void initWtvMatrix(MATRIX4 A, GLdouble wx, GLdouble wy, GLdouble ww, GLdouble wh,
	GLdouble vx, GLdouble vy, GLdouble vw, GLdouble vh)
{
	initIdentityMatrix(A);

	A[0][0] = vw / ww;
	A[1][1] = vh / wh;
	A[0][3] = -wx*(vw / ww) + vx;
	A[1][3] = -wy*(vh / wh) + vy;
}

void initViewMatrix(MATRIX4 A, VECTOR3 eye, VECTOR3 center, VECTOR3 up) {
	initIdentityMatrix(A);

	// a kamerabol az iranyt meghatarozo pontba mutato vektor
	// center az a pont,amely fele a kamerat tartjuk, eye a kamera hely�t adja
	VECTOR3 centerMinusEye = initVector3(center.x - eye.x, center.y - eye.y, center.z - eye.z);

	// a fenti vektor -1 szeres�nek egyseg hosszra norm�ltja, ebbol lesz a kamera rendszer�nek z-tengelye
	VECTOR3 f = normalizeVector(initVector3(-centerMinusEye.x, -centerMinusEye.y, -centerMinusEye.z));

	// az up vektor es a leendo z tengelyirany vektorialis szorzata adja a kamera x-tengelyiranyat
	VECTOR3 s = normalizeVector(crossProduct(up, f));

	// a kamera y tengelyiranya a mar elkeszult z irany es x irany vektorialis szorzatak�nt j�n ki
	VECTOR3 u = crossProduct(f, s);

	A[0][0] = s.x;
	A[0][1] = s.y;
	A[0][2] = s.z;
	A[1][0] = u.x;
	A[1][1] = u.y;
	A[1][2] = u.z;
	A[2][0] = f.x;
	A[2][1] = f.y;
	A[2][2] = f.z;
	A[0][3] = -dotProduct(s, eye);
	A[1][3] = -dotProduct(u, eye);
	A[2][3] = -dotProduct(f, eye);
}




/**
* visszaadja az A m�trix �s a v vektor szorzat�t, A*v-t
*/
VECTOR4 mulMatrixVector(MATRIX4 A, VECTOR4 v) {
	return initVector4(
		A[0][0] * v.x + A[0][1] * v.y + A[0][2] * v.z + A[0][3] * v.w,
		A[1][0] * v.x + A[1][1] * v.y + A[1][2] * v.z + A[1][3] * v.w,
		A[2][0] * v.x + A[2][1] * v.y + A[2][2] * v.z + A[2][3] * v.w,
		A[3][0] * v.x + A[3][1] * v.y + A[3][2] * v.z + A[3][3] * v.w);
}

/**
* felt�lti a C m�trixot az A �s B m�trixok szorzat�val, A*B-vel
*/
void mulMatrices(MATRIX4 A, MATRIX4 B, MATRIX4 C) {
	int i, j, k;

	GLdouble sum;
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++) {
			sum = 0;
			for (k = 0; k < 4; k++)
				sum = sum + A[i][k] * B[k][j];
			C[i][j] = sum;
		}
}

/*======================================*/

/**
* ablak m�rete
*/
GLdouble winWidth = 800.0f, winHeight = 600.0f;


MATRIX4 view;

/**
* mer�leges �s centr�lis vet�t�sek m�trixai
*/
MATRIX4 Vc;

/**
* Wtv m�trixok
*/
MATRIX4 wtv;

MATRIX4 TcTorusX,T;

GLdouble alpha = 0.0f, alphaFel = 0.0f, deltaAlpha = PI / 80.0f;
GLdouble forog = 0.0f;

/**
* egys�gkocka cs�csai
*/
VECTOR4 identityCube[8] =
{
	initVector4(0.5f, 0.5f,-0.5f, 1.0f),
	initVector4(-0.5f, 0.5f,-0.5f, 1.0f),
	initVector4(-0.5f,-0.5f,-0.5f, 1.0f),
	initVector4(0.5f,-0.5f,-0.5f, 1.0f),
	initVector4(0.5f, 0.5f, 0.5f, 1.0f),
	initVector4(-0.5f, 0.5f, 0.5f, 1.0f),
	initVector4(-0.5f,-0.5f, 0.5f, 1.0f),
	initVector4(0.5f,-0.5f, 0.5f, 1.0f),
};

/**
* egys�gkocka lapjainak indexei
*/
FACE faces[6] =
{
	initFace(0, 4, 7, 3),
	initFace(1, 5, 4, 0),
	initFace(4, 5, 6, 7),
	initFace(3, 2, 1, 0),
	initFace(6, 5, 1, 2),
	initFace(7, 6, 2, 3),
};

/**
* centr�lis vet�t�s k�z�ppontj�nak Z koordin�t�ja
*/
GLdouble center = 15.0f;
/**
* n�zet koordin�t�i
*/
GLdouble cX = (winWidth - winHeight) / 2.0f, cY = 0.0f, cW = winHeight, cH = winHeight;
/**
* el��ll�tja a sz�ks�ges m�trixokat
*/

VECTOR4 transformedCube[8];

void initTransformations()
{
	MATRIX4 tmp1, tmp2, rx;
	// vet�t�si m�trixok
	initPersProjMatrix(Vc, center);
	// Wtv m�trixok
	initWtvMatrix(wtv, -4.0f, -4.0f, 8.0f, 8.0f, cX, cY, cW, cH);


	for (int i = 0; i < 8; i++) {
		transformedCube[i] = mulMatrixVector(view, identityCube[i]);
	}

	// vet�t�s * kamera
	mulMatrices(Vc, view, tmp1);

	// Wtv * vet�t�s * kamera
	mulMatrices(wtv, tmp1, T);


}



VECTOR3 eye, up, centerVec;

double R = 37;

void init()
{
	glClearColor(1.0f, 1.0f, 1.0f, 0.0f);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0.0f, winWidth, 0.0f, winHeight, 0.0f, 1.0f);

	eye = initVector3(R*cos(alpha), alphaFel, R*sin(alpha)); //megadja a kamera poz�ci�j�t (Ez legyen most a z tengely pozit�v fel�n)
	centerVec = initVector3(0.0, 0.0, 0.0); //megadja, hogy merre n�z a kamera (Ez legyen most az orig�)
	up = initVector3(0.0, 1.0, 0.0); //megdja, hogy merre van a felfele ir�ny (Ez legyen most az y tengely)

	initViewMatrix(view, eye, centerVec, up);

	initTransformations();


}

GLdouble sulypont(FACE* face) {
	GLdouble sum = 0;
	VECTOR3 asd = initVector3(0, 0, 0);
	for (int i = 0; i < 4; i++)
	{
		sum += transformedCube[face->v[i]].z;
	}

	return sum / 4;
}

int comparePointsZ(const void *a, const void *b) {
	double az, bz;
	az = sulypont((FACE*)a);
	bz = sulypont((FACE*)b);

	if (az < bz) return -1;
	if (az == bz) return  0;

	return  1;
}


void drawCube(MATRIX4 T)
{
	int i, j, id = 0;

	VECTOR4 ph, pt;
	VECTOR3 pih;

	// v�gigmegy�nk a lapokon
	for (i = 0; i < 6; i++)
	{
		// meghat�rozzuk a lap megfelel�i �leit
		VECTOR3 edges[2] =
		{
			vecSub(convertToInhomogen(transformedCube[faces[i].v[0]]),
			convertToInhomogen(transformedCube[faces[i].v[1]])),
			vecSub(convertToInhomogen(transformedCube[faces[i].v[0]]),
				convertToInhomogen(transformedCube[faces[i].v[2]])),
		};

		// kisz�m�tjuk ezekb�l a a lap norm�lvektor�t
		VECTOR3 normal = normalize(crossProduct(edges[0], edges[1]));

		// camer�ba mutat� vektort
		VECTOR3 toCamera;


			toCamera = normalize(vecSub(convertToInhomogen(transformedCube[faces[i].v[0]]), initVector3(0.0f, 0.0f, center)));


		// l�that�s�g eld�nt�se skal�ris szorzat alapj�n
		if (dotProduct(normal, toCamera) > 0) {

			glColor3f(0.5,0.5,0.5);
			glBegin(GL_POLYGON);

			for (int j = 0; j < 4; j++) {
				ph = initVector4(identityCube[faces[i].v[j]].x, identityCube[faces[i].v[j]].y, identityCube[faces[i].v[j]].z, 1.0f);

				pt = mulMatrixVector(T, ph);

				pih = initVector3(pt.x / pt.w, pt.y / pt.w, pt.z / pt.w);

				glVertex2f(pih.x, pih.y);
			}

			glEnd();

			glColor3f(0.0, 1.0, 0.0);
			glBegin(GL_LINE_STRIP);

			for (int j = 0; j < 4; j++) {
				ph = initVector4(identityCube[faces[i].v[j]].x, identityCube[faces[i].v[j]].y, identityCube[faces[i].v[j]].z, 1.0f);

				pt = mulMatrixVector(T, ph);

				pih = initVector3(pt.x / pt.w, pt.y / pt.w, pt.z / pt.w);

				glVertex2f(pih.x, pih.y);
			}

			glEnd();
		}
	}
}

void draw()
{
	glClear(GL_COLOR_BUFFER_BIT);

	glPointSize(10.0);

	double now = glfwGetTime();
	if (now - lastUpdate >= updateFrequency) {
		forog += 3.14f / 760.0f;;
		initTransformations();
		lastUpdate = now;
	}

	drawCube(T);


	glFlush();
}


void keyPressed(GLFWwindow * windows, GLint key, GLint scanCode, GLint action, GLint mods) {
	if (action == GLFW_PRESS || GLFW_REPEAT) {
		switch (key) {
		case GLFW_KEY_LEFT:
			alpha += deltaAlpha;
			eye.x = R*cos(alpha);
			eye.z = R*sin(alpha);
			initViewMatrix(view, eye, centerVec, up);
			initTransformations();
			break;
		case GLFW_KEY_RIGHT:
			alpha -= deltaAlpha;
			eye.x = R*cos(alpha);
			eye.z = R*sin(alpha);
			initViewMatrix(view, eye, centerVec, up);
			initTransformations();
			break;


		case GLFW_KEY_W:
			alphaFel += PI / 5;
			eye.y = alphaFel;
			initViewMatrix(view, eye, centerVec, up);
			initTransformations();
			break;
		case GLFW_KEY_S:
			alphaFel -= PI / 5;
			eye.y = alphaFel;
			initViewMatrix(view, eye, centerVec, up);
			initTransformations();
			break;


		case GLFW_KEY_UP:
			R -= 0.4f;
			eye.x = R * cos(alpha);
			eye.z = R * sin(alpha);
			initViewMatrix(view, eye, centerVec, up);
			initTransformations();
			break;
		case GLFW_KEY_DOWN:
			R += 0.4f;
			eye.x = R* cos(alpha);
			eye.z = R * sin(alpha);
			initViewMatrix(view, eye, centerVec, up);
			initTransformations();
			break;
		}
	}

	glfwPollEvents();
}

int main(int argc, char ** argv)
{
	GLFWwindow* window;

	/* Initialize the library */
	if (!glfwInit())
		return -1;

	/* Create a windowed mode window and its OpenGL context */
	window = glfwCreateWindow(winWidth, winHeight, "Wonderful cube", NULL, NULL);
	if (!window)
	{
		glfwTerminate();
		return -1;
	}

	/* Make the window's context current */
	glfwMakeContextCurrent(window);

	glfwSetKeyCallback(window, keyPressed);

	init();
	lastUpdate = glfwGetTime();
	/* Loop until the user closes the window */
	while (!glfwWindowShouldClose(window))
	{
		draw();

		glfwSwapBuffers(window);

		glfwPollEvents();
	}

	glfwTerminate();

	return 0;
}

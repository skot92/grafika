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
* 4x4 mátrix, sorfolytonosan tárolva
*/
typedef GLdouble MATRIX4[4][4];

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

/**
* feltölti az A mátrixot az egységmátrixszal
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

/**
* feltölti az A mátrixot a (0,0,s) középpontú centrális vetítés mátrixszával
*/
void initPersProjMatrix(MATRIX4 A, GLdouble s)
{
	initIdentityMatrix(A);

	A[2][2] = 0.0f;
	A[3][2] = -1.0f / s;
}

/**
* feltölti az A mátrixot a (wx,wy):(wx+ww,wy+wh) ablakból (vx,vy):(vx+vw,vy+vh) nézetbe
* történõ transzformáció mátrixszával
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
	// center az a pont,amely fele a kamerat tartjuk, eye a kamera helyét adja
	VECTOR3 centerMinusEye = initVector3(center.x - eye.x, center.y - eye.y, center.z - eye.z);

	// a fenti vektor -1 szeresének egyseg hosszra normáltja, ebbol lesz a kamera rendszerének z-tengelye
	VECTOR3 f = normalizeVector(initVector3(-centerMinusEye.x, -centerMinusEye.y, -centerMinusEye.z));

	// az up vektor es a leendo z tengelyirany vektorialis szorzata adja a kamera x-tengelyiranyat
	VECTOR3 s = normalizeVector(crossProduct(up, f));

	// a kamera y tengelyiranya a mar elkeszult z irany es x irany vektorialis szorzataként jön ki
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

void initRotationMatrixX(MATRIX4 A, GLdouble alpha)
{
	GLdouble c = cos(alpha);
	GLdouble s = sin(alpha);

	initIdentityMatrix(A);

	A[1][1] = c;
	A[1][2] = -s;
	A[2][1] = s;
	A[2][2] = c;
}

void initRotationMatrixY(MATRIX4 A, GLdouble alpha)
{
	GLdouble c = cos(alpha);
	GLdouble s = sin(alpha);

	initIdentityMatrix(A);

	A[0][0] = c;
	A[0][2] = s;
	A[2][0] = -s;
	A[2][2] = c;
}

void initRotationMatrixZ(MATRIX4 A, GLdouble alpha)
{
	GLdouble c = cos(alpha);
	GLdouble s = sin(alpha);

	initIdentityMatrix(A);

	A[0][0] = c;
	A[0][1] = -s;
	A[1][0] = s;
	A[1][1] = c;
}

void initEltolasMatrix(MATRIX4 A, VECTOR3 P) {
	initIdentityMatrix(A);
	A[0][3] = P.x;
	A[1][3] = P.y;
	A[2][3] = P.z;
}

void initScaleMatrix(MATRIX4 A, VECTOR3 P) {
	initIdentityMatrix(A);
	A[0][0] = P.x;
	A[1][1] = P.y;
	A[2][2] = P.z;
}


/**
* visszaadja az A mátrix és a v vektor szorzatát, A*v-t
*/
VECTOR4 mulMatrixVector(MATRIX4 A, VECTOR4 v) {
	return initVector4(
		A[0][0] * v.x + A[0][1] * v.y + A[0][2] * v.z + A[0][3] * v.w,
		A[1][0] * v.x + A[1][1] * v.y + A[1][2] * v.z + A[1][3] * v.w,
		A[2][0] * v.x + A[2][1] * v.y + A[2][2] * v.z + A[2][3] * v.w,
		A[3][0] * v.x + A[3][1] * v.y + A[3][2] * v.z + A[3][3] * v.w);
}

/**
* feltölti a C mátrixot az A és B mátrixok szorzatával, A*B-vel
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
* ablak mérete
*/
GLdouble winWidth = 800.0f, winHeight = 600.0f;


MATRIX4 view;

/**
* merõleges és centrális vetítések mátrixai
*/
MATRIX4 Vc;

/**
* Wtv mátrixok
*/
MATRIX4 Wtv;

/**
* a fenti mátrixokból elõállított két transzformációs mátrix
*/
MATRIX4 TcTorusX, TcTorusY, TcTorusZ, TcGrid,TcCube1, TcCube2, TcCube3, TcCube4;



GLdouble alpha = PI / 2, alphaFel = 0.0f, deltaAlpha = PI / 80.0f;
GLdouble forog = 0.0f;


/**
* centrális vetítés középpontjának Z koordinátája
*/
GLdouble center = 6.0f;

/**
* egységkocka csúcsai
*/
VECTOR4 identityCube[8] =
{

	
	//c5
	initVector4(0.5, 0.5,-0.5, 1.0),
	//c6
	initVector4(-0.5, 0.5,-0.5, 1.0),
	
	//c2
	initVector4(-0.5,-0.5,-0.5, 1.0),
	//c1
	initVector4(0.5,-0.5,-0.5, 1.0),

	//c4
	initVector4(0.5, 0.5, 0.5, 1.0),
	//c7
	initVector4(-0.5, 0.5, 0.5, 1.0),

	//c3
	initVector4(-0.5,-0.5, 0.5, 1.0),
	//c0
	initVector4(0.5,-0.5, 0.5, 1.0),
};

/**
* egységkocka lapjainak indexei
*/
GLuint faces[24] =
{
	0,1,2,3, //alsó
	0,1,5,4, //jobb
	1,2,6,5, //hátsó
	2,3,7,6, //bal
	3,0,4,7, //elsõ
	4,5,6,7, //felsõ
};


/**
* nézet koordinátái
*/
GLdouble cX = (winWidth - winHeight) / 2.0f, cY = 0.0f, cW = winHeight, cH = winHeight;
/**
* elõállítja a szükséges mátrixokat
*/
void initTransformations()
{

	initPersProjMatrix(Vc, center);
	initWtvMatrix(Wtv, -4.0f, -4.0f, 8.0f, 8.0f, cX, cY, cW, cH);

	MATRIX4 rX, rY, rZ, el,Tmp, tmp1, tmp2,s;
	initRotationMatrixX(rX, forog);
	initRotationMatrixY(rY, forog );
	initRotationMatrixZ(rZ, forog );
	initEltolasMatrix(el, initVector3(0, 9, 0));


	//torus
	mulMatrices(el, rX, Tmp);
	mulMatrices(view, Tmp, tmp1);
	mulMatrices(Vc, tmp1, tmp2);
	mulMatrices(Wtv, tmp2, TcTorusX);

	mulMatrices(el, rY, Tmp);
	mulMatrices(view, Tmp, tmp1);
	mulMatrices(Vc, tmp1, tmp2);
	mulMatrices(Wtv, tmp2, TcTorusY);


	mulMatrices(el, rZ, Tmp);
	mulMatrices(view, Tmp, tmp1);
	mulMatrices(Vc, tmp1, tmp2);
	mulMatrices(Wtv, tmp2, TcTorusZ);

	//racs
	mulMatrices(Vc, view, Tmp);
	mulMatrices(Wtv, Tmp, TcGrid);

	//hasab
	initScaleMatrix(s, initVector3(5, 18, 5));
	initEltolasMatrix(el, initVector3(9.5, 9, 9.5));
	mulMatrices(el, s, Tmp);
	mulMatrices(view, Tmp, tmp1);
	mulMatrices(Vc, tmp1, tmp2);
	mulMatrices(Wtv, tmp2, TcCube1);

	initEltolasMatrix(el, initVector3(-9.5, 9, -9.5));
	mulMatrices(el, s, Tmp);
	mulMatrices(view, Tmp, tmp1);
	mulMatrices(Vc, tmp1, tmp2);
	mulMatrices(Wtv, tmp2, TcCube2);

	initEltolasMatrix(el, initVector3(9.5, 9, -9.5));
	mulMatrices(el, s, Tmp);
	mulMatrices(view, Tmp, tmp1);
	mulMatrices(Vc, tmp1, tmp2);
	mulMatrices(Wtv, tmp2, TcCube3);

	initEltolasMatrix(el, initVector3(-9.5, 9, 9.5));
	mulMatrices(el, s, Tmp);
	mulMatrices(view, Tmp, tmp1);
	mulMatrices(Vc, tmp1, tmp2);
	mulMatrices(Wtv, tmp2, TcCube4);

}

/*======================================*/

VECTOR3 eye, up, centerVec;

double R = 37;

void init()
{
	glClearColor(1.0f, 1.0f, 1.0f, 0.0f);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0.0f, winWidth, 0.0f, winHeight, 0.0f, 1.0f);

	eye = initVector3(R*cos(alpha), alphaFel, R*sin(alpha)); //megadja a kamera pozícióját (Ez legyen most a z tengely pozitív felén)
	centerVec = initVector3(0.0, 0.0, 0.0); //megadja, hogy merre néz a kamera (Ez legyen most az origó)
	up = initVector3(0.0, 1.0, 0.0); //megdja, hogy merre van a felfele irány (Ez legyen most az y tengely)

	initViewMatrix(view, eye, centerVec, up);

	initTransformations();
}


void drawCube(VECTOR3 color, MATRIX4 T)
{
	int i, j, id = 0;
	VECTOR4 ph, pt;
	VECTOR3 pih;

	glLineWidth(2.0f);
	glColor3d(color.x, color.y, color.z);

	for (i = 0;i < 6;i++)
	{
		glBegin(GL_LINE_LOOP);
		for (j = 0;j < 4;j++)
		{
			ph = initVector4(identityCube[faces[id]].x, identityCube[faces[id]].y, identityCube[faces[id]].z, 1.0f);
			pt = mulMatrixVector(T, ph);
			pih = initVector3(pt.x / pt.w, pt.y / pt.w, pt.z / pt.w);

			glVertex2d(pih.x, pih.y);

			id++;
		}
		glEnd();
	}
}


void drawSphere(VECTOR3 color, MATRIX4 T)
{
	GLdouble c = 8;
	GLdouble a = 1;
	VECTOR4 ph, pt;
	VECTOR3 pih;

	glLineWidth(2.0f);
	glColor3d(color.x, color.y, color.z);

	for (double u = 0; u <= 2 * PI + 0.001; u += PI / 6) {
		glBegin(GL_LINE_STRIP);
		for (double v = 0; v <= 2 * PI + 0.001; v += PI / 6) {
			ph = initVector4((c + a * cos(v)) * cos(u), a * sin(v), (c + a * cos(v)) * sin(u), 1.0);
			pt = mulMatrixVector(T, ph);
			pih = initVector3(pt.x / pt.w, pt.y / pt.w, pt.z / pt.w);
			glVertex2d(pih.x, pih.y);
		}
		glEnd();
	}
	for (double v = 0; v <= 2 * PI + 0.001; v += PI / 6) {
		glBegin(GL_LINE_STRIP);
		for (double u = 0; u <= 2 * PI + 0.001; u += PI / 6) {
			ph = initVector4((c + a * cos(v)) * cos(u), a * sin(v), (c + a * cos(v)) * sin(u), 1.0);
			pt = mulMatrixVector(T, ph);
			pih = initVector3(pt.x / pt.w, pt.y / pt.w, pt.z / pt.w);
			glVertex2d(pih.x, pih.y);
		}
		glEnd();
	}
}

void drawGrid(VECTOR3 color, MATRIX4 T) {
	glColor3d(color.x, color.y, color.z);
	glLineWidth(2.0f);

	int  id = 0;
	VECTOR4 ph, pt;
	VECTOR3 pih;

	for (double u = -12; u <= 12; u += 1) {
		glBegin(GL_LINE_STRIP);
		for (double v = -12; v <= 12; v += 1) {
			ph = initVector4(v, 0, u, 1.0);
			pt = mulMatrixVector(T, ph);
			pih = initVector3(pt.x / pt.w, pt.y / pt.w, pt.z / pt.w);
			glVertex2d(pih.x, pih.y);
		}
		glEnd();
	}


	for (double v = -12; v <= 12; v += 1) {
		glBegin(GL_LINE_STRIP);
		for (double u = -12; u <= 12; u += 1) {
			ph = initVector4(v, 0, u, 1.0);
			pt = mulMatrixVector(T, ph);
			pih = initVector3(pt.x / pt.w, pt.y / pt.w, pt.z / pt.w);
			glVertex2d(pih.x, pih.y);
		}
		glEnd();
	}
}


void draw()
{
	glClear(GL_COLOR_BUFFER_BIT);

	glPointSize(10.0);

	double now = glfwGetTime();
	if (now - lastUpdate >= updateFrequency) {
		forog += 3.14f / 180.0f;;
		initTransformations();
		lastUpdate = now;
	}

	drawGrid(initVector3(0.0f, 0.0f, 0.0f), TcGrid);
	//piros
	drawSphere(initVector3(1.0f, 0.0f, 0.0f), TcTorusX);
	//zold
	drawSphere(initVector3(0.0f, 1.0f, 0.0f), TcTorusY);
	//kek
	drawSphere(initVector3(0.0f, 0.0f, 1.0f), TcTorusZ);
	drawCube(initVector3(0.0f, 0.0f, 1.0f), TcCube1);
	drawCube(initVector3(0.0f, 0.0f, 1.0f), TcCube2);
	drawCube(initVector3(0.0f, 0.0f, 1.0f), TcCube3);
	drawCube(initVector3(0.0f, 0.0f, 1.0f), TcCube4);
	

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

#include <GLFW\glfw3.h>		// (or others, depending on the system in use)
#include <math.h>
#include <stdio.h>
#include <iostream>

GLFWwindow* window; // GLFW window

GLdouble updateFrequency = 0.001, lastUpdate;

//inhomogen koordinata
typedef struct vector2 { GLdouble x, y; } VECTOR2;

//homogen koordinata
typedef struct { GLdouble x, y, z; } VECTOR3;

typedef GLdouble MATRIX4[4][4];

typedef GLdouble MATRIX3[3][3];

//kell-e forgatni
GLint forgasVan = 0;

//hany pont van lerakva az indulasnal
GLint pontokszama = 0;

//melyik pont korul forog
GLint forogXPontKorul = -1;

GLsizei winWidth = 800, winHeight = 600;

#define bezierPontokSzama 5
#define hermitePontokSzama 3

VECTOR2 bezierPoints[bezierPontokSzama] = { -100, -300, -200, -100, -300, -200, -400, -100, -400, -400 };

//ez lesz a G, utolsó pont az érintõ , +1, põedig az erinto miatt kell
VECTOR2 hermitePoints[hermitePontokSzama + 1] = { bezierPoints[bezierPontokSzama - 1], -500, -300, -600, -350, -400, -500 };


//G*Mvesszo eredmenye
VECTOR2 C[4];

/*
0.000 1.000 8.000 0.000
0.000 1.000 4.000 0.000
0.000 1.000 2.000 1.000
1.000 1.000 1.000 0.000
*/

MATRIX4 Mvesszo = {
	0.750, -1.750,  0.000,  1.000,
	-1.000,  2.000,  0.000,  0.000,
	0.250,-0.250 , 0.000,  0.000,
	0.500, -1.500 , 1.000,  0.000 };

GLint Round(GLfloat n) { return (GLint)(n + 0.5); }

bool movePoint = false;
GLint selectedPointIndex = -1;
GLint selectedPointIndex2 = -2;


VECTOR2 initPoint2D(GLfloat x, GLfloat y) {
	VECTOR2 P;
	P.x = x;
	P.y = y;
	return P;
}

VECTOR3 initVector3(GLfloat x, GLfloat y, GLfloat z) {
	VECTOR3 P;
	P.x = x;
	P.y = y;
	P.z = z;
	return P;
}


/*
*  Ket pont tavolsaganak negyzetet adja vissza.
*/
GLfloat dist2(VECTOR2 P1, VECTOR2 P2) {
	GLfloat t1 = P1.x - P2.x;
	GLfloat t2 = P1.y - P2.y;
	return t1 * t1 + t2 * t2;
}



//egysegmatrixot csinal
void initIdentityMatrix(MATRIX3 A)
{
	int i, j;

	for (i = 0;i < 3;i++)
		for (j = 0;j < 3;j++)
			A[i][j] = 0.0f;

	for (i = 0;i < 3;i++)
		A[i][i] = 1.0f;
}

//forgatas matrixa alpha szoggel
void initForgatasMatrix(MATRIX3 A, float alpha) {
	float c = cos(alpha);
	float s = sin(alpha);

	initIdentityMatrix(A);

	A[0][0] = c;
	A[0][1] = -s;
	A[1][0] = s;
	A[1][1] = c;
}


//eltolas matrixa, P ponttal
void initEltolasMatrix(MATRIX3 A, VECTOR2 P) {
	initIdentityMatrix(A);
	A[0][2] = P.x;
	A[1][2] = P.y;

}

//3X3 matixok osszeszorzasa
void mulMatrices3x3(MATRIX3 A, MATRIX3 B, MATRIX3 C) {
	int i, j, k;

	float sum;
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++) {
			sum = 0;
			for (k = 0; k < 3; k++)
				sum = sum + A[i][k] * B[k][j];
			C[i][j] = sum;
		}
}



//A matrix a G vektorral valo osszeszorzasa, res-be kerul az eredmeny
void arrayMulMatrix4x4(MATRIX4 A, VECTOR2 G[], VECTOR2 res[]) {
	//-G[0] itt szamoljuk ki a megfelelo erintot, es ezzel kell szorozni
	res[0].y = A[0][0] * G[0].y + A[1][0] * G[1].y + A[2][0] * G[2].y + A[3][0] * (G[3].y - G[0].y);
	res[1].y = A[0][1] * G[0].y + A[1][1] * G[1].y + A[2][1] * G[2].y + A[3][1] * (G[3].y - G[0].y);
	res[2].y = A[0][2] * G[0].y + A[1][2] * G[1].y + A[2][2] * G[2].y + A[3][2] * (G[3].y - G[0].y);
	res[3].y = A[0][3] * G[0].y + A[1][3] * G[1].y + A[2][3] * G[2].y + A[3][3] * (G[3].y - G[0].y);

	res[0].x = A[0][0] * G[0].x + A[1][0] * G[1].x + A[2][0] * G[2].x + A[3][0] * (G[3].x - G[0].x);
	res[1].x = A[0][1] * G[0].x + A[1][1] * G[1].x + A[2][1] * G[2].x + A[3][1] * (G[3].x - G[0].x);
	res[2].x = A[0][2] * G[0].x + A[1][2] * G[1].x + A[2][2] * G[2].x + A[3][2] * (G[3].x - G[0].x);
	res[3].x = A[0][3] * G[0].x + A[1][3] * G[1].x + A[2][3] * G[2].x + A[3][3] * (G[3].x - G[0].x);
}


VECTOR3 mulMatrixVector3x3(MATRIX3 A, VECTOR3 v) {
	return initVector3(
		A[0][0] * v.x + A[0][1] * v.y + A[0][2] * v.z,
		A[1][0] * v.x + A[1][1] * v.y + A[1][2] * v.z,
		A[2][0] * v.x + A[2][1] * v.y + A[2][2] * v.z);
}







void init() {
	glClearColor(1.0, 1.0, 1.0, 0.0);
	glMatrixMode(GL_PROJECTION);		// Set projection parameters.
	glLoadIdentity();                   // replace the current matrix with the identity matrix
	glOrtho(0.f, winWidth, winHeight, 0.f, 0.f, 1.f); //multiply the current matrix with an orthographic matrix
	glShadeModel(GL_FLAT);
	glEnable(GL_POINT_SMOOTH);
	//glEnable(GL_LINE_STIPPLE);
	glPointSize(10.0);
	glLineWidth(5.0);
	glLineStipple(1, 0xFF00);

}



//bezier gorbe pontjai iszamitasa
VECTOR2 getCasteljauPoint(int r, int i, double t) {
	if (r == 0)
		return bezierPoints[i];

	VECTOR2 p1 = getCasteljauPoint(r - 1, i, t);
	VECTOR2 p2 = getCasteljauPoint(r - 1, i + 1, t);
	VECTOR2 res = initPoint2D(0, 0);
	res.x = (1 - t) * p1.x + t * p2.x;
	res.y = (1 - t) * p1.y + t * p2.y;
	return res;
}

//hermite iv kiszamitasa, kirajzolasa
void hermite() {
	VECTOR2 C[4] = { 0,0,0,0,0,0,0,0 };
	GLdouble tveszzo[4];
	VECTOR2 draw = initPoint2D(0, 0);

	glLineWidth(2.0);
	glBegin(GL_LINE_STRIP);
	//G * Mvesszo kiszamitasa
	arrayMulMatrix4x4(Mvesszo, hermitePoints, C);
	for (GLdouble t = 0.0; t <= 2; t += 0.001) {
		tveszzo[0] = pow(t, 3);
		tveszzo[1] = pow(t, 2);
		tveszzo[2] = t;
		tveszzo[3] = 1;

		//C*t
		draw.x = C[0].x * tveszzo[0] + C[1].x * tveszzo[1] + C[2].x * tveszzo[2] + C[3].x * tveszzo[3];
		draw.y = C[0].y * tveszzo[0] + C[1].y * tveszzo[1] + C[2].y * tveszzo[2] + C[3].y * tveszzo[3];
		glColor3f(0.0, 0.0, 1.0);
		glVertex2f(draw.x, draw.y);

	}
	glEnd();
}




//forgatas utani ponto kiszamítasa
VECTOR2 forgatas(VECTOR2 centralPoint, VECTOR2 otherPoint) {

	MATRIX3 elXforgatas, elTolasOrigoba, forgatas, visszatolas, vegso;

	//eltolas az origoba
	VECTOR2 tmp = initPoint2D(0, 0);
	tmp.x = -centralPoint.x;
	tmp.y = -centralPoint.y;

	//eltolas matrix
	initEltolasMatrix(elTolasOrigoba, tmp);
	//forgatas matrix
	initForgatasMatrix(forgatas, 0.01);
	//visszatolas matrix
	initEltolasMatrix(visszatolas, centralPoint);

	mulMatrices3x3(forgatas, elTolasOrigoba, elXforgatas);
	mulMatrices3x3(visszatolas, elXforgatas, vegso);

	VECTOR3 ph, pt;
	vector2 pih;
	//homogen alak
	ph = initVector3(otherPoint.x, otherPoint.y, 1);
	//pont es matrix osszeszorzasa
	pt = mulMatrixVector3x3(vegso, ph);
	//vissza inhomogen alakba
	pih = initPoint2D(pt.x / pt.z, pt.y / pt.z);
	return pih;
}

//hogy forogjanak a pontok 
void forgatasIndul(VECTOR2 adottPontKorulForog) {
	//+1 hogy a erinto pont is forogjon
	for (int i = 0; i < hermitePontokSzama + 1; i++)
		hermitePoints[i] = forgatas(adottPontKorulForog, hermitePoints[i]);
	for (int i = 0; i < bezierPontokSzama; i++)
		bezierPoints[i] = forgatas(adottPontKorulForog, bezierPoints[i]);
}

void display() {
	glClear(GL_COLOR_BUFFER_BIT);

	// ha nics meg osszes pont indulasnal
	if (pontokszama < 7) {
		glColor3f(1.0, 0.0, 0.0);
		glPointSize(5.0);
		glBegin(GL_POINTS);
		for (int i = 0; i < bezierPontokSzama; i++)
			glVertex2f(bezierPoints[i].x, bezierPoints[i].y);
		glEnd();

		glBegin(GL_POINTS);
		for (int i = 0; i < hermitePontokSzama; i++)
			glVertex2f(hermitePoints[i].x, hermitePoints[i].y);
		glEnd();
	}
	else {
		double now = glfwGetTime();
		if (now - lastUpdate >= updateFrequency) {
			//kell-e forogni
			if (forgasVan == 1) {
				//hermite iv 2. pont kurok forog
				if (forogXPontKorul == 5)
					forgatasIndul(hermitePoints[1]);
				//hermite iv 3. pont kurok forog
				else if (forogXPontKorul == 6)
					forgatasIndul(hermitePoints[2]);
				//tobbi pont kkorul forog
				else
					forgatasIndul(bezierPoints[forogXPontKorul]);
			}
			lastUpdate = now;
		}

		GLint i;


		glColor3f(1.0, 0.0, 0.0);
		glPointSize(5.0);
		glBegin(GL_POINTS);
		for (i = 0; i < bezierPontokSzama; i++)
			glVertex2f(bezierPoints[i].x, bezierPoints[i].y);
		glEnd();

		glLineWidth(1.0);
		glColor3f(0.0, 0.0, 0.0);
		glBegin(GL_LINE_STRIP);
		for (i = 0; i < bezierPontokSzama; i++)
			glVertex2f(bezierPoints[i].x, bezierPoints[i].y);
		glEnd();

		//bezier gorbe
		glLineWidth(2.0);
		glBegin(GL_LINE_STRIP);
		for (GLdouble t = 0.0; t <= 1.0; t += 0.001) {
			glColor3f(0.0, 1.0, 0.0);
			glPointSize(1.0);
			VECTOR2 tmp = getCasteljauPoint(bezierPontokSzama - 1, 0, t);
			glVertex2f(tmp.x, tmp.y);

		}
		glEnd();


		//hermite
		glColor3f(1.0, 0.0, 0.0);
		glPointSize(5.0);
		glBegin(GL_POINTS);
		for (i = 0; i < hermitePontokSzama; i++)
			glVertex2f(hermitePoints[i].x, hermitePoints[i].y);
		glEnd();

		glColor3f(1.0, 0.0, 1.0);
		glBegin(GL_POINTS);
		glVertex2f(hermitePoints[hermitePontokSzama].x, hermitePoints[hermitePontokSzama].y);
		glEnd();

		glLineWidth(1.0);
		glBegin(GL_LINE_STRIP);
		glVertex2f(hermitePoints[0].x, hermitePoints[0].y);
		glVertex2f(hermitePoints[hermitePontokSzama].x, hermitePoints[hermitePontokSzama].y);
		glEnd();

		hermite();


	}

	glFlush();

}


/*
select the point the use mouse coordinate
*/
GLint getActivePoint1(VECTOR2 p[], GLint size, GLint sens, GLint x, GLint y) {
	GLint i, s = sens * sens;
	VECTOR2 P = initPoint2D(x, y);

	for (i = 0; i < size; i++)
		if (dist2(p[i], P) < s)
			return i;
	return -1;
}

//hermite pontok megfogasa
GLint getActivePoint2(VECTOR2 p[], GLint size, GLint sens, GLint x, GLint y) {
	GLint i, s = sens * sens;
	VECTOR2 P = initPoint2D(x, y);

	for (i = 0; i < size; i++)
		if (dist2(p[i], P) < s)
			return i;
	return -2;
}


//erointo ujraszamolasa
void lastHermitePoint()
{
	GLdouble dx, dy;
	dx = bezierPoints[bezierPontokSzama - 1].x - bezierPoints[bezierPontokSzama - 2].x;
	dy = bezierPoints[bezierPontokSzama - 1].y - bezierPoints[bezierPontokSzama - 2].y;
	hermitePoints[3].x = bezierPoints[bezierPontokSzama - 1].x + 4 * dx;  //ratio * dx;
	hermitePoints[3].y = bezierPoints[bezierPontokSzama - 1].y + 4 * dy;  //ratio * dy;
}


//erinto mozgatasa
void penultBezierPoint()
{
	GLdouble ratio = 3;
	GLdouble dx, dy;
	dx = hermitePoints[3].x - hermitePoints[0].x;
	dy = hermitePoints[3].y - hermitePoints[0].y;
	bezierPoints[bezierPontokSzama - 2].x = hermitePoints[0].x - 1 * dx / 4; // ratio;
	bezierPoints[bezierPontokSzama - 2].y = hermitePoints[0].y - 1 * dy / 4; // ratio;
}

void kozosPontMozgatasa(GLdouble xMouse, GLdouble yMouse) {
	hermitePoints[selectedPointIndex].x = xMouse; //set the selected point new coordinates
	hermitePoints[selectedPointIndex].y = yMouse;
	bezierPoints[bezierPontokSzama - 1].x = xMouse; //set the selected point new coordinates
	bezierPoints[bezierPontokSzama - 1].y = yMouse;
	//itt kellene valtoznia az erintonek is 
	lastHermitePoint();
}

void mouseMove(GLFWwindow* window, GLdouble xMouse, GLdouble yMouse) {
	if (movePoint) {

		//hermite ív pontját kell mozgatni
		if (selectedPointIndex != -1) {
			//bezier pontjat nem kell mozgatni
			selectedPointIndex2 = -2;
			//hermite elso, bezier utolso pont mozgatasa
			if (selectedPointIndex == 0) {
				kozosPontMozgatasa(xMouse, yMouse);
			}
			//erinto mozgatasa
			else if (selectedPointIndex == hermitePontokSzama) {
				hermitePoints[selectedPointIndex].x = xMouse; //set the selected point new coordinates
				hermitePoints[selectedPointIndex].y = yMouse;

				penultBezierPoint();
			}
			//tobbi pont
			else {
				hermitePoints[selectedPointIndex].x = xMouse; //set the selected point new coordinates
				hermitePoints[selectedPointIndex].y = yMouse;
			}
		}
		else {
			selectedPointIndex = getActivePoint1(hermitePoints, hermitePontokSzama + 1, 8, xMouse, yMouse); //which point was clicked by user 
		}


		//bezier gorbe pontjait mozgatas
		if (selectedPointIndex2 != -2) {
			//hermite pontjait nem kell mozgatni
			selectedPointIndex = -1;
			//hermite elso, bezier utolso pont mozgatasa
			if (selectedPointIndex2 == bezierPontokSzama - 1) {
				kozosPontMozgatasa(xMouse, yMouse);
			}
			else if (selectedPointIndex2 == bezierPontokSzama - 2) {
				bezierPoints[selectedPointIndex2].x = xMouse; //set the selected point new coordinates
				bezierPoints[selectedPointIndex2].y = yMouse;
				lastHermitePoint();
			}
			//tobbi pont
			else {
				bezierPoints[selectedPointIndex2].x = xMouse; //set the selected point new coordinates
				bezierPoints[selectedPointIndex2].y = yMouse;
			}
		}
		else {
			selectedPointIndex2 = getActivePoint2(bezierPoints, bezierPontokSzama, 8, xMouse, yMouse); //which point was clicked by user 
		}
	}
}


//pont lerakasa indulasnal
void mouseButton(GLFWwindow* window, GLint button, GLint action, GLint mods) {
	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
		movePoint = true; //left button was presseb, so we need to do some interact inside the mouseMove function

		GLdouble x = 0, y = 0;
		glfwGetCursorPos(window, &x, &y);
		printf("%lf %lf\n", x, y);
		//bezier pontok
		if (pontokszama < 5) {
			bezierPoints[pontokszama].x = x;
			bezierPoints[pontokszama].y = y;
		}
		else if (pontokszama == 5) {
			//hermite iv elso pontja, bezier gorbe utolso pontja ugyan az lesz
			hermitePoints[0].x = bezierPoints[bezierPontokSzama - 1].x;
			hermitePoints[0].y = bezierPoints[bezierPontokSzama - 1].y;
			//hermite ív 2. pontja
			hermitePoints[1].x = x;
			hermitePoints[1].y = y;
		}
		else if (pontokszama == 6) {
			//hermite ív utolso pontja
			hermitePoints[2].x = x;
			hermitePoints[2].y = y;
			lastHermitePoint();

		}
		if (pontokszama < 8) {
			pontokszama++;
		}

	}
	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
		movePoint = false; //left mouse button was released so we can take a rest
		selectedPointIndex = -1; //no point is selected currently
		selectedPointIndex2 = -2;
	}


}

//melyik pontnak is kell forognia
void melyikPontforogjon(GLint index) {
	if (forgasVan == 0) {
		forgasVan = 1;
		forogXPontKorul = index;
	}
	else {
		forgasVan = 0;
	}
}


void keyPressed(GLFWwindow* windows, GLint key, GLint scanCode, GLint action, GLint mods) {
	if (action == GLFW_PRESS || action == GLFW_REPEAT) {
		switch (key) {
		case GLFW_KEY_1:
			melyikPontforogjon(0);
			break;
		case GLFW_KEY_2:
			melyikPontforogjon(1);
			break;
		case GLFW_KEY_3:
			melyikPontforogjon(2);
			break;
		case GLFW_KEY_4:
			melyikPontforogjon(3);
			break;
		case GLFW_KEY_5:
			melyikPontforogjon(4);
			break;
		case GLFW_KEY_6:
			melyikPontforogjon(5);
			break;
		case GLFW_KEY_7:
			melyikPontforogjon(6);
			break;

		}


	}

	glfwPollEvents();
}

int main(int argc, char** argv)
{

	/* Initialize the library */
	if (!glfwInit())
		return -1;

	/* Create a windowed mode window and its OpenGL context */
	window = glfwCreateWindow(winWidth, winHeight, "Hello World", NULL, NULL);
	if (!window)
	{
		glfwTerminate();
		return -1;
	}

	/* Make the window's context current */
	glfwMakeContextCurrent(window);

	//callbacks for mouse
	glfwSetCursorPosCallback(window, mouseMove);
	glfwSetMouseButtonCallback(window, mouseButton);

	glfwSetKeyCallback(window, keyPressed);

	init();

	/* Loop until the user closes the window */
	while (!glfwWindowShouldClose(window))
	{
		display(); /* Render here */

				   /* Swap front and back buffers */
		glfwSwapBuffers(window);

		/* Poll for and process events */
		glfwPollEvents();
	}

	glfwTerminate();

	return 0;
}

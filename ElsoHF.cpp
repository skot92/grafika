#include <GLFW\glfw3.h>				// (or others, depending on the system in use)
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>

using namespace std;


#define PI 3.141592653589793

#define körökSzama 18

#define sugar 20.0

GLdouble updateFrequency = 0.015, lastUpdate;

GLsizei winWidth = 800.f, winHeight = 600.f;



typedef struct kor2d { GLdouble x, y, vX, vY, tomeg, r, pont; } KOR2D;

typedef struct point2d { GLdouble x, y; } POINT2D;

GLdouble balA, balB, balC;
GLdouble jobA, jobB, jobC;
GLdouble fentA, fentB, fentC;
GLdouble lentA, lentB, lentC;


//kepernyo frissites
double now = glfwGetTime();

KOR2D initKor2D(GLdouble x, GLdouble y) {
	KOR2D P;
	P.x = x;
	P.y = y;
	return P;
}

POINT2D initPoint2D(GLdouble x, GLdouble y) {
	POINT2D P;
	P.x = x;
	P.y = y;
	return P;
}

//körök kezdõpontja
KOR2D circles[körökSzama] = {
	initKor2D(600,500) ,	initKor2D(50,100), initKor2D(450,350), initKor2D(400,300) , initKor2D(50,200) , initKor2D(200,200)
	, initKor2D(200,50) , initKor2D(50,300) , initKor2D(300,300)
	, initKor2D(300,50) , initKor2D(50,400) , initKor2D(400,400)
	, initKor2D(400,50) , initKor2D(50,500) , initKor2D(500,500) ,initKor2D(100,500),initKor2D(300,500),400,500 };


void teglalepegyenese()
{
	//y1-y2
	//x2-x1
	//x1*y2 - x2*y1
	//noemalvektorok 
	//x1,y1		x2, y2
	//(800,0) (800,600)
	//jobb egyenes
	jobA = (0.0 - 600.0);
	jobB = (800.0 - 800.0);
	jobC = (800.0 * 600) - (800.0 * 0.0);

	//x1,y1	 x2, y2
	//(0,0) (0,600)
	//ball egyenes
	balA = (0.0 - 600.0);
	balB = (0.0 - 0.0);
	balC = (0.0 * 600.0) - (0, 0 * 0.0);

	//x1,y1	 x2, y2
	//(0,600) (800,600)
	//also egyenes
	lentA = (600.0 - 600.0);
	lentB = (800.0 - 0.0);
	lentC = (0.0 * 600.0) - (800.0*600.0);

	//x1,y1	 x2, y2
	//(0,0) (800,0)
	//felso egyenes
	fentA = (0.0 - 0.0);
	fentB = (800.0 - 0.0);
	fentC = (0.0 * 0.0) - (800.0 * 0.0);
}

double fRand(double fMin, double fMax){
	double f = (double)rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}

void init()
{
	glClearColor(0.7, 0.7, 0.5, 0.0);
	glLineWidth(2.0);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0.f, winWidth, winHeight, 0.f, 0.f, 1.f);

	teglalepegyenese();

	//beállítjuk a kezdõ sebességvektort, illetve a kezdotomeget
	GLdouble kezdoTomeg = 20;
	for (int i = 0; i < körökSzama; i++) {

		circles[i].vX = fRand(-2, 2);
		circles[i].vY = fRand(-2, 2); 
		circles[i].tomeg = kezdoTomeg;
		circles[i].r = sugar;
		circles[i].pont = 0;
	}
	circles[0].tomeg = kezdoTomeg / 2;
	circles[1].tomeg = kezdoTomeg / 2;
}


//pont egyenes tavolsaga
GLdouble pontEgyenesTavolsaga(GLdouble A, GLdouble B, GLdouble C, KOR2D labda)
{
	GLdouble tavolsagx = fabs(A * labda.x + B * labda.y + C) / sqrt(A*A + B*B);
	return tavolsagx;
}


//Ket pont tavolsaganak negyzetet adja vissza.
GLdouble ketPontTavolsaga(KOR2D P1, KOR2D P2) {
	GLdouble t1 = P1.x - P2.x;
	GLdouble t2 = P1.y - P2.y;
	return t1 * t1 + t2 * t2;
}


//körök kirajzolasa
void korRajz() {

	//1. jatekos
	glBegin(GL_LINE_LOOP);
	for (GLdouble t = 0; t <= 2 * PI; t += 0.01)
	{
		glColor3f(1, 0, 0);
		glVertex2d(circles[0].x + circles[0].r * cos(t), circles[0].y + circles[0].r * sin(t));
	}
	glEnd();

	//2. jatekos
	glBegin(GL_LINE_LOOP);
	for (GLdouble t = 0; t <= 2 * PI; t += 0.01)
	{
		glColor3f(0, 0, 1);
		glVertex2d(circles[1].x + circles[1].r * cos(t), circles[1].y + circles[1].r * sin(t));
	}
	glEnd();

	//fekete labda
	glBegin(GL_LINE_LOOP);
	for (GLdouble t = 0; t <= 2 * PI; t += 0.01)
	{
		glColor3f(0, 0, 0);
		glVertex2d(circles[2].x + circles[2].r * cos(t), circles[2].y + circles[2].r * sin(t));
	}
	glEnd();

	//feher labda
	glBegin(GL_LINE_LOOP);
	for (GLdouble t = 0; t <= 2 * PI; t += 0.01)
	{
		glColor3f(1, 1, 1);
		glVertex2d(circles[3].x + circles[3].r * cos(t), circles[3].y + circles[3].r * sin(t));
	}
	glEnd();

	//tobbi labda
	for (int i = 4; i < körökSzama; i++)
	{
		glBegin(GL_LINE_LOOP);
		for (GLdouble t = 0; t <= 2 * PI; t += 0.01)
		{
			glColor3f(0.5, 0.5, 0.5);
			glVertex2d(circles[i].x + circles[i].r * cos(t), circles[i].y + circles[i].r * sin(t));
		}
		glEnd();
	}
}

void labdaNullaz(int labda) {
	circles[labda].r = 0.0;
	circles[labda].x = -1.0;
	circles[labda].y = -1.0;
	circles[labda].vX = 0.0;
	circles[labda].vY = 0.0;
}


void pontSzamNoveles(int jatekos, int labda) {

	////fekete labda, vege a jateknak,elso jatekos nyer
	if ((labda == 2 && jatekos == 0) || (jatekos == 0 && labda == 2)) {
		glClearColor(0.0, 0.0, 1.0, 0.0);
		now = 10.0;
	}
	////fekete labda, vege a jateknak,masodik jatekos nyer
	if ((labda == 2 && jatekos == 1) || (jatekos == 1 && labda == 2)) {
		glClearColor(1.0, 0.0, 0.0, 0.0);
		now = 10.0;
	}

	//feher labda
	if ((labda == 3 && jatekos == 0) || (jatekos == 0 && labda == 3)) {
		circles[0].pont += 5;
		labdaNullaz(3);
	}
	//feher labda
	if ((labda == 3 && jatekos == 1) || (jatekos == 1 && labda == 3)) {
		circles[1].pont += 5;
		labdaNullaz(3);
	}

	if ((labda > 3 && jatekos == 1) || (jatekos == 1 && labda > 3)) {
		circles[1].pont += 1;
		labdaNullaz(labda);
	}

	if ((labda > 3 && jatekos == 0) || (jatekos == 0 && labda > 3)) {
		circles[0].pont += 1;
		labdaNullaz(labda);
	}

	//elfogyak-e a labdak
	if (circles[0].pont + circles[1].pont == körökSzama + 1) {
		cout<<circles[0].pont << "  " << circles[1].pont << endl;
		if (circles[0].pont > circles[1].pont) {
			glClearColor(1.0, 0.0, 0.0, 0.0);
			now = 10.0;
		}
		else {
			glClearColor(0.0, 0.0, 1.0, 0.0);
			now = 10.0;
		}
	}
}

//A és B az egyenes pontjai
KOR2D pattanasFaltol(KOR2D kör, GLdouble A, GLdouble B) {
	kör.vX *= -1;
	kör.vY *= -1;
	kör.vX = 2 * ((kör.vX * A + kör.vY * B) / (A * A + B * B)) * A - kör.vX;
	kör.vY = 2 * ((kör.vX * A + kör.vY * B) / (A * A + B * B)) * B - kör.vY;
	return kör;
}

void kellePattanniAfalaktol()
{
	KOR2D kor1 = initKor2D(0, 0);
	for (int i = 0; i < körökSzama; i++)
	{
		kor1 = circles[i];
		kor1.x += kor1.vX;
		kor1.y += kor1.vY;
		//jobb faltol valo pattanas
		if (pontEgyenesTavolsaga(jobA, jobB, jobC, kor1) < sugar) {
			circles[i] = pattanasFaltol(circles[i], jobA, jobB);
		}

		//ball faltol valo pattanas
		if (pontEgyenesTavolsaga(balA, balB, balC, kor1) < sugar) {
			circles[i] = pattanasFaltol(circles[i], balA, balB);
		}

		//also faltol valo pattanas
		if (pontEgyenesTavolsaga(lentA, lentB, lentC, kor1) < sugar) {
			circles[i] = pattanasFaltol(circles[i], lentA, lentB);
		}

		//felso faltol valo pattanas
		if (pontEgyenesTavolsaga(fentA, fentB, fentC, kor1) < sugar) {
			circles[i] = pattanasFaltol(circles[i], fentA, fentB);
		}
	}
}

void pattanasEgymastol(KOR2D km, KOR2D ka, int i, int j) {

	GLdouble r = km.tomeg / ka.tomeg;

	//n
	POINT2D w = initPoint2D(0.0, 0.0);
	w.x = km.x - ka.x;
	w.y = km.y - ka.y;

	//u
	POINT2D m1 = initPoint2D(0.0, 0.0);
	m1.x = km.vX - ka.vX;
	m1.y = km.vY - ka.vY;

	//un
	POINT2D m1Vesszo = initPoint2D(0.0, 0.0);
	//m1Vesszo = componentVector(m1, w);
	m1Vesszo.x = (m1.x * w.x + m1.y * w.y) / (w.x * w.x + w.y * w.y) * w.x;
	m1Vesszo.y = (m1.x * w.x + m1.y * w.y) / (w.x * w.x + w.y * w.y) * w.y;

	//ut
	POINT2D u1 = initPoint2D(0.0, 0.0);
	u1.x = m1.x - m1Vesszo.x;
	u1.y = m1.y - m1Vesszo.y;

	//vn
	POINT2D m2 = initPoint2D(0.0, 0.0);
	m2.x = (r - 1) * m1Vesszo.x / (r + 1) + u1.x;
	m2.y = (r - 1) * m1Vesszo.y / (r + 1) + u1.y;

	//wn
	POINT2D a2 = initPoint2D(0.0, 0.0);
	a2.x = (m1Vesszo.x * 2 * r) / (r + 1) ;
	a2.y = (m1Vesszo.y * 2 * r) / (r + 1);

	circles[i].vX = m2.x +ka.vX;
	circles[i].vY = m2.y +ka.vY;
	circles[j].vX = a2.x +ka.vX;
	circles[j].vY = a2.y +ka.vY;
}

void kellEEgymastolPattanni() {
	for (int i = 0; i < körökSzama; i++)
	{
		//azt vizsgaljuk hogy a kovetkezo rajzoloasnal hol lenne a korRajz
		KOR2D kor1 = initKor2D(0, 0);
		kor1 = circles[i];
		kor1.x += kor1.vX;
		kor1.y += kor1.vY;
 		for (int j = i + 1; j < körökSzama; j++)
		{
			KOR2D kor2 = initKor2D(0, 0);
			kor2 = circles[j];
			kor2.x += kor2.vX;
			kor2.y += kor2.vY;
			if (ketPontTavolsaga(kor1, kor2) <= ( 2 * (sugar * sugar + sugar * sugar)) )
			//if (sqrt(pow(kor2.x - kor1.x, 2) + pow(kor2.y - kor1.y, 2)) <= sugar + sugar)
			{
				pattanasEgymastol(circles[i], circles[j], i, j);
				pontSzamNoveles(i,j);
			}
		}
	}
}


void display()
{
	if (now != 10.0)
		now = glfwGetTime();
	if (now - lastUpdate >= updateFrequency) {
		for (int i = 0; i < körökSzama; i++)
		{
			circles[i].x += circles[i].vX;
			circles[i].y += circles[i].vY;
			//if(circles[i].vX > 2 || circles[i].vY > 2)
			//std::cout<< circles[i].vX << "	" << circles[i].vY<< std::endl;
		}
		kellePattanniAfalaktol();
		kellEEgymastolPattanni();
		lastUpdate = now;
		
	}
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3f(0.0, 0.4, 0.2);
	if (now != 10.0)
		korRajz();

	glFlush();
}


void keyPressed(GLFWwindow* windows, GLint key, GLint scanCode, GLint action, GLint mods) {
	if (action == GLFW_PRESS || action == GLFW_REPEAT ) {
		switch (key) {
		case GLFW_KEY_I: 
			if (circles[0].vY < 0)
				circles[0] = pattanasFaltol(circles[0], fentA, fentB);
			else
				circles[0] = pattanasFaltol(circles[0], lentA, lentB);
			break;
		case GLFW_KEY_L: 
			if (circles[0].vX > 0)
				circles[0] = pattanasFaltol(circles[0], jobA, jobB);
			else
				circles[0] = pattanasFaltol(circles[0], balA, balB);
			break; 

		case GLFW_KEY_W: 
			if (circles[1].vY < 0)
				circles[1] = pattanasFaltol(circles[1], fentA, fentB);
			else 
				circles[1] = pattanasFaltol(circles[1], lentA, lentB);
			break;
		case GLFW_KEY_D: 
			if (circles[1].vX > 0)
				circles[1] = pattanasFaltol(circles[1], jobA, jobB);
			else 
				circles[1] = pattanasFaltol(circles[1], balA, balB);
			break;
		}
	}

	glfwPollEvents();
}


int main(int argc, char** argv)
{
	GLFWwindow* window;

	/* Initialize the library */
	if (!glfwInit())
		return -1;

	/* Create a windowed mode window and its OpenGL context */
	window = glfwCreateWindow(winWidth, winHeight, "Elso", NULL, NULL);
	if (!window)
	{
		glfwTerminate();
		return -1;
	}

	/* Make the window's context current */
	glfwMakeContextCurrent(window);

	init();

	glfwSetKeyCallback(window, keyPressed);

	/* Store the current time*/
	lastUpdate = glfwGetTime();

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
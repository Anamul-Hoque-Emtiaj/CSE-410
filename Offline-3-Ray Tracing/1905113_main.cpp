#include <GL/glut.h>
#include <cmath>
#include <bits/stdc++.h>
using namespace std;
#include "1905113_classes.h"
#include "bitmap_image.hpp"

// Global variables
struct Point pos;
struct Point look;
struct Point rightVec;
struct Point up;
struct Point ws;
int recursionLevel;
int imageDimention;
int imageCount;
double rotRate;
double movRate;
vector<Object *> objects;
vector<PointLight *> pointlights;
vector<SpotLight *> spotlights;

double windowWidth = 500, windowHeight = 500;
double viewAngle = 80;

void loadData()
{
    ifstream fin("scene.txt");
    int object_count;
    fin >> recursionLevel >> imageDimention >> object_count;
    while (object_count--)
    {
        string obj_type;
        fin >> obj_type;
        Object *obj;
        if (obj_type == "general")
        {
            obj = new General();
            fin >> *((General *)obj);
        }
        else if (obj_type == "triangle")
        {
            obj = new Triangle();
            fin >> *((Triangle *)obj);
        }
        else if (obj_type == "sphere")
        {
            obj = new Sphere();
            fin >> *((Sphere *)obj);
        }
        else
        {
            cout << obj_type << " is not valid" << endl;
        }
        objects.push_back(obj);
    }

    int pointlight_count;
    fin >> pointlight_count;
    while (pointlight_count--)
    {
        PointLight *pl = new PointLight();
        fin >> *((PointLight *)pl);
        pointlights.push_back(pl);
    }

    int spotlight_count;
    fin >> spotlight_count;
    while (spotlight_count--)
    {
        SpotLight *sl = new SpotLight();
        fin >> *((SpotLight *)sl);
        spotlights.push_back(sl);
    }

    Object *floor;
    floor = new Floor();
    floor->setColor(0.5, 0.5, 0.5);
	floor->setCoEfficients(0.4, 0.2, 0.2, 0.2);
    objects.push_back(floor);
}

void capture()
{
    bitmap_image image(imageDimention, imageDimention);
    image.set_all_channels(0, 0, 0); // Set all pixels to black

    double planeDistance = (windowHeight / 2.0) / tan((pi * viewAngle) / (360.0));

	Point topLeft = pos + (look * planeDistance) + (up * (windowHeight / 2.0)) - (rightVec * (windowWidth / 2.0));

	double du = windowWidth / (imageDimention*1.0);
	double dv = windowHeight / (imageDimention*1.0);

	topLeft = topLeft + (rightVec * du / 2.0) - (up * dv / 2.0);
    

    for (int i = 0; i < imageDimention; i++){
        for (int j = 0; j < imageDimention; j++){
            Point pixel = topLeft + (rightVec * j * du) - (up * i * dv);
            Ray ray(pos, pixel - pos);
            Color color(0, 0, 0);

            double tMin = -1;
			int nearestObjIndex = -1;
			for (int k = 0; k < objects.size(); k++)
			{
				double t = objects[k]->intersect(ray, color, 0);
				if(t>0 && (nearestObjIndex == -1 || t<tMin)){
					tMin = t;
					nearestObjIndex = k;
				}
			}

			if (nearestObjIndex != -1)
			{
				Color color(0, 0, 0); 
				double t = objects[nearestObjIndex]->intersect(ray, color,  1);

				if(color.r < 0) color.r = 0;
				if(color.g < 0) color.g = 0;
				if(color.b < 0) color.b = 0;

                if(color.r > 1) color.r = 1;
				if(color.g > 1) color.g = 1;
				if(color.b > 1) color.b = 1;
                image.set_pixel(j, i, 255*color.r, 255*color.g, 255*color.b);
			}
            
        }
    }

    string filename = "output_1" + to_string(imageCount) + ".bmp";
    image.save_image(filename);
    imageCount++;	
}

void init()
{
    loadData();

    pos = {200,0,10};
    look = {-1 / sqrt(2), -1 / sqrt(2), 0};
    rightVec = {-1 / sqrt(2), 1 / sqrt(2), 0};
    up = {0, 0, 1};
    rotRate = pi / 180.0;
    movRate = 3;
    imageCount = 1;

    //clear the screen
	glClearColor(0,0,0,0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(80,	1,	1,	1000.0);
}

void drawAxes() {
    glLineWidth(3);
    glBegin(GL_LINES);
        glColor3f(2,0,0);   // X - Red
        // X axis
        glVertex3f(-100,0,0);
        glVertex3f(100,0,0);

        glColor3f(0,2,0);   // Y - Green
        // Y axis
        glVertex3f(0,-100,0);
        glVertex3f(0,100,0);

        glColor3f(0,0,2);   // Z - Blue
        // Z axis
        glVertex3f(0,0,-100);
        glVertex3f(0,0,100);
    glEnd();
}

void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    gluLookAt(pos.x, pos.y, pos.z,
              pos.x + look.x, pos.y + look.y, pos.z + look.z,
              up.x, up.y, up.z);

    //drawAxes();

    for (int i = 0; i < objects.size(); i++)
    {
        Object *obj = objects[i];
        obj->draw();
    }

    for (int i = 0; i < pointlights.size(); i++)
    {
        pointlights[i]->draw();
    }

    for (int i = 0; i < spotlights.size(); i++)
    {
        spotlights[i]->draw();
    }
    glutSwapBuffers();
}

void keyboard(unsigned char key, int x, int y)
{
    double s;
    switch (key)
    {

    case '0':
        capture();
        break;
    case '1':
        rotate(up, look, rotRate);
        rotate(up, rightVec, rotRate);
        break;
    case '2':
        rotate(up, look, -rotRate);
        rotate(up, rightVec, -rotRate);
        break;
    case '3':
        rotate(rightVec, look, rotRate);
        rotate(rightVec, up, rotRate);
        break;
    case '4':
        rotate(rightVec, look, -rotRate);
        rotate(rightVec, up, -rotRate);
        break;
    case '5':
        rotate(look, rightVec, rotRate);
        rotate(look, up, rotRate);
        break;
    case '6':
        rotate(look, rightVec, -rotRate);
        rotate(look, up, -rotRate);
        break;
    }
    glutPostRedisplay();
}

void special_keyboard(int key, int x, int y)
{
    switch (key)
    {
    case GLUT_KEY_UP:
        pos = pos + look * movRate;
        break;
    case GLUT_KEY_DOWN:
        pos = pos - look * movRate;
        break;
    case GLUT_KEY_LEFT:
        pos = pos - rightVec * movRate;
        break;
    case GLUT_KEY_RIGHT:
        pos = pos + rightVec * movRate;
        break;
    case GLUT_KEY_PAGE_UP:
        pos = pos + up * movRate;
        break;
    case GLUT_KEY_PAGE_DOWN:
        pos = pos - up * movRate;
        break;
    }
    glutPostRedisplay();
}

void reshapeListener(GLsizei width, GLsizei height)
{ // GLsizei for non-negative integer
    // Compute aspect ratio of the new window
    if (height == 0)
        height = 1; // To prevent divide by 0
    GLfloat aspect = (GLfloat)width / (GLfloat)height;

    // Set the viewport to cover the new window
    glViewport(0, 0, width, height);

    // Set the aspect ratio of the clipping area to match the viewport
    glMatrixMode(GL_PROJECTION); // To operate on the Projection matrix
    glLoadIdentity();            // Reset the projection matrix
    gluPerspective(45.0f, aspect, 0.1f, 100.0f);
}

int main(int argc, char **argv)
{
    glutInit(&argc,argv);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);

	glutCreateWindow("1905113: Ray Tracing");

	init();

	glEnable(GL_DEPTH_TEST);

	glutDisplayFunc(display);

	glutKeyboardFunc(keyboard);
	glutSpecialFunc(special_keyboard); 
    glutMainLoop();

    objects.clear();
    objects.shrink_to_fit();

    pointlights.clear();
    pointlights.shrink_to_fit();

    spotlights.clear();
    spotlights.shrink_to_fit();
    return 0;
}
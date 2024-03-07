#include <GL/glut.h>  
#include <vector>
#include <cmath>
#include <iostream>


struct Point
{
	Point() {}
	double x,y,z;

    Point(double x, double y, double z) : x(x), y(y), z(z) {}
    Point(const Point &p) : x(p.x), y(p.y), z(p.z) {}

	// arithemtic operations
	Point operator +(Point b)  {return Point(x+b.x,y+b.y, z+b.z);}
    Point operator -(Point b)  {return Point(x-b.x,y-b.y, z-b.z);}
	Point operator *(double b)  {return Point(x*b,y*b, z*b);}
	Point operator /(double b)  {return Point(x/b,y/b, z/b);}
};

// Global variables
struct Point pos;   
struct Point look;     
struct Point right;     
struct Point up;    
struct Point ws;
struct Point center;
std::vector<std::vector<struct Point>> vertices;

double rotRate = 0.05;
double movRate = 0.1;
double rotAng = 0;
double scale = 1;
double scaleRate = 0.01;
double distFromCenter;


double dihedralAng = 109.47;
double cylAng = 70.5287794;
double cylMaxR = (1.0 / 3) / sin((cylAng / 2) * (M_PI / 180));
double cylMinR = 0;
double cylMaxDist = (1.0/3) / tan((dihedralAng / 2) * (M_PI / 180)) + (1./3) / tan((cylAng / 2) * (M_PI / 180));
double sphereMaxR = 1.0/sqrt(3);

const GLfloat planeCol[][3] = {
    {0.25, 0.96, 0.98},
    {0.94, 0.25, 0.96}
};

const GLfloat sphereCol[][3] = {
    {1, 0, 0},
    {0, 1, 0},
    {0, 0, 1}
};

const GLfloat cylCol[3] = {1, 1, 0};

// Vector Operation
struct Point crossProduct(struct Point p1, struct Point p2){
    return {p1.y * p2.z - p2.y * p1.z, -p1.x * p2.z + p2.x * p1.z, p1.x * p2.y - p2.x * p1.y};
}
double dotProduct(struct Point p1, struct Point p2) { 
    return p1.x * p2.x + p1.y * p2.y + p1.z*p2.z; 
}

void rotate(struct Point& axis, struct Point& vector, double angle){
    struct Point p = crossProduct(axis, vector);
    vector = p * sin(angle) + vector * cos(angle);
}

double magnitude(struct Point p){
    return sqrt(dotProduct(p, p));
}

Point normalize(struct Point p){
    double mag = magnitude(p);
    return p/mag;
}

// OpenGl

void drawAxes() {
    glLineWidth(3);
    glBegin(GL_LINES);
        glColor3f(2,0,0);   // X - Red
        glVertex3f(0,0,0);
        glVertex3f(2,0,0);

        glColor3f(0,2,0);   // Y - Green
        glVertex3f(0,0,0);
        glVertex3f(0,2,0);

        glColor3f(0,0,2);   // Z - Blue
        glVertex3f(0,0,0);
        glVertex3f(0,0,2);
    glEnd();
}
void drawTriangle(int colInd){
    glBegin(GL_TRIANGLES);{
        glColor3fv(planeCol[colInd]);
        glVertex3f(1,0,0);
        glVertex3f(0,1,0);
        glVertex3f(0,0,1);
    }
    glEnd();
}
void drawPlane(int angle, int colInd){
    struct Point t = center * (1 - scale);
    glPushMatrix();
        glRotatef(angle, 0, 1, 0);
        glTranslatef(t.x, t.y, t.z);
        glScaled(scale,scale,scale);
        drawTriangle(colInd);
    glPopMatrix();
}

void drawOctaHedral(){

    for(int i = 0; i < 4; ++i){
        if(i==0 || i==2){
            drawPlane(i * 90, 0);
        }
        else
            drawPlane(i * 90, 1);
    }
        
    
    glPushMatrix();
        glRotatef(180, 1, 0, 0);
        for(int i = 0; i < 4; ++i){
            if(i==0 || i==2){
                drawPlane(i * 90, 0);
            }
            else
                drawPlane(i * 90, 1);
        }
        
    glPopMatrix();

}

std::vector<std::vector<struct Point>> buildUnitPositiveX(int subdivision)
{
    const float DEG2RAD = acos(-1) / 180.0f;

    std::vector<std::vector<struct Point>> vertices;
    float n1[3];        // normal of longitudinal plane rotating along Y-axis
    float n2[3];        // normal of latitudinal plane rotating along Z-axis
    float v[3];         // direction vector intersecting 2 planes, n1 x n2
    float a1;           // longitudinal angle along Y-axis
    float a2;           // latitudinal angle along Z-axis

    // compute the number of vertices per row, 2^n + 1
    int pointsPerRow = (int)pow(2, subdivision) + 1;

    // rotate latitudinal plane from 45 to -45 degrees along Z-axis (top-to-bottom)
    for(unsigned int i = 0; i < pointsPerRow; ++i)
    {
        // normal for latitudinal plane
        // if latitude angle is 0, then normal vector of latitude plane is n2=(0,1,0)
        // therefore, it is rotating (0,1,0) vector by latitude angle a2
        a2 = DEG2RAD * (45.0f - 90.0f * i / (pointsPerRow - 1));
        n2[0] = -sin(a2);
        n2[1] = cos(a2);
        n2[2] = 0;

        std::vector<struct Point> row;
        // rotate longitudinal plane from -45 to 45 along Y-axis (left-to-right)
        for(unsigned int j = 0; j < pointsPerRow; ++j)
        {
            // normal for longitudinal plane
            // if longitude angle is 0, then normal vector of longitude is n1=(0,0,-1)
            // therefore, it is rotating (0,0,-1) vector by longitude angle a1
            a1 = DEG2RAD * (-45.0f + 90.0f * j / (pointsPerRow - 1));
            n1[0] = -sin(a1);
            n1[1] = 0;
            n1[2] = -cos(a1);

            // find direction vector of intersected line, n1 x n2
            v[0] = n1[1] * n2[2] - n1[2] * n2[1];
            v[1] = n1[2] * n2[0] - n1[0] * n2[2];
            v[2] = n1[0] * n2[1] - n1[1] * n2[0];

            // normalize direction vector
            float scale = 1 / sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
            v[0] *= scale;
            v[1] *= scale;
            v[2] *= scale;

            // add a point into vector
            row.push_back({v[0], v[1], v[2]});
        }
        vertices.push_back(row);
    }

    return vertices;
}

void drawUnitSphere(){
    glBegin(GL_QUADS);
        int n = vertices.size();
        for(int i = 0; i < n-1; i++){
            int m = vertices[i].size();
            for(int j = 0; j < m-1; j++){
                glVertex3f(vertices[i][j].x, vertices[i][j].y,vertices[i][j].z);
                glVertex3f(vertices[i+1][j].x, vertices[i+1][j].y,vertices[i+1][j].z);
                glVertex3f(vertices[i+1][j+1].x, vertices[i+1][j+1].y,vertices[i+1][j+1].z);
                glVertex3f(vertices[i][j+1].x, vertices[i][j+1].y,vertices[i][j+1].z);
            }
        }
    glEnd();
}

void drawSpherePositiveX(int colInd){
    glPushMatrix();
        glColor3fv(sphereCol[colInd]);
        glTranslatef(scale, 0, 0);
        glScalef((1-scale) / sqrt(3) , (1-scale) / sqrt(3), (1-scale) / sqrt(3));
        drawUnitSphere();
    glPopMatrix();
}

void drawSpheres(){
    glPushMatrix();
        drawSpherePositiveX(0);
        glRotatef(90, 0, 1, 0);

        drawSpherePositiveX(1);
        glRotatef(90, 0, 1, 0);

        drawSpherePositiveX(0);
        glRotatef(90, 0, 1, 0);

        drawSpherePositiveX(1);
        glRotatef(90, 0, 1, 0);
        
        glRotatef(90, 0, 0, 1);
        drawSpherePositiveX(2);

        glRotatef(-180, 0, 0, 1);
        drawSpherePositiveX(2);
    glPopMatrix();
}

void drawCylinderPortion(double h, double r, int segments){
    
    double startAngle = (-cylAng / 2) * (M_PI / 180.0);
    double endAngle = (cylAng / 2) * (M_PI / 180.0);

    double delta = (endAngle - startAngle) / segments;
    double x1, z1, x2, z2;
    glBegin(GL_QUADS);
        while(startAngle < endAngle){
            x1 = r * cos(startAngle);
            z1 = r * sin(startAngle);
            x2 = r * cos(startAngle + delta);
            z2 = r * sin(startAngle + delta);

            glVertex3d(x1, -h / 2, z1);
            glVertex3d(x1, h / 2, z1);
            glVertex3d(x2, h / 2, z2);
            glVertex3d(x2, -h / 2, z2);
            startAngle += delta;
        }
    glEnd();
}

void drawCylinder(){
    double r = (1 - scale) * cylMaxR;
    double offset = ((1 - scale) * cylMaxDist) / sqrt(2);
    glPushMatrix();
        glColor3fv(cylCol);
        glTranslatef((0.5 - offset),(0.5 - offset), 0);
        glRotatef(45, 0,0,1);
        glPushMatrix();
            glScalef(1, (scale) * sqrt(2),1);
            glScalef(r, 1, r);
            drawCylinderPortion(1,1,100);
        glPopMatrix();
    glPopMatrix();
}


void drawSides(){
    glPushMatrix();

        drawCylinder();
        glRotatef(90, 0, 1, 0);
        drawCylinder();
        glRotatef(90, 0, 1, 0);
        drawCylinder();
        glRotatef(90, 0, 1, 0);
        drawCylinder();
        glRotatef(90, 0, 1, 0);

        glRotatef(90, 1, 0, 0);

        drawCylinder();
        glRotatef(90, 0, 0, 1);
        drawCylinder();
        glRotatef(90, 0, 0, 1);
        drawCylinder();
        glRotatef(90, 0, 0, 1);
        drawCylinder();
        glRotatef(90, 0, 0, 1);

        glRotatef(90, 0, 0, 1);
        glRotatef(90, 1, 0, 0);

        drawCylinder();
        glRotatef(90, 0, 1, 0);
        drawCylinder();
        glRotatef(90, 0, 1, 0);
        drawCylinder();
        glRotatef(90, 0, 1, 0);
        drawCylinder();
        glRotatef(90, 0, 1, 0);

    glPopMatrix();
}

void init() {
    // Set "clearing" or background color
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);   // Black and opaque
    glEnable(GL_DEPTH_TEST);   // Enable depth testing for z-culling

    pos = {2.56357, 1.91455, 2.92836};
    look = {-0.615372 ,-0.337015 ,-0.712558};
    right = {0.759506 ,-0.0116736, -0.650395};
    up = {-0.210875 ,0.941427 ,-0.263148};
    ws = {0, 0.2, 0};
    center = {1./3, 1./3, 1./3};
    distFromCenter = magnitude(pos);
    vertices = buildUnitPositiveX(5);
}

void display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    gluLookAt(pos.x, pos.y, pos.z,
              pos.x+look.x*distFromCenter, pos.y+look.y*distFromCenter, pos.z+look.z*distFromCenter,
              up.x, up.y, up.z);
    
    glPushMatrix();
        //drawAxes();
        glRotatef(rotAng, 0, 1, 0);
        drawOctaHedral();
        drawSides();
        drawSpheres();
    glPopMatrix();

    glutSwapBuffers();
}
void keyboard(unsigned char key, int x, int y) {
    switch (key) {
        case ' ' :
            printf("pos: %f %f %f\n",pos.x,pos.y,pos.z);
            printf("look: %f %f %f\n",look.x,look.y,look.z);
            printf("right: %f %f %f\n",right.x,right.y,right.z);
            printf("up: %f %f %f\n",up.x,up.y,up.z);
            printf("angle: %f\n",rotAng);
            printf("scale: %f\n",scale);
            break;
        case '.':
            scale = fmin(scale + scaleRate, 1);
            break;
        case ',':
            scale = fmax(0, scale-scaleRate);
            break;
        case 'a':
            rotAng+=3;
            break;
        case 'd':
            rotAng-=3;
            break;
        case '1':
            rotate(up, look, rotRate);
            rotate(up, right, rotRate);
            break;
        case '2':
            rotate(up, look, -rotRate);
            rotate(up, right, -rotRate);
            break;
        case '3':
            rotate(right, look, rotRate);
            rotate(right, up, rotRate);
            break;
        case '4':
            rotate(right, look, -rotRate);
            rotate(right, up, -rotRate);
            break;
        case '5':
            rotate(look, right, rotRate);
            rotate(look, up, rotRate);
            break;
        case '6':
            rotate(look, right, -rotRate);
            rotate(look, up, -rotRate);
            break;
        
        case 'w':
    
            look = look * distFromCenter - ws;
            distFromCenter = magnitude(look);
            look = normalize(look);
            up = crossProduct(right, look);
            pos = pos + ws;
            break;
        
        case 's':
            look = look * distFromCenter + ws;
            distFromCenter = magnitude(look);
            look = normalize(look);
            up = crossProduct(right, look);
            pos = pos - ws;
            break;
    }
    glutPostRedisplay();
}

void special_keyboard(int key, int x, int y) {
    switch (key) {
        case GLUT_KEY_UP:
            pos = pos + look*movRate;
            break;
        case GLUT_KEY_DOWN:
            pos = pos - look*movRate;
            break;
        case GLUT_KEY_LEFT:
            pos = pos - right*movRate;
            break;
        case GLUT_KEY_RIGHT:
            pos = pos + right*movRate;
            break;
        case GLUT_KEY_PAGE_UP:
            pos = pos + up*movRate;
            break;
        case GLUT_KEY_PAGE_DOWN:
            pos = pos - up*movRate;
            break;
    }
    glutPostRedisplay();
}


void reshapeListener(GLsizei width, GLsizei height) {  // GLsizei for non-negative integer
    // Compute aspect ratio of the new window
    if (height == 0) height = 1;                // To prevent divide by 0
    GLfloat aspect = (GLfloat)width / (GLfloat)height;

    // Set the viewport to cover the new window
    glViewport(0, 0, width, height);

    // Set the aspect ratio of the clipping area to match the viewport
    glMatrixMode(GL_PROJECTION);  // To operate on the Projection matrix
    glLoadIdentity();             // Reset the projection matrix
    gluPerspective(45.0f, aspect, 0.1f, 100.0f);
}


int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutCreateWindow("Magic Cube");
    glutInitWindowSize(840, 840);
    init();
    glutDisplayFunc(display);
    glutReshapeFunc(reshapeListener); 
    glutSpecialFunc(special_keyboard);
    glutKeyboardFunc(keyboard);
    glutMainLoop();
    return 0;
}



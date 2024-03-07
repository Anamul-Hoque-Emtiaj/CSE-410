#include <GL/glut.h>  
#include <vector>
#include <cmath>
#include <queue>
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

struct Event 
{
    double time;
    bool isX;

    // Comparison operator for priority queue
    bool operator>(const Event& other) const {
        return time > other.time;
    }
};

// Global variables
struct Point pos;   
struct Point look;     
struct Point right;     
struct Point up;    
struct Point ws;
struct Point ballPosition;
struct Point ballDirection;
struct Point zAxes;


double rotRate = 0.3;
double movRate = 0.1;
double ballRadius = 0.5;
double ballSpeed = 0.1;
double boundary_lenght = 10;
const int sectorCount = 20;
const int stackCount = 20;
bool manualControl = true;
double simulationTime = 0;
double rollingAngle = 0;
int rollingDir = 0;
double distFromCenter;
std::vector<float> vertices;
std::vector<float> texCoords;
std::priority_queue<Event, std::vector<Event>, std::greater<Event>> eventQueue;


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

double getAngleInXYPlane(const Point& direction) {

    double angle = atan2(direction.y, direction.x);
    if (angle < 0) {
        angle += 2.0 * M_PI;
    }

    return angle*180.0/M_PI;
}

// OpenGl
void drawAxes() {
    glLineWidth(3);
    glBegin(GL_LINES);
        glColor3f(2,0,0);   // X - Red
        // X axis
        glVertex3f(0,0,0);
        glVertex3f(2,0,0);

        glColor3f(0,2,0);   // Y - Green
        // Y axis
        glVertex3f(0,0,0);
        glVertex3f(0,2,0);

        glColor3f(0,0,2);   // Z - Blue
        // Z axis
        glVertex3f(0,0,0);
        glVertex3f(0,0,2);
    glEnd();
}

void drawWall(){
    glBegin(GL_QUADS);{
        glColor3f(1,0,0);
        glVertex3f(boundary_lenght,boundary_lenght,0);
        glVertex3f(boundary_lenght,boundary_lenght,1);
        glVertex3f(-boundary_lenght,boundary_lenght,1);
        glVertex3f(-boundary_lenght,boundary_lenght,0);
    }
    glEnd();
}

void drawSide(){
    glBegin(GL_QUADS);{
        glColor3f(1,0,0);
        glVertex3f(ballRadius*2,0.025,0.025);
        glVertex3f(ballRadius*2,-0.025,0.025);
        glVertex3f(0,-0.025,0.025);
        glVertex3f(0,0.025,0.025);
    }
    glEnd();
}



void drawCheckerboard() {

    int numSquares = 200;
    float squareSize = 2.0;

    for (int i = -numSquares; i < numSquares; ++i) {
        for (int j = -numSquares; j < numSquares; ++j) {
            glBegin(GL_QUADS);
            float color = ((i + j) % 2 == 0) ? 0.8 : 0.2; // Alternating colors
            glColor3f(color, color, color);
            glVertex3f(j * squareSize, i * squareSize, 0.0);
            glVertex3f((j + 1) * squareSize, i * squareSize, 0.0);
            glVertex3f((j + 1) * squareSize, (i + 1) * squareSize, 0.0);
            glVertex3f(j * squareSize, (i + 1) * squareSize, 0.0);
            glEnd();
        }
    }

    //wall
    for(int i=0; i<4; i++){
        glPushMatrix();
            glRotatef(90*i,0,0,1);
            drawWall();
        glPopMatrix();
    }
}

void generateBall() {
    // clear memory of previous arrays
    vertices.clear();
    texCoords.clear();

    float x, y, z, xy;
    float s, t;
    float sectorStep = 2 * M_PI / sectorCount;
    float stackStep = M_PI / stackCount;
    float sectorAngle, stackAngle;

    for (int i = 0; i <= stackCount; ++i) {
        stackAngle = M_PI / 2 - i * stackStep;
        xy = ballRadius * cosf(stackAngle);
        z = ballRadius * sinf(stackAngle);

        for (int j = 0; j <= sectorCount; ++j) {
            sectorAngle = j * sectorStep;

            x = xy * cosf(sectorAngle);
            y = xy * sinf(sectorAngle);

            vertices.push_back(x);
            vertices.push_back(y);
            vertices.push_back(z);

            s = (float)j / sectorCount;
            t = (float)i / stackCount;
            texCoords.push_back(s);
            texCoords.push_back(t);
        }
    }
}

void drawPyramid() {
    glBegin(GL_QUADS);
    // Base
    glColor3f(1.0, 0.0, 0.0); // Red
    glVertex3f(ballRadius*2, 0.05, 0.05);
    glVertex3f(ballRadius*2, -0.05, 0.05);
    glVertex3f(ballRadius*2, -0.05, -0.05);
    glVertex3f(ballRadius*2, 0.05, -0.05);
    glEnd();

    
    glBegin(GL_TRIANGLES);
    // Front face
    glColor3f(0.0, 1.0, 0.0); // Green
    glVertex3f(ballRadius*2.2, 0.0, 0.0);
    glVertex3f(ballRadius*2, 0.05, 0.05);
    glVertex3f(ballRadius*2, -0.05, 0.05);
    

    // Right face
    glColor3f(0.0, 0.0, 1.0); // Blue
    glVertex3f(ballRadius*2.2, 0.0, 0.0);
    glVertex3f(ballRadius*2, 0.05, 0.05);
    glVertex3f(ballRadius*2, 0.05, -0.05);
    
    // Back face
    glColor3f(1.0, 1.0, 0.0); // Yellow
    glVertex3f(ballRadius*2.2, 0.0, 0.0);
    glVertex3f(ballRadius*2, -0.05, -0.05);
    glVertex3f(ballRadius*2, 0.05, -0.05);
    

    // Left face
    glColor3f(1.0, 0.0, 1.0); // Magenta
    glVertex3f(ballRadius*2.2, 0.0, 0.0);
    glVertex3f(ballRadius*2, -0.05, 0.05);
    glVertex3f(ballRadius*2, -0.05, -0.05);
    
    glEnd();
}

void drawArrow() {
    double angle = getAngleInXYPlane(ballDirection);
    glPushMatrix();
    glRotatef(angle,0,0,1);
    drawPyramid();
    for(int i=0; i<4; i++){
        glRotatef(90*i,1,0,0);
        drawSide();
    }
    glPopMatrix();
}



void drawBall() {

    
    struct Point normal;
    normal = crossProduct(ballDirection,zAxes);
    
    // generate sphere vertices
    generateBall();
    // draw the sphere using quads
    glPushMatrix();
    glTranslatef(ballPosition.x, ballPosition.y, ballPosition.z);
    drawArrow();
    glRotatef(rollingAngle,normal.x,normal.y,normal.z);
    glBegin(GL_QUADS);

    int k1, k2;
    
    for (int i = 0; i < stackCount; ++i) {
        k1 = i * (sectorCount + 1);
        k2 = k1 + sectorCount + 1;

        for (int j = 0; j < sectorCount; ++j, ++k1, ++k2) {
            if ((i+j)%2==0) {
                glColor3f(1.0, 1.0, 0.0);  // yellow
            } else {
                glColor3f(0.0, 1.0, 1.0);  // cyan
            }

            // quad
            glVertex3f(vertices[k1 * 3], vertices[k1 * 3 + 1], vertices[k1 * 3 + 2]);
            glVertex3f(vertices[k2 * 3], vertices[k2 * 3 + 1], vertices[k2 * 3 + 2]);
            glVertex3f(vertices[(k2 + 1) * 3], vertices[(k2 + 1) * 3 + 1], vertices[(k2 + 1) * 3 + 2]);
            glVertex3f(vertices[(k1 + 1) * 3], vertices[(k1 + 1) * 3 + 1], vertices[(k1 + 1) * 3 + 2]);

        }
    }

    glEnd();
    
    glPopMatrix();
}
void init() {
    // Set "clearing" or background color
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);   // Black and opaque
    glEnable(GL_DEPTH_TEST);   // Enable depth testing for z-culling

    pos = {0, -5, 2.5};
    look = {0 ,1 ,-0.1};
    right = {1 ,0, 0};
    up = {0 ,0 ,1};
    ws = {0, 0.2, 0};
    distFromCenter = magnitude(pos);
    ballPosition = {0,0,ballRadius};
    ballDirection = {0.5,0.5,0};
    zAxes = {0,0,1};
}

void display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    gluLookAt(pos.x, pos.y, pos.z,
              pos.x+look.x*distFromCenter, pos.y+look.y*distFromCenter, pos.z+look.z*distFromCenter,
              up.x, up.y, up.z);


    drawCheckerboard();
    drawBall();
    //drawAxes();
    glutSwapBuffers();
}

void handleCollision(){
    if(!eventQueue.empty() && eventQueue.top().isX){
        ballDirection.x *=-1;
    }
    else if(!eventQueue.empty()){
        ballDirection.y *=-1;
    }

    double velocity = ballRadius*rotRate/50;
    double dirAngle = getAngleInXYPlane(ballDirection)*M_PI/180;

    double vx = velocity*cos(dirAngle);
    double vy = velocity*sin(dirAngle);

    double dx;
    double dy;

    if(ballDirection.x>0)
        dx = (boundary_lenght - ballPosition.x-ballRadius);
    else
        dx = ballPosition.x+boundary_lenght-ballRadius;


    if(ballDirection.y>0)
        dy = (boundary_lenght - ballPosition.y-ballRadius);
    else
        dy = ballPosition.y+boundary_lenght-ballRadius;

    

    double tx = abs(dx*1.00/vx)+simulationTime;
    double ty = abs(dy*1.00/vy)+simulationTime;

    while (!eventQueue.empty()) {
        eventQueue.pop();
    }

    struct Event e1,e2;
    e1.time = tx;
    e1.isX = true;
    eventQueue.push(e1);

    e2.time = ty;
    e2.isX = false;
    eventQueue.push(e2);
}

void eventDrivenSimulation(int value){

    if(!manualControl){
        simulationTime+=value;

        if (eventQueue.empty() || eventQueue.top().time-simulationTime<=1.0){
            value = 50;
            handleCollision();
        }
        else if(eventQueue.top().time-simulationTime < 50.0){
            value = eventQueue.top().time-simulationTime;
        }
        
        rollingAngle-=(rotRate*180/M_PI);
        if(rollingAngle<0)
            rollingAngle+=360.0;
        
        double s = ballRadius*rotRate*value/50;
        ballPosition = ballPosition + (ballDirection*s);
        glutPostRedisplay();
        
    }
    glutTimerFunc(50, eventDrivenSimulation, value);
}

void timeDrivenSimulation(int value){

    if(!manualControl){
        rollingAngle-=(rotRate*180/M_PI);
        if(rollingAngle<0)
            rollingAngle+=360.0;
        
        double s = ballRadius*rotRate;
        ballPosition = ballPosition + (ballDirection*s);
        if(ballPosition.x+ballRadius>=boundary_lenght || ballPosition.x-ballRadius<=(-boundary_lenght)){
            ballDirection.x *=-1;
        }
        else if(ballPosition.y+ballRadius>=boundary_lenght || ballPosition.y-ballRadius<=(-boundary_lenght)){
            ballDirection.y *=-1;
        }
        glutPostRedisplay();
    }
    glutTimerFunc(100, timeDrivenSimulation, 0);
}

void keyboard(unsigned char key, int x, int y) {
    double s;
    switch (key) {
        
        case ' ':
            manualControl = !manualControl;
            if(!manualControl){
                simulationTime = 0;
                while (!eventQueue.empty()) {
                    eventQueue.pop();
                }
                handleCollision();
                glutTimerFunc(50, eventDrivenSimulation, 50);
            }
            break;
        case 'j':
            rotate(zAxes,ballDirection,rotRate);
            ballDirection = normalize(ballDirection);
            rollingAngle = 0;
            if(!manualControl){
                while (!eventQueue.empty()) {
                    eventQueue.pop();
                }
                handleCollision();
            }
            break;
        case 'l':
            rotate(zAxes,ballDirection,-rotRate);
            ballDirection = normalize(ballDirection);
            rollingAngle = 0;
            if(!manualControl){
                while (!eventQueue.empty()) {
                    eventQueue.pop();
                }
                handleCollision();
            }
            break;
        case 'i':
            if(manualControl){
                rollingAngle-=(rotRate*180/M_PI);
                if(rollingAngle<0)
                    rollingAngle+=360.0;
                
                s = ballRadius*rotRate;
                ballPosition = ballPosition + (ballDirection*s);
                if(ballPosition.x+ballRadius>=boundary_lenght || ballPosition.x-ballRadius<=(-boundary_lenght)){
                    ballDirection.x *=-1;
                }
                else if(ballPosition.y+ballRadius>=boundary_lenght || ballPosition.y-ballRadius<=(-boundary_lenght)){
                    ballDirection.y *=-1;
                }
            }
            break;
        case 'k':
            if(manualControl){
                rollingAngle+=(rotRate*180/M_PI);;
                if(rollingAngle>360)
                    rollingAngle-=360;
                s = ballRadius*rotRate;
                ballPosition = ballPosition + (ballDirection*(-s));
                if(ballPosition.x+ballRadius>=boundary_lenght || ballPosition.x-ballRadius<=(-boundary_lenght)){
                    ballDirection.x *=-1;
                }
                else if(ballPosition.y+ballRadius>=boundary_lenght || ballPosition.y-ballRadius<=(-boundary_lenght)){
                    ballDirection.y *=-1;
                }
            }
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
    glutCreateWindow("Rolling Ball");
    glutInitWindowSize(1200, 1200);
    init();
    glutDisplayFunc(display);
    //glutTimerFunc(50, timeDrivenSimulation, 50);
    glutReshapeFunc(reshapeListener); 
    glutSpecialFunc(special_keyboard);
    glutKeyboardFunc(keyboard);
    glutMainLoop();
    return 0;
}



#include<bits/stdc++.h>
#include "bitmap_image.hpp"

#define PI 2.0*acos(0.0)

using namespace std;

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

    friend istream& operator>>(istream& is, Point& p){
        is >> p.x >> p.y >> p.z;
        return is;
    }
    friend ofstream& operator<<(ofstream& os, Point& p){
        os << fixed << setprecision(7);
        os << p.x << " " << p.y << " " << p.z << endl;
        return os;
    }
};

struct Point crossProduct(struct Point p1, struct Point p2){
    return {p1.y * p2.z - p2.y * p1.z, -p1.x * p2.z + p2.x * p1.z, p1.x * p2.y - p2.x * p1.y};
}
double dot(struct Point p1, struct Point p2) { 
    return p1.x * p2.x + p1.y * p2.y + p1.z*p2.z; 
}
double magnitude(struct Point p){
    return sqrt(dot(p, p));
}

struct Point normalize(struct Point p){
    double mag = magnitude(p);
    return p/mag;
}

struct Point scaleDown(Point p, double n){
    return {p.x/n,p.y/n,p.z/n};
}

static unsigned long int g_seed = 1; 
inline int getrandom()
{
    g_seed = (214013 * g_seed + 2531011); 
    return (g_seed >> 16) & 0x7FFF;
}

struct Triangle{
    struct Point p1,p2,p3;
    int rgb[3];

    void setColor()
    {
        rgb[0] = getrandom();
        rgb[1] = getrandom();
        rgb[2] = getrandom();
    }

    friend istream& operator>>(istream& is, Triangle& t){
        is >> t.p1 >> t.p2 >> t.p3;
        return is;
    }

    friend ofstream& operator<<(ofstream& os, Triangle& t){
        os << t.p1 << t.p2 << t.p3;
        return os;
    }
};

struct Matrix{

    vector<vector<double> > mat;
    int dim;

    Matrix()
    {
        mat.resize(4,vector<double>(4,0));
        dim = 4;
    }

    Matrix(int dim)
    {
        mat.resize(dim,vector<double>(dim,0));
        this->dim = dim;
    }
    Matrix operator*(Matrix other)
    {
        Matrix res(dim);
        for(int i=0; i<dim; i++)
        {
            for(int j=0; j<dim; j++)
            {
                for(int k=0; k<dim; k++)
                {
                    res.mat[i][j] += mat[i][k]*other.mat[k][j];
                }
            }
        }
        return res;
    }

    Point operator*(Point p){

        double res[4]={0.0};
        double p2[4]={p.x, p.y, p.z, 1};

        for(int i=0; i<4; i++){
            for(int j=0; j<4; j++){
                res[i] += mat[i][j] * p2[j];
            }
        }
        Point ret = scaleDown(Point(res[0], res[1], res[2]),res[3]);
        return ret;
    }

    Triangle operator*(Triangle t){
        Triangle ret;
        ret.p1 = (*this)*t.p1;
        ret.p2 = (*this)*t.p2;
        ret.p3 = (*this)*t.p3;
        return ret;
    }
};

struct Matrix getIdentityMat(int dim=4)
{
    struct Matrix I(dim);
    for(int i=0; i<dim; i++)
    {
        for(int j=0; j<dim; j++)
        {
            if(i==j)
                I.mat[i][j] = 1;
            else
                I.mat[i][j] = 0;
        }
    }
    return I;
}

struct Matrix getTranslationMat(double ax, double ay, double az, int dim=4)
{
    struct Matrix res = getIdentityMat(dim);
    res.mat[0][3] = ax;
    res.mat[1][3] = ay;
    res.mat[2][3] = az;
    return res;
}

struct Matrix getScalingMat(double ax, double ay, double az, int dim=4)
{
    struct Matrix res = getIdentityMat(dim);
    res.mat[0][0] = ax;
    res.mat[1][1] = ay;
    res.mat[2][2] = az;
    return res;
}

struct Point applyRodriguesFormula(struct Point x, struct Point a, double theta) 
{   
    double rad = theta*PI/180;
    return x*cos(rad)+(a*dot(a,x))*(1-cos(rad))+crossProduct(a,x)*sin(rad);
}

struct Matrix getRotationMat(Point a, double angle, int dim=4)
{
    a = normalize(a);
    struct Point c1 = applyRodriguesFormula(Point(1.0, 0.0, 0.0), a, angle);
    struct Point c2 = applyRodriguesFormula(Point(0.0, 1.0, 0.0), a, angle);
    struct Point c3 = applyRodriguesFormula(Point(0.0, 0.0, 1.0), a, angle);

    struct Matrix res = getIdentityMat(dim);
    res.mat[0][0] = c1.x;
    res.mat[1][0] = c1.y;
    res.mat[2][0] = c1.z;

    res.mat[0][1] = c2.x;
    res.mat[1][1] = c2.y;
    res.mat[2][1] = c2.z;

    res.mat[0][2] = c3.x;
    res.mat[1][2] = c3.y;
    res.mat[2][2] = c3.z;
    return res;
}

struct Matrix getViewMat(struct Point eye, struct Point look, struct Point up)
{
    struct Point z = (look-eye);
    z = normalize(z);
    struct Point x = crossProduct(z,up);
    x = normalize(x);
    struct Point y = crossProduct(x,z);
    y = normalize(y);

    struct Matrix res = getIdentityMat(4);

    res.mat[0][0] = x.x;
    res.mat[0][1] = x.y;
    res.mat[0][2] = x.z;

    res.mat[1][0] = y.x;
    res.mat[1][1] = y.y;
    res.mat[1][2] = y.z;

    res.mat[2][0] = -z.x;
    res.mat[2][1] = -z.y;
    res.mat[2][2] = -z.z;

    struct Matrix t = getTranslationMat(-eye.x,-eye.y,-eye.z);

    res = t*res;
    return res;
}

struct Matrix getProjectionMat(double fovY, double aspect, double near, double far)
{
    double fovX = fovY*aspect;
    double angleR = (fovX/2)*PI/180.0;
    double r = near * tan(angleR);
    double angleT = (fovY/2)*PI/180.0;
    double t = near * tan(angleT);

    struct Matrix res = getIdentityMat(4);
    
    res.mat[0][0] = near/r;
    res.mat[1][1] = near/t;
    res.mat[2][2] = -(far+near)/(far-near);
    res.mat[3][2] = -1;
    res.mat[2][3] = -2.0*far*near/(far-near);
    res.mat[3][3] = 0;

    return res;
}

int main()
{
    ifstream fin("scene.txt");
    ofstream fout("stage1.txt");
    
    struct Point eye, look, up;
    double fovY, aspectRatio, near, far;
    fin>>eye>>look>>up;
    fin>>fovY>>aspectRatio>>near>>far;

    /*
        stage 1: Modeling Transformation
    */

    stack<struct Matrix> st;
    st.push(getIdentityMat());

    string command;
    while(fin>>command)
    {
        if( command == "triangle" ){
            struct Triangle t;
            fin>>t;
            t = st.top()*t;
            fout<<t;
            fout<<endl;
        } else if( command == "translate" ){
            struct Point p;
            fin>>p;
            struct Matrix m = getTranslationMat(p.x,p.y,p.z);
            st.top() = st.top()*m;
        } else if( command == "scale" ){
            struct Point p;
            fin>>p;
            struct Matrix m = getScalingMat(p.x,p.y,p.z);
            st.top() = st.top()*m;
        } else if( command == "rotate" ){
            struct Point p;
            double theta;
            fin>>theta>>p;
            struct Matrix m = getRotationMat(p,theta);
            st.top() = st.top()*m;
        } else if( command == "push" ){
            st.push(st.top());
        } else if( command == "pop" ){
            if(st.empty())
            {
                cout<<"Error: stack is empty"<<endl;
                return 0;
            }
            st.pop();
        } else if( command == "end" ){
            break;
        }else{
            cout<<"Invalid Command"<<endl;
            break;
        }
    }
    fin.close();
    fout.close();

    /*
        stage 2: View Transformation
    */
    fin.open("stage1.txt");
    fout.open("stage2.txt");

    struct Matrix viewTransform = getViewMat(eye,look,up);
    struct Triangle triangle;
    while(fin>>triangle){
        triangle = viewTransform*triangle;
        fout<<triangle<<endl;
    }

    fin.close();
    fout.close();

    /*
        stage 3: Projection Transformation
    */
    fin.open("stage2.txt");
    fout.open("stage3.txt");

    Matrix projTransform = getProjectionMat(fovY,aspectRatio,near,far);
    while(fin>>triangle){
        triangle = projTransform*triangle;
        fout<<triangle<<endl;
    }

    fin.close();
    fout.close();

    /*
        stage 4: Clipping and Scan Conversion using Z-Buffer
    */

    fin.open("config.txt");
    fout.open("z_buffer.txt");

    int screenWidth, screenHeight;
    fin>>screenWidth>>screenHeight;

    double boxLeft=-1, boxRight=1, boxUp=1, boxDown=-1;
    double zMin=-1, zMax=1;

    double dx = (boxRight - boxLeft)/screenWidth;
    double dy = (boxUp - boxDown)/screenHeight;
    double topY = boxUp - dy/2; 
    double bottomY = boxDown + dy/2;
    double leftX = boxLeft + dx/2;
    double rightX = boxRight - dx/2;

    // initialize z-buffer and frame buffer
    vector<vector<double>> z_buffer(screenHeight,vector<double>(screenWidth,zMax));
    bitmap_image image(screenWidth, screenHeight);
    image.set_all_channels(0, 0, 0);


    fin.close();
    fin.open("stage3.txt");

    while(fin>>triangle){
        triangle.setColor();

        double minX, maxX, minY, maxY;

        minX = min(min(triangle.p1.x,triangle.p2.x),triangle.p3.x);
        maxX = max(max(triangle.p1.x,triangle.p2.x),triangle.p3.x);
        minY = min(min(triangle.p1.y,triangle.p2.y),triangle.p3.y);
        maxY = max(max(triangle.p1.y,triangle.p2.y),triangle.p3.y);

        // clipping
        minX = max(minX,leftX);
        maxX = min(maxX,rightX);

        minY = max(minY,bottomY);
        maxY = min(maxY,topY);
        
        // top_scanline and bottom_scanline
        int bottom_scanline = round((topY-minY)/dy);
        int top_scanline = round((topY-maxY)/dy);

        for(int scanY = top_scanline;scanY <= bottom_scanline; scanY++){

            double y_s = topY - scanY*dy;
            vector<double>x_ab(2),z_ab(2);
            int cnt = 0;

            if(triangle.p1.y != triangle.p2.y && (y_s >= min(triangle.p1.y,triangle.p2.y) && y_s <= max(triangle.p1.y,triangle.p2.y))){
                x_ab[cnt] = triangle.p1.x - (triangle.p1.x - triangle.p2.x)*(triangle.p1.y - y_s)/(triangle.p1.y - triangle.p2.y); //formula
                z_ab[cnt] = triangle.p1.z - (triangle.p1.z - triangle.p2.z)*(triangle.p1.y - y_s)/(triangle.p1.y - triangle.p2.y);
                cnt++;
            }

            if(triangle.p2.y != triangle.p3.y && (y_s >= min(triangle.p2.y,triangle.p3.y) && y_s <= max(triangle.p2.y,triangle.p3.y))){
                x_ab[cnt] = triangle.p2.x - (triangle.p2.x - triangle.p3.x)*(triangle.p2.y - y_s)/(triangle.p2.y - triangle.p3.y);
                z_ab[cnt] = triangle.p2.z - (triangle.p2.z - triangle.p3.z)*(triangle.p2.y - y_s)/(triangle.p2.y - triangle.p3.y);
                cnt++;
            }

            if(triangle.p3.y != triangle.p1.y && (y_s >= min(triangle.p3.y,triangle.p1.y) && y_s <= max(triangle.p3.y,triangle.p1.y))){
                x_ab[cnt] = triangle.p3.x - (triangle.p3.x - triangle.p1.x)*(triangle.p3.y - y_s)/(triangle.p3.y - triangle.p1.y);
                z_ab[cnt] = triangle.p3.z - (triangle.p3.z - triangle.p1.z)*(triangle.p3.y - y_s)/(triangle.p3.y - triangle.p1.y);
                cnt++;
            }

            vector<double>tx_ab(2);
            tx_ab = x_ab;

            // Clipping
            for(int i=0;i<2;i++){
                if(x_ab[i]<minX) x_ab[i] = minX;
                if(x_ab[i]>maxX) x_ab[i] = maxX;
            }

            z_ab[0] = z_ab[1] - (z_ab[1] - z_ab[0])*(tx_ab[1] - x_ab[0])/(tx_ab[1] - tx_ab[0]);
            z_ab[1] = z_ab[1] - (z_ab[1] - z_ab[0])*(tx_ab[1] - x_ab[1])/(tx_ab[1] - tx_ab[0]);


            double xa,za,xb,zb;
            xa = x_ab[0];
            za = z_ab[0];
            xb = x_ab[1];
            zb = z_ab[1];

            if(xa >= xb){
                swap(xa,xb);
                swap(za,zb);
            }

            int left_intersecting_column = round((xa-leftX)/dx);
            int right_intersecting_column = round((xb-leftX)/dx);
            for(int scanX=left_intersecting_column;scanX<=right_intersecting_column;scanX++){
                double xp = leftX + scanX*dx;
                double zp = zb - (zb-za)*((xb-xp)/(xb-xa)); // formula
                if(zp>=zMin){
                    if(zp < z_buffer[scanY][scanX]){
                        z_buffer[scanY][scanX] = zp;
                        image.set_pixel(scanX,scanY,triangle.rgb[0],triangle.rgb[1],triangle.rgb[2]);
                    }
                }
            }

        }

    }

    // save
    image.save_image("out.bmp");
    for (int i = 0; i < screenHeight; i++) 
    {
        for (int j = 0; j < screenWidth; j++)
        {
            if (z_buffer[i][j] < zMax)
            {
                fout << setprecision(6) << fixed << z_buffer[i][j] << "\t";
            }
        }
        fout << endl;
    }
    fin.close();
    fout.close();

    // Free memory
    z_buffer.clear();
    z_buffer.shrink_to_fit();
}
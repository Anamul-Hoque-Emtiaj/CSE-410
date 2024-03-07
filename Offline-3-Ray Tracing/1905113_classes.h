#include <GL/glut.h>
#include <cmath>
#include <bits/stdc++.h>
#define EPSILON 0.00001
#define pi (2 * acos(0.0))

using namespace std;

class Object;
struct PointLight;
struct SpotLight;

extern vector<Object *> objects;
extern vector<PointLight *> pointlights;
extern vector<SpotLight *> spotlights;
extern int recursionLevel;

struct Point
{
	Point() {}
	double x, y, z;

	Point(double x, double y, double z) : x(x), y(y), z(z) {}
	Point(const Point &p) : x(p.x), y(p.y), z(p.z) {}

	// arithemtic operations
	Point operator+(Point b) { return Point(x + b.x, y + b.y, z + b.z); }
	Point operator-(Point b) { return Point(x - b.x, y - b.y, z - b.z); }
	Point operator*(double b) { return Point(x * b, y * b, z * b); }
	Point operator/(double b) { return Point(x / b, y / b, z / b); }

	friend istream &operator>>(istream &is, Point &p)
	{
		is >> p.x >> p.y >> p.z;
		return is;
	}
	friend ofstream &operator<<(ofstream &os, Point &p)
	{
		os << fixed << setprecision(7);
		os << p.x << " " << p.y << " " << p.z << endl;
		return os;
	}
};

// Vector Operation
struct Point crossProduct(struct Point p1, struct Point p2)
{
	return {p1.y * p2.z - p2.y * p1.z, -p1.x * p2.z + p2.x * p1.z, p1.x * p2.y - p2.x * p1.y};
}
double dotProduct(struct Point p1, struct Point p2)
{
	return p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;
}

void rotate(struct Point &axis, struct Point &vector, double angle)
{
	struct Point p = crossProduct(axis, vector);
	vector = p * sin(angle) + vector * cos(angle);
}

double magnitude(struct Point p)
{
	return sqrt(dotProduct(p, p));
}

Point normalize(struct Point p)
{
	double mag = magnitude(p);
	return p / mag;
}

struct Point scaleDown(Point p, double n)
{
	return {p.x / n, p.y / n, p.z / n};
}

static unsigned long int g_seed = 1;
inline int getrandom()
{
	g_seed = (214013 * g_seed + 2531011);
	return (g_seed >> 16) & 0x7FFF;
}

double determinant(double arr[3][3])
{
	return (arr[0][0] * (arr[1][1] * arr[2][2] - arr[1][2] * arr[2][1])) - (arr[0][1] * (arr[1][0] * arr[2][2] - arr[1][2] * arr[2][0])) + (arr[0][2] * (arr[1][0] * arr[2][1] - arr[1][1] * arr[2][0]));
}

struct Color
{
	double r, g, b;
	Color()
	{
		r = g = b = 0.0;
	}
	Color(double r, double g, double b) : r(r), g(g), b(b) {}

	friend istream &operator>>(istream &is, Color &color)
	{
		is >> color.r >> color.g >> color.b;
		return is;
	}

	friend ostream &operator<<(ostream &os, Color &color)
	{
		os << fixed << setprecision(7);
		os << color.r << " " << color.g << " " << color.b << endl;
		return os;
	}
};

// pointlight
struct PointLight
{
	Point pos;
	Color color;

	void draw()
	{
		glPointSize(5);
		glBegin(GL_POINTS);
		glColor3f(color.r, color.g, color.b);
		glVertex3f(pos.x, pos.y, pos.z);
		glEnd();
	}

	// input stream
	friend istream &operator>>(istream &in, PointLight &l)
	{
		in >> l.pos.x >> l.pos.y >> l.pos.z;
		in >> l.color;
		return in;
	}
};

// spotlight
struct SpotLight
{
	Point pos;
	Color color;
	Point dir;
	double cutoffAngle; // this is different from the spotlight

	void draw()
	{
		glPointSize(15);
		glBegin(GL_POINTS);
		glColor3f(color.r, color.g, color.b);
		glVertex3f(pos.x, pos.y, pos.z);
		glEnd();
	}

	// input stream
	friend istream &operator>>(istream &in, SpotLight &l)
	{
		in >> l.pos;
		in >> l.color;
		in >> l.dir;
		in >> l.cutoffAngle;
		return in;
	}
};

struct Ray
{
	Point start, dir;

	Ray(Point start, Point dir)
	{
		this->start = start;
		this->dir = normalize(dir);
	}
};

class Object
{
public:
	Point reference_point;
	double height, width, length;
	Color color;
	double coEfficients[4]; // ambient, diffuse, specular, reflection coefficients
	int shine;				// exponent term of specular component

	Object()
	{
		reference_point = {0, 0, 0};
		height = 0;
		width = 0;
		length = 0;
		color = {0, 0, 0};
		coEfficients[0] = 0;
		coEfficients[1] = 0;
		coEfficients[2] = 0;
		coEfficients[3] = 0;
		shine = 0;
	}
	void setCoEfficients(double a, double b, double c, double d)
	{
		coEfficients[0] = a;
		coEfficients[1] = b;
		coEfficients[2] = c;
		coEfficients[3] = d;
	}
	void setShine(int s)
	{
		shine = s;
	}
	void setColor(double r, double g, double b)
	{
		color = {r, g, b};
	}
	virtual void draw() = 0;
	virtual Color getIntersectPointColor(Point point)
	{
		return Color(this->color.r, this->color.g, this->color.b);
	}
	virtual double levelZeroIntersect(Ray ray, Color &color, int level) = 0;
	virtual Ray getNormalAtIP(Point point, Ray incidentRay) = 0;
	virtual double intersect(Ray ray, Color &color, int level) // intersect with illumination and recursive reflection
	{
		double intersectionDistance = levelZeroIntersect(ray, color, level);

		if (intersectionDistance < 0)
			return -1;
		if (level == 0)
			return intersectionDistance;

		Point intersectPoint = ray.start + ray.dir * intersectionDistance;
		Color intersectPointColor = getIntersectPointColor(intersectPoint);

		// ambient
		color.r = intersectPointColor.r * coEfficients[0];
		color.g = intersectPointColor.g * coEfficients[0];
		color.b = intersectPointColor.b * coEfficients[0];

		// diffuse and specular of PointLight
		for (PointLight *pointlight : pointlights)
		{

			Ray lightRay = Ray(pointlight->pos, intersectPoint - pointlight->pos);
			double lightDistance = magnitude(intersectPoint - pointlight->pos);
			if (lightDistance < EPSILON)
				continue;

			bool obscured = false;
			for (Object *object : objects)
			{
				double tempDistance = object->levelZeroIntersect(lightRay, color, 0);
				if (tempDistance > 0 && tempDistance + EPSILON < lightDistance)
				{
					obscured = true;
					break;
				}
			}

			if (obscured)
				continue;

			// diffuse
			Ray normal = getNormalAtIP(intersectPoint, lightRay);
			double lambertVal = max(0.0, dotProduct(normal.dir, lightRay.dir));
			color.r += intersectPointColor.r * coEfficients[1] * pointlight->color.r * lambertVal;
			color.g += intersectPointColor.g * coEfficients[1] * pointlight->color.g * lambertVal;
			color.b += intersectPointColor.b * coEfficients[1] * pointlight->color.b * lambertVal;

			// specular
			Ray reflectedRay = Ray(intersectPoint, lightRay.dir - normal.dir * 2 * dotProduct(lightRay.dir, normal.dir));
			double specularVal = max(0.0, dotProduct(reflectedRay.dir, ray.dir * -1));
			specularVal = pow(specularVal, shine);
			color.r += intersectPointColor.r * coEfficients[2] * pointlight->color.r * specularVal;
			color.g += intersectPointColor.g * coEfficients[2] * pointlight->color.g * specularVal;
			color.b += intersectPointColor.b * coEfficients[2] * pointlight->color.b * specularVal;
		}

		// diffuse and specular of SpotLight
		for (SpotLight *spotlight : spotlights)
		{
			Ray lightRay = Ray(spotlight->pos, intersectPoint - spotlight->pos);
			double lightDistance = magnitude(intersectPoint - spotlight->pos);
			if (lightDistance < EPSILON)
				continue;

			double angle = acos(dotProduct(lightRay.dir, spotlight->dir) / (magnitude(lightRay.dir) * magnitude(spotlight->dir))) * (180.0 / pi);
			if (angle > spotlight->cutoffAngle)
				continue;

			bool obscured = false;
			for (Object *object : objects)
			{
				double tempDistance = object->levelZeroIntersect(lightRay, color, 0);
				if (tempDistance > 0 && tempDistance + EPSILON < lightDistance)
				{
					obscured = true;
					break;
				}
			}

			if (obscured)
				continue;

			// diffuse
			Ray normal = getNormalAtIP(intersectPoint, lightRay);
			double lambertVal = max(0.0, dotProduct(normal.dir, lightRay.dir));
			color.r += intersectPointColor.r * coEfficients[1] * spotlight->color.r * lambertVal;
			color.g += intersectPointColor.g * coEfficients[1] * spotlight->color.g * lambertVal;
			color.b += intersectPointColor.b * coEfficients[1] * spotlight->color.b * lambertVal;

			// specular
			Ray reflectedRay = Ray(intersectPoint, lightRay.dir - normal.dir * 2 * dotProduct(lightRay.dir, normal.dir));
			double specularVal = max(0.0, dotProduct(reflectedRay.dir, ray.dir * -1));
			specularVal = pow(specularVal, shine);
			color.r += intersectPointColor.r * coEfficients[2] * spotlight->color.r * specularVal;
			color.g += intersectPointColor.g * coEfficients[2] * spotlight->color.g * specularVal;
			color.b += intersectPointColor.b * coEfficients[2] * spotlight->color.b * specularVal;
		}

		// reflection
		if (recursionLevel > level)
		{
			Ray normal = getNormalAtIP(intersectPoint, ray);
			Ray reflectedRay = Ray(intersectPoint, ray.dir - normal.dir * 2 * dotProduct(ray.dir, normal.dir));
			reflectedRay.start = reflectedRay.start + reflectedRay.dir * EPSILON; // to avoid self intersection

			int nearestObjIndex = -1;
			double t = -1, tMin = 1e9;
			for (int i = 0; i < objects.size(); i++)
			{
				t = objects[i]->intersect(reflectedRay, color, 0);
				if (t > 0 && t < tMin)
				{
					tMin = t;
					nearestObjIndex = i;
				}
			}

			if (nearestObjIndex != -1)
			{
				Color colorTemp(0, 0, 0);
				double t = objects[nearestObjIndex]->intersect(reflectedRay, colorTemp, level + 1);

				color.r += colorTemp.r * coEfficients[3];
				color.g += colorTemp.g * coEfficients[3];
				color.b += colorTemp.b * coEfficients[3];
			}
		}
		return intersectionDistance;
	}
};

class Floor : public Object
{
public:
	int tileCount;
	Floor()
	{
		// floorwidth = 1000, tilewidth = 20
		reference_point = {-500.0, -500.0, 0.0};
		length = 20;
		tileCount = 50;
	}
	Floor(double floorwidth, double tilewidth)
	{
		reference_point = {-floorwidth / 2, -floorwidth / 2, 0.0};
		length = tilewidth;
		tileCount = floorwidth / tilewidth;
	}
	virtual void draw()
	{
		for (int i = 0; i < tileCount; i++)
		{
			for (int j = 0; j < tileCount; j++)
			{
				if (((i + j) % 2) == 0)
					glColor3f(1, 1, 1);
				else
					glColor3f(0, 0, 0);
				glBegin(GL_QUADS);
				{
					glVertex3f(reference_point.x + i * length, reference_point.y + j * length, 0);
					glVertex3f(reference_point.x + (i + 1) * length, reference_point.y + j * length, 0);
					glVertex3f(reference_point.x + (i + 1) * length, reference_point.y + (j + 1) * length, 0);
					glVertex3f(reference_point.x + i * length, reference_point.y + (j + 1) * length, 0);
				}
				glEnd();
			}
		}
	}
	virtual double levelZeroIntersect(Ray ray, Color &color, int level)
	{
		Point normal = {0, 0, 1};
		double denominator = dotProduct(normal, ray.dir);
		if (round(denominator * 100) == 0)
			return -1;
		double t = -1 * dotProduct(normal, ray.start) / denominator;
		Point p = ray.start + ray.dir * t;
		if (p.x <= reference_point.x || p.x >= fabs(reference_point.x) && p.y <= reference_point.y && p.y >= fabs(reference_point.y))
			return -1;
		return t;
	}
	virtual Ray getNormalAtIP(Point point, Ray incidentRay)
	{
		if (incidentRay.dir.z > 0)
			return Ray(point, {0, 0, 1});
		return Ray(point, {0, 0, -1});
	}
	virtual Color getIntersectPointColor(Point point)
	{
		int x = (point.x - reference_point.x) / length;
		int y = (point.y - reference_point.y) / length;

		if (x < 0 || x >= tileCount || y < 0 || y >= tileCount)
		{
			return Color(0, 0, 0);
		}

		if (((x + y) % 2) == 0)
		{
			return Color(1, 1, 1);
		}
		else
		{
			return Color(0, 0, 0);
		}
	}
};

class Sphere : public Object
{
public:
	vector<float> vertices;
	int stackCount = 20;
	int sectorCount = 20;
	Sphere()
	{
	}

	~Sphere()
	{
		vertices.clear();
		vertices.shrink_to_fit();
	}

	Sphere(Point center, double radius)
	{
		reference_point = center;
		length = radius;
	}
	void generateVertices()
	{
		// clear memory of previous arrays
		vertices.clear();

		float x, y, z, xy;
		float sectorStep = 2 * M_PI / sectorCount;
		float stackStep = M_PI / stackCount;
		float sectorAngle, stackAngle;

		for (int i = 0; i <= stackCount; ++i)
		{
			stackAngle = M_PI / 2 - i * stackStep;
			xy = length * cosf(stackAngle);
			z = length * sinf(stackAngle);

			for (int j = 0; j <= sectorCount; ++j)
			{
				sectorAngle = j * sectorStep;

				x = xy * cosf(sectorAngle);
				y = xy * sinf(sectorAngle);

				vertices.push_back(x);
				vertices.push_back(y);
				vertices.push_back(z);
			}
		}
	}
	virtual void draw()
	{
		int k1, k2;

		generateVertices();

		glPushMatrix();
		glTranslatef(reference_point.x, reference_point.y, reference_point.z);
		glColor3f(color.r, color.g, color.b);
		glBegin(GL_QUADS);
		for (int i = 0; i < stackCount; ++i)
		{
			k1 = i * (sectorCount + 1);
			k2 = k1 + sectorCount + 1;

			for (int j = 0; j < sectorCount; ++j, ++k1, ++k2)
			{
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

	virtual double levelZeroIntersect(Ray ray, Color &color, int level)
	{
		ray.start = ray.start - reference_point; // adjust ray origin

		double a = 1;
		double b = 2 * dotProduct(ray.dir, ray.start);
		double c = dotProduct(ray.start, ray.start) - (length * length);

		double d = b * b - 4 * a * c;
		double t = -1;
		if (d < 0)
			return -1;

		if (fabs(a) < 1e-5)
			return -c / b;

		double t1 = (-b - sqrt(d)) / (2 * a);
		double t2 = (-b + sqrt(d)) / (2 * a);

		if (t2 < t1)
			swap(t1, t2);

		if (t1 > 0)
			return t1;
		else if (t2 > 0)
			return t2;
		return -1;
	}

	virtual Ray getNormalAtIP(Point point, Ray incidentRay)
	{
		Point normal = point - reference_point;
		normal = normalize(normal);
		return Ray(point, normal);
	}

	friend istream &operator>>(istream &is, Sphere &sphere)
	{
		is >> sphere.reference_point >> sphere.length;
		is >> sphere.color;
		is >> sphere.coEfficients[0] >> sphere.coEfficients[1] >> sphere.coEfficients[2] >> sphere.coEfficients[3];
		is >> sphere.shine;
		return is;
	}
};

class Triangle : public Object
{
public:
	Point a, b, c;

	Triangle()
	{
	}

	Triangle(Point a, Point b, Point c)
	{
		this->a = a;
		this->b = b;
		this->c = c;
	}

	virtual void draw()
	{
		glColor3f(color.r, color.g, color.b);
		glBegin(GL_TRIANGLES);
		{
			glVertex3f(a.x, a.y, a.z);
			glVertex3f(b.x, b.y, b.z);
			glVertex3f(c.x, c.y, c.z);
		}
		glEnd();
	}

	virtual double levelZeroIntersect(Ray ray, Color &color, int level)
	{
		double A[3][3] = {
			{a.x - b.x, a.x - c.x, ray.dir.x},
			{a.y - b.y, a.y - c.y, ray.dir.y},
			{a.z - b.z, a.z - c.z, ray.dir.z}};
		double detA = determinant(A);
		if (fabs(detA) < EPSILON)
			return -1;

		double beta[3][3] = {
			{a.x - ray.start.x, a.x - c.x, ray.dir.x},
			{a.y - ray.start.y, a.y - c.y, ray.dir.y},
			{a.z - ray.start.z, a.z - c.z, ray.dir.z}};
		double detBeta = determinant(beta);
		double betaVal = detBeta / detA;
		if (betaVal < 0 || betaVal > 1)
			return -1;

		double gamma[3][3] = {
			{a.x - b.x, a.x - ray.start.x, ray.dir.x},
			{a.y - b.y, a.y - ray.start.y, ray.dir.y},
			{a.z - b.z, a.z - ray.start.z, ray.dir.z}};
		double detGamma = determinant(gamma);
		double gammaVal = detGamma / detA;
		if (gammaVal < 0 || gammaVal > 1 - betaVal)
			return -1;

		double t[3][3] = {
			{a.x - b.x, a.x - c.x, a.x - ray.start.x},
			{a.y - b.y, a.y - c.y, a.y - ray.start.y},
			{a.z - b.z, a.z - c.z, a.z - ray.start.z}};
		double detT = determinant(t);
		double tVal = detT / detA;
		if (tVal < 0)
			return -1;
		return tVal;
	}

	virtual Ray getNormalAtIP(Point point, Ray incidentRay)
	{
		Point normal = crossProduct(b - a, c - a);
		normal = normalize(normal);
		if (dotProduct(normal, incidentRay.dir) < 0)
			normal = normal * -1;
		return Ray(point, normal);
	}

	friend istream &operator>>(istream &is, Triangle &triangle)
	{
		is >> triangle.a >> triangle.b >> triangle.c;
		is >> triangle.color;
		is >> triangle.coEfficients[0] >> triangle.coEfficients[1] >> triangle.coEfficients[2] >> triangle.coEfficients[3];
		is >> triangle.shine;
		return is;
	}
};

class General : public Object
{
public:
	double A, B, C, D, E, F, G, H, I, J;

	General()
	{
	}

	virtual void draw()
	{
		return;
	}

	virtual Ray getNormalAtIP(Point point, Ray incidentRay)
	{
		double x, y, z;
		x = 2 * A * point.x + D * point.y + E * point.z + G;
		y = 2 * B * point.y + D * point.x + F * point.z + H;
		z = 2 * C * point.z + E * point.x + F * point.y + I;

		Point dir(x, y, z);

		return Ray(point, dir);
	}

	bool isInside(Point point)
	{
		if (fabs(length) > EPSILON)
		{
			if (point.x < reference_point.x)
				return false;
			if (point.x > reference_point.x + length)
				return false;
		}

		if (fabs(width) > EPSILON)
		{
			if (point.y < reference_point.y)
				return false;
			if (point.y > reference_point.y + width)
				return false;
		}

		if (fabs(height) > EPSILON)
		{
			if (point.z < reference_point.z)
				return false;
			if (point.z > reference_point.z + height)
				return false;
		}

		return true;
	}

	virtual double levelZeroIntersect(Ray ray, Color &color, int level)
	{

		double X0 = ray.start.x;
		double X1 = ray.dir.x;
		double Y0 = ray.start.y;
		double Y1 = ray.dir.y;
		double Z0 = ray.start.z;
		double Z1 = ray.dir.z;

		double a = A * X1 * X1 + B * Y1 * Y1 + C * Z1 * Z1 + D * X1 * Y1 + E * X1 * Z1 + F * Y1 * Z1;
		double b = 2 * A * X0 * X1 + 2 * B * Y0 * Y1 + 2 * C * Z0 * Z1 + D * (X0 * Y1 + X1 * Y0) + E * (X0 * Z1 + X1 * Z0) + F * (Y0 * Z1 + Y1 * Z0) + G * X1 + H * Y1 + I * Z1;
		double c = A * X0 * X0 + B * Y0 * Y0 + C * Z0 * Z0 + D * X0 * Y0 + E * X0 * Z0 + F * Y0 * Z0 + G * X0 + H * Y0 + I * Z0 + J;

		double d = b * b - 4 * a * c;
		if (d < 0)
			return -1;
		if (fabs(a) < 1e-5)
		{
			return -b / c;
		}
		double t1 = (-b - sqrt(d)) / (2 * a);
		double t2 = (-b + sqrt(d)) / (2 * a);

		if (t1 < 0 && t2 < 0)
			return -1;

		if (t2 < t1)
			swap(t1, t2);

		if (t1 > 0)
		{
			Point intersectionPoint = ray.start + ray.dir * t1;
			if (isInside(intersectionPoint))
			{
				return t1;
			}
		}
		if (t2 > 0)
		{
			Point intersectionPoint = ray.start + ray.dir * t2;
			if (isInside(intersectionPoint))
			{
				return t2;
			}
		}

		return -1;
	}

	// input stream
	friend istream &operator>>(istream &in, General &g)
	{
		in >> g.A >> g.B >> g.C >> g.D >> g.E >> g.F >> g.G >> g.H >> g.I >> g.J;
		in >> g.reference_point >> g.length >> g.width >> g.height;

		in >> g.color.r >> g.color.g >> g.color.b; // color
		for (int i = 0; i < 4; i++)
			in >> g.coEfficients[i];
		in >> g.shine;
		return in;
	}
};
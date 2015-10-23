/*!
*
* \author VaHiD AzIzI
*
*/


#include "obstacles/GJK_EPA.h"

SteerLib::GJK_EPA::GJK_EPA()
{
}

//initialize vector d here so the functions intersect and contains origin use same updating d value
Util::Vector d;

//Look at the GJK_EPA.h header file for documentation and instructions
bool SteerLib::GJK_EPA::intersect(float& return_penetration_depth, Util::Vector& return_penetration_vector, const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB)
{
	// run decomposition
	std::vector<std::vector<Util::Vector>> triangleListA = decompose(_shapeA);
	std::vector<std::vector<Util::Vector>> triangleListB = decompose(_shapeB);

	bool collision = false;

	// check GJK for every triangle in B with every triangle in A
	for (int i = 0; i < triangleListA.size(); i++) {
		for (int j = 0; j < triangleListB.size(); j++) {
			if (GJK(triangleListA.at(i), triangleListB.at(j))) {
				collision = true;
			}
		}
	}

	// return result of GJK over decomposed convex sets
	return collision;
}

bool SteerLib::GJK_EPA::GJK(const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB)
{
	// choose any initial d vector
	d = Util::Vector(1, 0, 0);

	// get the first Minkowski Difference point
	Util::Vector minkowskiDiffPoint = support(_shapeA, d) - support(_shapeB, -d);

	std::vector<Util::Vector> simplex;
	// add the point to the simplex
	simplex.push_back(minkowskiDiffPoint);

	// negate d for next point
	d = -d;

	// keep looping, exit conditions are inside loop
	// for convex objects and no rounding errors, this loop should always terminate
	while (true) {

		// add a new point to the simplex because we haven't terminated yet
		minkowskiDiffPoint = support(_shapeA, d) - support(_shapeB, -d);
		simplex.push_back(minkowskiDiffPoint);

		// make sure that the last point added is past the origin
		if (dotProduct3d(minkowskiDiffPoint, d) < 0) {
			// if the point added last is not past the origin in the direction of d, then the Minkowski Sum can't contain the origin 
			// since the last point added is on the edge of the Minkowski Difference
			return false; // There is no collision
		}
		else {
			// since the last point is past the origin, check if the origin is in the current simplex
			if (containsOrigin(simplex)) {
				// if it does then there is a collision
				return true;
			}
			// otherwise containsOrigin updates the simplex and d vector
			// the simplex is updated to be the line segment closest to the origin
			// the d vector is set to be the normal of that line segment in the direction of the origin
		}
	}
}

// support function, returns point on shape with highest projection on vector d
// in other words, the farthest point on the shape in the d direction
Util::Vector SteerLib::GJK_EPA::support(const std::vector<Util::Vector>& _shape, Util::Vector d)
{
	double max = dotProduct3d(_shape[0], d);
	int pos = 0;
	for (int i = 0; i < _shape.size(); i++) {
		if (dotProduct3d(_shape[i], d) > max) {
			max = dotProduct3d(_shape[i], d);
			pos = i;
		}
	}
	return _shape[pos];
}

// returns dot product of two 3d vectors
float SteerLib::GJK_EPA::dotProduct3d(Util::Vector vector1, Util::Vector vector2)
{
	return vector1.x*vector2.x + vector1.y*vector2.y + vector1.z*vector2.z;
}

// returns cross product of two 3d vectors
Util::Vector SteerLib::GJK_EPA::crossProduct3d(Util::Vector vector1, Util::Vector vector2)
{
	float x = vector1.y*vector2.z - vector1.z*vector2.y;
	float y = vector1.z*vector2.x - vector1.x*vector2.z;
	float z = vector1.x*vector2.y - vector1.y*vector2.x;
	return Util::Vector(x,y,z);
}

// check if the simplex contains the origin and updates the simplex and d value if not
bool SteerLib::GJK_EPA::containsOrigin(std::vector<Util::Vector> simplex)
{
	// get the last point added to the simplex
	Util::Vector a = simplex.back();
	
	// compute a0 (same thing as -A)
	Util::Vector a0 = -a;

	if (simplex.size() == 3) {
		// then its the triangle case
		// since in this case the simplex can only have 3 values, b and c are the first 2 values
		Util::Vector b = simplex.at(1);
		Util::Vector c = simplex.at(0);
		
		// compute the edges
		Util::Vector ab = b - a;
		Util::Vector ac = c - a;
		
		// compute the normals to those edges
		Util::Vector abNormal = crossProduct3d(crossProduct3d(ac, ab), ab);
		Util::Vector acNormal = crossProduct3d(crossProduct3d(ab, ac), ac);
		
		// check location of origin using a series of plane tests
		// is the origin on ab side opposite to c?
		if (dotProduct3d(abNormal, a0) > 0) {

			// remove point c
			simplex.erase(simplex.begin());
			
			// set the new direction to abNormal
			d = abNormal;
		}
		else {
			// is the origin on ac side opposite to b?
			if (dotProduct3d(acNormal, a0) > 0) {

				// remove point b
				simplex.erase(simplex.begin()+1);

				// set the new direction to acNormal
				d = acNormal;
			}
			else {
				// since the last direction pointed towards a, and all directions after the initial are made to point towards the origin 
				// the origin cannot be on bc side opposite to a
				// therefore the origin is on ab side towards c, ac side towards b, and bc side towards a
				// which is inside the triangle
				return true;
			}
		}
	}
	else {
		// then its the line segment case, where we simply update the d vector to be normal to the line and pointing towards the origin
		// since in this case the simplex can only have 2 values, b is the first value
		Util::Vector b = simplex.at(0);
		
		// compute AB
		Util::Vector ab = b - a;
		
		// get the Normal to AB in the direction of the origin
		Util::Vector abNormal = crossProduct3d(crossProduct3d(ab, a0), ab);
		
		// set the direction to abNormal
		d = abNormal;
	}
	return false;
}

// uses ear-clipping method, so won't work on polygons with holes, and is relatively slow at O(n^2)
// assumes size of _shape is at least 3
std::vector<std::vector<Util::Vector>> SteerLib::GJK_EPA::decompose(std::vector<Util::Vector> _shape)
{
	// list of triangles that _shape decomposed into
	std::vector<std::vector<Util::Vector>> triangleList;

	// temporary shape that can be modified
	std::vector<Util::Vector> tempShape = _shape;

	// start the search at any point
	int shapePos = 0;

	// loop until there are only 3 points left in tempShape
	// the loop will find ears in tempShape, and remove the ear, so the number of points in tempShape slowly decreases
	while (tempShape.size() > 3) {
		// get the first point, and the two points adjacent to it
		Util::Vector point = tempShape.at(shapePos);

		// get the predecessor, if shapePos is 0, get the back of tempShape
		Util::Vector predecessor;
		if (shapePos == 0) {
			predecessor = tempShape.back();
		}
		else {
			predecessor = tempShape.at(shapePos - 1);
		}

		// get the predecessor, if shapePos is the end, get the front of tempShape
		Util::Vector successor;
		if (shapePos == tempShape.size()) {
			successor = tempShape.front();
		}
		else {
			successor = tempShape.at(shapePos + 1);
		}

		// find the angle of the first point
		Util::Vector line1 = predecessor - point;
		Util::Vector line2 = successor - point;
		float angle = findAngle(line1, line2);
		
		// if the angle is less than 180, it is possible for an ear to exist
		if (angle < M_PI && angle > -M_PI) {
			// make another temporary shape that does not include the point and it's adjacent points
			std::vector<Util::Vector> triangleShape = tempShape;

			int remove = indexOf(triangleShape, predecessor);
			triangleShape.erase(triangleShape.begin() + remove);
			remove = indexOf(triangleShape, point);
			triangleShape.erase(triangleShape.begin() + remove);
			remove = indexOf(triangleShape, successor);
			triangleShape.erase(triangleShape.begin() + remove);

			// check if any points in tempShape are in the triangle formed by the point and it's adjacent points
			// if no points are in the triangle, then the point and it's adjacent points form an ear
			if (!checkTriangle(triangleShape, predecessor, point, successor)) {
				// construct the ear
				std::vector<Util::Vector> ear;
				ear.push_back(predecessor);
				ear.push_back(point);
				ear.push_back(successor);

				// add the ear to the list of decomposed triangles
				triangleList.push_back(ear);

				// remove the point from tempShape since no other triangles in the decomposition can use it
				remove = indexOf(tempShape, point);
				tempShape.erase(tempShape.begin() + remove);
			}
			// there is another point inside the triangle, so it's not an ear
			// move on to the next position
			else {
				shapePos++;
				// if shapePos reaches the end, loop back to the beginning
				if (shapePos == tempShape.size()) {
					shapePos = 0;
				}
			}

		}
		// if the angle is 180, then the point and it's adjacent points forms a line
		// the point can be removed from tempShape
		else if (angle == M_PI || angle == -M_PI) {

			int remove = indexOf(tempShape, point);
			tempShape.erase(tempShape.begin() + remove);
		}
		// if the angle is greater than 180, it can't be an ear
		// move on to the next position
		else {
			shapePos++;
			// if shapePos reaches the end, loop back to the beginning
			if (shapePos == tempShape.size()) {
				shapePos = 0;
			}
		}

	}
	// since there are only 3 points left in tempShape, it must form a triangle in the decomposition
	triangleList.push_back(tempShape);

	return triangleList;
}

// find the angle formed between 2 vectors
float SteerLib::GJK_EPA::findAngle(Util::Vector vector1, Util::Vector vector2)
{
	float dot = dotProduct3d(vector1, vector2);
	float det = vector1.x * vector2.z - vector1.z *vector2.x;
	float theta = std::atan2(det, dot);
	return theta;
}

// returns the index of point in _shape
int SteerLib::GJK_EPA::indexOf(std::vector<Util::Vector> _shape, Util::Vector point)
{
	int pos = -1;
	for (int i = 0; i < _shape.size(); i++) {
		if (_shape.at(i) == point) {
			pos = i;
		}
	}
	return pos;
}

// checks if any point in _shape is inside the triangle formed by point and it's adjacent points
bool SteerLib::GJK_EPA::checkTriangle(std::vector<Util::Vector> _shape, Util::Vector predecessor, Util::Vector point, Util::Vector successor)
{
	for (int i = 0; i < _shape.size(); i++) {
		// checks if shapePoint at i is inside the triangle formed by point and it's adjacent points
		if (pointInTriangle(_shape.at(i), predecessor, point, successor)) {
			return true;
		}
	}
	return false;
}

// checks if shapePoint is inside the triangle formed by point and it's adjacent points
bool SteerLib::GJK_EPA::pointInTriangle(Util::Vector shapePoint, Util::Vector predecessor, Util::Vector point, Util::Vector successor)
{
	bool b1, b2, b3;

	// which side is point on for vector point-predecessor?
	b1 = sign(shapePoint, predecessor, point) < 0.0f;

	// which side is point on for vector successor-point?
	b2 = sign(shapePoint, point, successor) < 0.0f;

	// which side is point on for vector predecessor-successor?
	b3 = sign(shapePoint, successor, predecessor) < 0.0f;

	// if point is to the same side of all 3 vectors, the point is inside the triangle formed by those 3 vectors
	return ((b1 == b2) && (b2 == b3));
}

// use the sign of the determinant of vectors AB and AP to see which side P is on
float SteerLib::GJK_EPA::sign(Util::Vector p, Util::Vector a, Util::Vector b)
{
	return (p.x - a.x) * (b.z - a.z) - (b.x - a.x) * (p.z - a.z);
}

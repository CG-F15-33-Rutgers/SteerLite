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
	std::vector<Util::Vector> masterSimplex;
	// run decomposition
	/*
	std::vector<std::vector<Util::Vector>> triangleListA = decompose(_shapeA);
	std::vector<std::vector<Util::Vector>> triangleListB = decompose(_shapeB);
	*/
	bool collision = false;
	/*
	// check GJK for every triangle in B with every triangle in A
	for (int i = 0; i < triangleListA.size(); i++) {
		for (int j = 0; j < triangleListB.size(); j++) {
			if (GJK(masterSimplex, triangleListA.at(i), triangleListB.at(j))) {
				collision = true;
			}
		}
	}
	*/
	collision = GJK(masterSimplex, _shapeA, _shapeB);
	
	// return result of EPA over total decomposed convex sets
	if (collision) {
		return_penetration_depth = EPA(masterSimplex, return_penetration_vector, _shapeA, _shapeB);
	}
	else {
		return_penetration_depth = 0.0;
	}
	// return result of GJK over decomposed convex sets
	return collision;
}

bool SteerLib::GJK_EPA::GJK(std::vector<Util::Vector>& simplex, const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB)
{
	// choose any initial d vector
	d = Util::Vector(1, 0, 0);

	// get the first Minkowski Difference point
	Util::Vector minkowskiDiffPoint = support(_shapeA, d) - support(_shapeB, -d);

	//std::vector<Util::Vector> simplex;
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

	int count = crossProductCount(_shape);
	
	// if absolute value of count is same as size of shape, then shape is convex
	if (std::abs(count) == _shape.size()) {
		triangleList.push_back(_shape);
		return triangleList;
	}

	bool isCCW;

	// if count > 0, shape is ordered ccw
	if (count > 0) {
		isCCW = true;
	}
	// if count < 0, shape is ordered cw
	else if (count < 0) {
		isCCW = false;
	}

	// start the search at any point
	int shapePos = 0;

	// loop until there are only 3 points left in tempShape
	// the loop will find ears in tempShape, and remove the ear, so the number of points in tempShape slowly decreases
	while (tempShape.size() > 0) {
		
		// if shapePos reaches the end, loop back to the beginning
		if (shapePos == tempShape.size()) {
			shapePos = 0;
		}

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
		if (shapePos == tempShape.size()-1) {
			successor = tempShape.front();
		}
		else {
			successor = tempShape.at(shapePos + 1);
		}

		// find the angle of the first point
		Util::Vector line1 = predecessor - point;
		Util::Vector line2 = successor - point;
		float angle = findAngle(line1, line2);
		
		bool angleBool;
		// if shape is ccw, angles less than pi are positive
		if (isCCW) {
			angleBool = angle < M_PI && angle > 0;
		}
		// if shape is cw, angles less than pi are negative
		else {
			angleBool = angle > -M_PI && angle < 0;
		}

		// if the angle is less than pi, it is possible for an ear to exist
		if (angleBool) {
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
			}

		}
		// if the angle is pi, then the point and it's adjacent points forms a line
		// the point can be removed from tempShape
		else if (angle == M_PI || angle == -M_PI) {
			int remove = indexOf(tempShape, point);
			tempShape.erase(tempShape.begin() + remove);
		}
		// if the angle is greater than pi, it can't be an ear
		// move on to the next position
		else {
			shapePos++;
		}

		if (tempShape.size() == 3) {
			// since there are only 3 points left in tempShape, it must form a triangle in the decomposition
			triangleList.push_back(tempShape);
			tempShape.clear();
		}

	}
	

	
	for (int i = 0; i < triangleList.size(); i++) {
		std::cout << "triangle " << i << "\n";
		for (int j = 0; j < triangleList.at(i).size(); j++) {
			std::cout << "point " << triangleList.at(i).at(j) << "\n";
		}
	}
	

	return triangleList;
}

// check how shape is ordered, either counterclockwise or clockwise, by using cross product on adjacent edges
// check whether shape is convex by using cross product on adjacent edges
int SteerLib::GJK_EPA::crossProductCount(std::vector<Util::Vector> _shape)
{
	int counter = 0;

	Util::Vector point;
	Util::Vector predecessor;
	Util::Vector successor;

	Util::Vector cross;

	// check all pairs of adjacent edges
	for (int i = 0; i < _shape.size(); i++) {
		// get the point
		point = _shape.at(i);

		// get the predecessor
		// if i is 0, then the predecessor is the last element
		if (i == 0) {
			predecessor = _shape.back();
		}
		else {
			predecessor = _shape.at(i - 1);
		}

		// get the successor
		// if i is the last element, then the successor is the first element
		if (i == _shape.size() - 1) {
			successor = _shape.front();
		}
		else {
			successor = _shape.at(i + 1);
		}

		// get the cross product of the two adjacent edges
		Util::Vector edge1 = point - predecessor;
		Util::Vector edge2 = successor - point;
		cross = crossProduct3d(edge1, edge2);

		// if the cross product is positive, increase counter
		if (cross.y > 0) {
			counter++;
		}
		// if the cross product is negative, decrease counter
		else if(cross.y < 0){
			counter--;
		}
		// if cross product is 0, then vectors are parallel or at least one vector is 0
		// counter does not move in this case
	}
	// counter is number of positive cross products - negative cross products
	// if more positive cross products, shape is counterclockwise
	// if more negative cross products, shape is clockwise
	// if all cross products have same sign, shape is convex
	return counter;
}

// find the angle formed between 2 vectors
// range is from -pi to pi
// angle is positive for ccw angles, and negative for cw angles
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

Util::Vector SteerLib::GJK_EPA::getNormal(Util::Vector vec)
{
	Util::Vector normal;
	normal.x = vec.z;
	normal.z = -vec.x;
	normal.y = vec.y;

	
	return normal;
}

Util::Vector SteerLib::GJK_EPA::normalize(Util::Vector vec)
{
	Util::Vector normalized;
	float length = vec.length();


	normalized = *(new Util::Vector((vec.x / length), (vec.y / length), (vec.z / length)));

	return normalized;
}

float SteerLib::GJK_EPA::EPA(std::vector<Util::Vector>& simplex, Util::Vector& return_penetration_vector, const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB)
{
	
	float prevDist = FLT_MAX;

	
	const float THRESHOLD = 0.01;
	while (true)
	{
		return_penetration_vector = closestMinkEdge(simplex, prevDist);
		

		Util::Vector supportPoint = support(_shapeA, return_penetration_vector) - support(_shapeB, -return_penetration_vector);

		// find distance from origin to support point
		float dist = dotProduct3d(supportPoint, return_penetration_vector);
		
		// now, we compare the difference between distances by the THRESHOLD
		if (std::abs(dist - prevDist) < THRESHOLD) {
			return dist;
		}
		else {
			prevDist = dist;
			simplex.push_back(supportPoint);
		}
	}
}
Util::Vector SteerLib::GJK_EPA::closestMinkEdge(std::vector<Util::Vector>& simplex, float& closestDist)
{ 
	closestDist = FLT_MAX;
	Util::Vector closestEdgeNormal;

	for (int i = 0; i < simplex.size(); i++) {

		int j;
		float dist;
		if (i + 1 == simplex.size()) {
			j = 0;
		}
		else {
			j = i + 1;
		}

		

		// get an edge from mink diff by getting next two points
		Util::Vector veci = simplex.at(i);
		Util::Vector vecj = simplex.at(j);
	
		Util::Vector edge = vecj - veci;
		
		Util::Vector edgeVec = crossProduct3d(crossProduct3d(edge, veci), edge);

		// try and avoid problems with origin being too close to edge
		if (edgeVec.length() <= 0.1)
		{
			edgeVec = getNormal(edge);
			
		}

		edgeVec = normalize(edgeVec);
		
		
		
		
		
		
		// get distance of the edgeNormal
		dist = dotProduct3d(veci, edgeVec);
		
	

		if (std::abs(dist) < closestDist) {
			closestDist = dist;
			closestEdgeNormal = edgeVec;
		}


	}
	return closestEdgeNormal;
} 
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
	return GJK(_shapeA, _shapeB);
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

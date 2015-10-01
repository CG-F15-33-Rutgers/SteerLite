//
// Copyright (c) 2015 Mahyar Khayatkhoei
// Copyright (c) 2009-2014 Shawn Singh, Glen Berseth, Mubbasir Kapadia, Petros Faloutsos, Glenn Reinman
// See license.txt for complete license.
//

#include <algorithm>
#include <vector>
#include <util/Geometry.h>
#include <util/Curve.h>
#include <util/Color.h>
#include <util/DrawLib.h>
#include "Globals.h"

using namespace Util;

Curve::Curve(const CurvePoint& startPoint, int curveType) : type(curveType)
{
	controlPoints.push_back(startPoint);
}

Curve::Curve(const std::vector<CurvePoint>& inputPoints, int curveType) : type(curveType)
{
	controlPoints = inputPoints;
	sortControlPoints();
}

// Add one control point to the vector controlPoints
void Curve::addControlPoint(const CurvePoint& inputPoint)
{
	controlPoints.push_back(inputPoint);
	sortControlPoints();
}

// Add a vector of control points to the vector controlPoints
void Curve::addControlPoints(const std::vector<CurvePoint>& inputPoints)
{
	for (int i = 0; i < inputPoints.size(); i++)
		controlPoints.push_back(inputPoints[i]);
	sortControlPoints();
}

void Curve::drawCurve(Color curveColor, float curveThickness, int window)
{
#ifdef ENABLE_GUI

	//================DELETE THIS PART AND THEN START CODING===================
	/*
	static bool flag = false;
	if (!flag)
	{
	std::cerr << "ERROR>>>>Member function drawCurve is not implemented!" << std::endl;
	flag = true;
	}
	*/
	//=========================================================================

	// Robustness: make sure there is at least two control point: start and end points
	if (!checkRobust())
		return;

	float endTime = Curve::controlPoints.back().time;
	Util::Point newPosition;

	// Move on the curve from t=0 to t=finalPoint, using window as step size, and linearly interpolate the curve points
	for (int i = 0; i < endTime + window; i = i + window) {
		Util::Point firstPoint = newPosition;
		Curve::calculatePoint(newPosition, i);
		DrawLib::drawLine(firstPoint, newPosition, curveColor, curveThickness);
	}


	return;
#endif
}

// helper function for sortControlPoints()
// Author: Michael Jiao (mj498)
bool compareControlPoints(CurvePoint i, CurvePoint j)
{
	return (i.time < j.time);
}

// Sort controlPoints vector in ascending order: min-first
void Curve::sortControlPoints()
{
	//================DELETE THIS PART AND THEN START CODING===================
	/*
	static bool flag = false;
	if (!flag)
	{
		std::cerr << "ERROR>>>>Member function sortControlPoints is not implemented!" << std::endl;
		flag = true;
	}
	*/
	//=========================================================================

	// sort using a lambda which defines how to order the points with respect to time.
	std::sort(Curve::controlPoints.begin(), Curve::controlPoints.end(), compareControlPoints);

	return;
}

// Calculate the position on curve corresponding to the given time, outputPoint is the resulting position
bool Curve::calculatePoint(Point& outputPoint, float time)
{
	// Robustness: make sure there is at least two control point: start and end points
	if (!checkRobust())
		return false;

	// Define temporary parameters for calculation
	unsigned int nextPoint;
	float normalTime, intervalTime;

	// Find the current interval in time, supposing that controlPoints is sorted (sorting is done whenever control points are added)
	if (!findTimeInterval(nextPoint, time))
		return false;

	// Calculate position at t = time on curve
	if (type == hermiteCurve)
	{
		outputPoint = useHermiteCurve(nextPoint, time);
	}
	else if (type == catmullCurve)
	{
		outputPoint = useCatmullCurve(nextPoint, time);
	}

	// Return
	return true;
}

// Check Roboustness
bool Curve::checkRobust()
{
	//================DELETE THIS PART AND THEN START CODING===================
	/*
	static bool flag = false;
	if (!flag)
	{
		std::cerr << "ERROR>>>>Member function checkRobust is not implemented!" << std::endl;
		flag = true;
	}
	*/
	//=========================================================================
	if (Curve::controlPoints.size() < 2) {
		return false;
	}

	return true;
}

// Find the current time interval (i.e. index of the next control point to follow according to current time)
// unsigned int& nextPoint is the index of next point in vector 
bool Curve::findTimeInterval(unsigned int& nextPoint, float time)
{
	//================DELETE THIS PART AND THEN START CODING===================
	/*
	static bool flag = false;
	if (!flag)
	{
		std::cerr << "ERROR>>>>Member function findTimeInterval is not implemented!" << std::endl;
		flag = true;
	}*/
	//=========================================================================
	for (int i = 0; i < Curve::controlPoints.size() - 1; i++) {
		if (time >= Curve::controlPoints.at(i).time && time <= Curve::controlPoints.at(i + 1).time) {
			nextPoint = (i + 1);
			return true;
		}
	}

	return false;
}

// Implement Hermite curve
// const unsigned int nextPoint is index in vector that contains nextPoint
Point Curve::useHermiteCurve(const unsigned int nextPoint, const float time)
{
	Point newPosition;
	float normalTime, intervalTime;

	//================DELETE THIS PART AND THEN START CODING===================
	/*
	static bool flag = false;
	if (!flag)
	{
	std::cerr << "ERROR>>>>Member function useHermiteCurve is not implemented!" << std::endl;
	flag = true;
	}
	*/
	//=========================================================================

	Point p1 = Curve::controlPoints.at(nextPoint - 1).position;
	Vector t1 = Curve::controlPoints.at(nextPoint - 1).tangent;
	Point p2 = Curve::controlPoints.at(nextPoint).position;
	Vector t2 = Curve::controlPoints.at(nextPoint).tangent;

	// Calculate time interval, and normal time required for later curve calculations
	float time1 = Curve::controlPoints.at(nextPoint - 1).time;
	float time2 = Curve::controlPoints.at(nextPoint).time;

	normalTime = time - time1;
	intervalTime = time2 - time1;

	// scale s to go from 0 to 1
	float s = normalTime / intervalTime;

	// calculate basis functions
	float h1 = 2 * s*s*s - 3 * s*s + 1;
	float h2 = -2 * s*s*s + 3 * s*s;
	float h3 = s*s*s - 2 * s*s + s;
	float h4 = s*s*s - s*s;

	// Calculate position at t = time on Hermite curve
	newPosition = h1*p1 + h2*p2 + h3*t1 + h4*t2;

	// Return result
	return newPosition;
}

// Implement Catmull-Rom curve
Point Curve::useCatmullCurve(const unsigned int nextPoint, const float time)
{
	Point newPosition;
	float normalTime, intervalTime;
	//================DELETE THIS PART AND THEN START CODING===================
	/*
	static bool flag = false;
	if (!flag)
	{
		std::cerr << "ERROR>>>>Member function useCatmullCurve is not implemented!" << std::endl;
		flag = true;
	}
	*/
	//=========================================================================

	Vector t1;
	Vector t2;
	Point p1 = Curve::controlPoints.at(nextPoint - 1).position;
	Point p2 = Curve::controlPoints.at(nextPoint).position;
	Point p3;
	// Calculate time interval, and normal time required for later curve calculations
	float time1 = Curve::controlPoints.at(nextPoint - 1).time;
	float time2 = Curve::controlPoints.at(nextPoint).time;

	normalTime = time - time1;
	intervalTime = time2 - time1;

	// scale s to go from 0 to 1
	float s = normalTime / intervalTime;
	// Calculate position at t = time on Catmull-Rom curve

	// need to check if we are at boundary conditions
	if (nextPoint == 1) {
		p3 = Curve::controlPoints.at(nextPoint + 1).position;
		t1 = 2 * ((p2 - p1) / s) - ((p3 - p2) / (2 * s));
		t2 = (p3 - p1) / (2 * s);
	}
	else if (nextPoint == (Curve::controlPoints.size() - 1)) {
		p3 = Curve::controlPoints.at(nextPoint - 2).position;
		t1 = (p2 - p3) / (2 * s);
		t2 = t1 - 2 * ((p1 - p3) / s);
	}
	else {
		p3 = Curve::controlPoints.at(nextPoint - 2).position;
		t1 = (p2 - p3) / (2 * s);
		p3 = Curve::controlPoints.at(nextPoint + 1).position;
		t2 = (p3 - p1) / (2 * s);
	}

	// Calculate hermite curve again using computed tangents
	// calculate basis functions
	float h1 = 2 * s*s*s - 3 * s*s + 1;
	float h2 = -2 * s*s*s + 3 * s*s;
	float h3 = s*s*s - 2 * s*s + s;
	float h4 = s*s*s - s*s;

	// Calculate position at t = time on Hermite curve
	newPosition = h1*p1 + h2*p2 + h3*t1 + h4*t2;

	// Return result
	return newPosition;
}
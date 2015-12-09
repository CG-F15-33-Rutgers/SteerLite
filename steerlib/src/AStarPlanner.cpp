//
// Copyright (c) 2009-2015 Glen Berseth, Mubbasir Kapadia, Shawn Singh, Petros Faloutsos, Glenn Reinman, Rahul Shome
// See license.txt for complete license.
//


#include <vector>
#include <stack>
#include <set>
#include <map>
#include <iostream>
#include <algorithm> 
#include <functional>
#include <queue>
#include <math.h>
#include "planning/AStarPlanner.h"


#define COLLISION_COST  1000
#define GRID_STEP  1
#define OBSTACLE_CLEARANCE 0
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

namespace SteerLib
{
	AStarPlanner::AStarPlanner() {}

	AStarPlanner::~AStarPlanner() {}

	bool AStarPlanner::canBeTraversed(int id)
	{
		double traversal_cost = 0;
		int current_id = id;
		unsigned int x, z;
		gSpatialDatabase->getGridCoordinatesFromIndex(current_id, x, z);
		int x_range_min, x_range_max, z_range_min, z_range_max;

		//std::cout << "\n X is " << x;
		//std::cout << "\n Z is " << z;
		

		x_range_min = MAX(x - OBSTACLE_CLEARANCE, 0);
		x_range_max = MIN(x + OBSTACLE_CLEARANCE, gSpatialDatabase->getNumCellsX());

		z_range_min = MAX(z - OBSTACLE_CLEARANCE, 0);
		z_range_max = MIN(z + OBSTACLE_CLEARANCE, gSpatialDatabase->getNumCellsZ());

		//std::cout << "\n numcells X is " << gSpatialDatabase->getNumCellsX();
		//std::cout << "\n numcells Z is " << gSpatialDatabase->getNumCellsZ();

		//std::cout << "\n X min is " << x_range_min;
		//std::cout << "\n Z min is " << z_range_min;

		//std::cout << "\n X max is " << x_range_max;
		//std::cout << "\n Z max is " << z_range_max;


		for (float i = x_range_min; i <= x_range_max; i += GRID_STEP)
		{
			for (float j = z_range_min; j <= z_range_max; j += GRID_STEP)
			{
				int index = gSpatialDatabase->getCellIndexFromGridCoords(i, j);

			//	std::cout << "\n index is " << index;

				traversal_cost += gSpatialDatabase->getTraversalCost(index);

				//std::cout << "\n cost at index is " << gSpatialDatabase->getTraversalCost(index);
				//std::cout << "\n total traversal cost is " << traversal_cost;


			}
		}

		if (traversal_cost > COLLISION_COST)
			return false;
		return true;
	}



	Util::Point AStarPlanner::getPointFromGridIndex(int id)
	{
		Util::Point p;
		gSpatialDatabase->getLocationFromIndex(id, p);
		return p;
	}

	bool AStarPlanner::computePath(std::vector<Util::Point>& agent_path, Util::Point start, Util::Point goal, SteerLib::GridDatabase2D * _gSpatialDatabase, bool append_to_path)
	{
		// Initialize the data structures before search, starting with start node in open set
		std::vector<SteerLib::AStarPlannerNode*> openset;
		std::vector<SteerLib::AStarPlannerNode*> closedset;
		

		std::vector<SteerLib::AStarPlannerNode*> nodeList;

		int explored = 0; // number of nodes expanded
		const int weight = 8;

		gSpatialDatabase = _gSpatialDatabase;

		double f_start = weight * heuristic(start, goal);
		SteerLib::AStarPlannerNode* startNode = new SteerLib::AStarPlannerNode(start, 0, f_start, NULL);
		openset.push_back(startNode);
		++explored;
		//TODO
		//std::cout << "\nIn A*";

		int count = 0;


		while (!openset.empty()) {

			/*
			count++;
			if (count > 3) {
				break;
			}
			*/


			int currentIndex = getLowestIndex(openset);
			SteerLib::AStarPlannerNode* current = openset.at(currentIndex);

			//std::cout << "\nCurrent is " << current->point << " with f value " << current->f;
			
			//std::cout << "\nGoal is " << goal;

			float current_diff_x = std::abs(current->point.x - goal.x);
			float current_diff_z = std::abs(current->point.z - goal.z);

			if (current_diff_x < GRID_STEP && current_diff_z < GRID_STEP) {
				// reconstruct path. 
				// IMPLEMENT PATH METHOD
				
				std::vector<Util::Point> new_path = reconstruct_path(start, current);
				new_path.push_back(goal);

				//std::cout << "\nConstructing path";

				if (append_to_path) {
					//std::cout << "\nAppending to path";
					agent_path.insert(agent_path.end(), new_path.begin(), new_path.end());
				}
				else {
					//std::cout << "\nNot appending to path";
					agent_path = new_path;
				}

				for (int i = 0; i < agent_path.size(); i++) {
					//std::cout << "\nMoving along path " << agent_path.at(i);
				}
				printf("\nTotal length of solution path(as number of nodes): %d", new_path.size());
				printf("\nTotal amount of expanded nodes: %d", explored);
				

				return true;
			}

			// move current to closedset
			openset.erase(openset.begin()+currentIndex);
			closedset.push_back(current);

			// next, we visit all neighbors of current
			
			std::vector<Util::Point> neighborList;
			bool isNeighborList = AStarPlanner::getNeighborNodes(current, neighborList, nodeList);

			//std::cout << "\nVisiting Neighbors";
			
			for (int i = 0; i < neighborList.size(); i++) {
				Util::Point neighbor = neighborList.at(i);
				double tentative_g_score = 0;
				
				//std::cout << "\nExploring neighbor " << neighbor.x << "," << neighbor.z;
				
				if (findNode(neighbor, closedset) != -1) {
					continue;
				}
				
				if (AStarPlanner::areDiagonal(current->point, neighbor)) {
					// change this to increase diagonal cost
					tentative_g_score = current->g + GRID_STEP;
				}
				else {
					tentative_g_score = current->g + GRID_STEP;
				}
				
				int nodeIndex = findNode(neighbor, nodeList);
				SteerLib::AStarPlannerNode* neighborNode = nodeList.at(nodeIndex);
				
				if (tentative_g_score < neighborNode->g) {

					neighborNode->parent = current;
					neighborNode->g = tentative_g_score;
					neighborNode->f = tentative_g_score + (weight * heuristic(neighbor, goal));

					openset.push_back(neighborNode);
					++explored;
					//std::cout << "\nUpdated neighbor is " << neighbor.x << "," << neighbor.z << "with f value " << neighborNode->f;

				}

				

			}
			

		}

		return false;
	}

	int AStarPlanner::getLowestIndex(std::vector<SteerLib::AStarPlannerNode*> openset) {
		double min = openset.front()->f;
		double gmin = openset.front()->g;
		int gIndex = 0;
		int index = 0;
		for (int i = 0; i < openset.size(); i++) {
			if (openset.at(i)->f < min) {
				min = openset.at(i)->f;
				gmin = openset.at(i)->g;
				gIndex = i;
				index = i;
			}
			else if (openset.at(i)->f == min) {
				// SWITCH comparison here
				if (gmin < openset.at(i)->g) {
					index = gIndex;
				}
				else {
					index = i;
				}
			}
		}
		return index;
	}

	std::vector<Util::Point> AStarPlanner::reconstruct_path(Util::Point start, SteerLib::AStarPlannerNode* current)
	{
		std::vector<Util::Point> totalPath;
		current->point = getPointFromGridIndex(gSpatialDatabase->getCellIndexFromLocation(current->point));

		totalPath.push_back(current->point);
		while (!(current->parent->point.x == start.x && current->parent->point.z == start.z)) {
			current = current->parent;
			current->point = getPointFromGridIndex(gSpatialDatabase->getCellIndexFromLocation(current->point));
			totalPath.push_back(current->point);
		}
		current->point = getPointFromGridIndex(gSpatialDatabase->getCellIndexFromLocation(current->point));
		totalPath.push_back(current->point);

		std::vector<Util::Point> reversePath;
		int size = totalPath.size();
		for (int i = 0; i < size; i++) {
			reversePath.push_back(totalPath.back());
			totalPath.pop_back();
		}

		return reversePath;
	}

	bool AStarPlanner::getNeighborNodes(SteerLib::AStarPlannerNode* current, std::vector<Util::Point>& neighborList, std::vector<SteerLib::AStarPlannerNode*>& nodeList)
	{
		// First, get coordinates of current point
		float x_range_min, x_range_max, z_range_min, z_range_max;
		int currIndex;
		currIndex = gSpatialDatabase->getCellIndexFromLocation(current->point);
		
		


		x_range_min = current->point.x - GRID_STEP;
		x_range_max = current->point.x + GRID_STEP;

		z_range_min = current->point.z - GRID_STEP;
		z_range_max = current->point.z + GRID_STEP;

		for (float i = x_range_min; i <= x_range_max; i += GRID_STEP) {
			for (float j = z_range_min; j <= z_range_max; j += GRID_STEP) {
				if (!(i == current->point.x && j == current->point.z)) {
					if (AStarPlanner::canBeTraversed(currIndex)) {
						Util::Point currentPoint = Util::Point(i, 0, j);
						
						//std::cout << "\nAdding neighbor point " << currentPoint;

						if (!containsNode(currentPoint, nodeList)) {
							SteerLib::AStarPlannerNode* neighbor = new AStarPlannerNode(currentPoint, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), current);
							nodeList.push_back(neighbor);
						}

						neighborList.push_back(currentPoint);

					}
				}
			}
		}
		return false;
	}

	bool AStarPlanner::containsNode(Util::Point point, std::vector<SteerLib::AStarPlannerNode*> nodeList) {
		for (int i = 0; i < nodeList.size(); i++) {
			if (nodeList.at(i)->point.x == point.x && nodeList.at(i)->point.z == point.z) {
				return true;
			}
		}
		return false;
	}

	int AStarPlanner::findNode(Util::Point point, std::vector<SteerLib::AStarPlannerNode*> nodeList) {
		for (int i = 0; i < nodeList.size(); i++) {
			if (nodeList.at(i)->point.x == point.x && nodeList.at(i)->point.z == point.z) {
				return i;
			}
		}
		return -1;
	}

	// Wrapper function that allows us to constantly set which heuristic function to employ
	double AStarPlanner::heuristic(Util::Point node, Util::Point goal)
	{
		// change which function this returns to retune A*
		// return AStarPlanner::manhattanHeuristic(node, goal);
		return AStarPlanner::euclideanHeuristic(node, goal);
	}

	double AStarPlanner::manhattanHeuristic(Util::Point node, Util::Point goal)
	{
		return (std::abs(node.x - goal.x) + std::abs(node.z - goal.z));
	}

	double AStarPlanner::euclideanHeuristic(Util::Point node, Util::Point goal)
	{
		return std::sqrt(std::pow((node.x - goal.x), 2.0) + std::pow((node.z - goal.z), 2.0));
	}

	bool AStarPlanner::areDiagonal(Util::Point curr, Util::Point neighbor)
	{
		int x = curr.x - neighbor.x;
		int z = curr.z - neighbor.z;

		if (std::abs(x) == std::abs(z)) {
			return true;
		}
		return false;
	}
}
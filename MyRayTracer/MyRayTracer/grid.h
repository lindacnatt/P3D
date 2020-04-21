//below two statements are for the h. file and implemented class
#ifndef GRID_H
#define GRID_H

#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>
#include "boundingBox.h"

class Object;

#define EPSILON			0.0001f


using namespace std;

double _bin_size = 1;

// defining classes to avoid C2061 by https://www.dreamincode.net/forums/topic/293640-problem-with-error-c2061-syntax-error-identifier-xxx/


class Grid
{
public:
	Grid(void);
	//~Grid(void);

	int getNumObjects();
	void addObject(Object* o);
	Object* getObject(unsigned int index);

	void Build(void);   // set up grid cells --> expanded in the grid.cpp file

	ShadeRec Traverse(const Ray& ray, double& tmin, ShadeRec& sr);
	
		// Don't get it why we should use a Bool for this? In the Original Algo they also return a list;
	std::vector<Vector> Traverse(Ray& ray);
	
	

private:
	vector<Object*> objects;
	vector<vector<Object*> > cells;

	int nx, ny, nz; // number of cells in the x, y, and z directions
	float m = 2.0f; // factor that allows to vary the number of cells

	Vector find_min_bounds(void);

	Vector find_max_bounds(void);
	/*//Setup function for Grid traversal ; Initialisation phase -> identifying Voxel which includes Ray Origin
	bool Init_Traverse(Ray& ray, int& ix, int& iy, int& iz, double& dtx, double& dty, double& dtz, double& tx_next, double& ty_next, double& tz_next,
		int& ix_step, int& iy_step, int& iz_step, int& ix_stop, int& iy_stop, int& iz_stop) {
	
		Vector ray_origin = ray.origin;
		Vector ray_direction = ray.direction;

		double dtx = dtx;
		double dty = dty;
		double dtz = dtz;

		double tx_next = tx_next;
		double ty_next = ty_next;
		double tz_next = tz_next;

		int ix_step = ix_step;
		int iy_step = iz_step;
		int iz_step = iz_step;

		int ix_stop = ix_stop;
		int iy_stop = iy_stop;
		int iz_stop = iz_stop;


		};*/

	AABB& bbox; // have to use '&' to avoid error c2079: uses undefined class
};
#endif
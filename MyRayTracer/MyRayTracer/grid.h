//below two statements are for the h. file and implemented class
#ifndef GRID_H
#define GRID_H

#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>
#include "scene.h"
#include "boundingBox.h";

#define EPSILON			0.0001f


using namespace std;

double _bin_size = 1;

class Grid
{
public:
	Grid(void);
	//~Grid(void);

	int getNumObjects();
	void addObject(Object* o);
	Object* getObject(unsigned int index);

	void Build(void);   // set up grid cells --> expanded in the grid.cpp file

	// bool Traverse(Ray& ray, Object** hitobject, Vector& hitpoint);  //(const Ray& ray, double& tmin, ShadeRec& sr)
	ShadeRec Traverse(const Ray& ray, double& tmin, ShadeRec& sr);
	
		// Don't get it why we should use a Bool for this? In the Original Algo they also return a list;
	std::vector<Vector> Traverse(Ray& ray);
	
	

private:
	vector<Object*> objects;
	vector<vector<Object*> > cells;

	int nx, ny, nz; // number of cells in the x, y, and z directions
	float m = 2.0f; // factor that allows to vary the number of cells

	Vector find_min_bounds(void) {
		AABB bbox;
		Vector pkHV;
		int num_objects = objects.size();
		for (int j = 0; j < num_objects; j++) {
			bbox = objects[j]->GetBoundingBox();

			if (bbox.min.x < pkHV.x); {
				pkHV.x = bbox.min.x;
			}
			if (bbox.min.y < pkHV.y); {
				pkHV.y = bbox.min.y;
			}
			if (bbox.min.y < pkHV.z); {
				pkHV.z = bbox.min.z;
			};
		};

		pkHV.x -= EPSILON; pkHV.y -= EPSILON; pkHV.z = EPSILON;

		return (pkHV);
	};

	Vector find_max_bounds(void) {
		AABB bbox;
		Vector pkHV;
		int num_objects = objects.size();
		for (int j = 0; j < num_objects; j++) {
			bbox = objects[j]->GetBoundingBox();

			if (bbox.max.x > pkHV.x); {
				pkHV.x = bbox.max.x;
			}
			if (bbox.max.y > pkHV.y); {
				pkHV.y = bbox.max.y;
			}
			if (bbox.min.y > pkHV.z); {
				pkHV.z = bbox.max.z;
			};
		};

		pkHV.x -= EPSILON; pkHV.y -= EPSILON; pkHV.z = EPSILON;

		return (pkHV);
	};

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

	AABB bbox;
};
#endif
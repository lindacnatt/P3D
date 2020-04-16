#ifndef GRID_H
#define GRID_H

#include <vector>
#include <cmath>
#include "scene.h"

using namespace std;

class Grid
{
public:
	Grid(void);
	//~Grid(void);

	int getNumObjects();
	void addObject(Object* o);
	Object* getObject(unsigned int index);

	void Build(void);   // set up grid cells

	bool Traverse(Ray& ray, Object** hitobject, Vector& hitpoint);  //(const Ray& ray, double& tmin, ShadeRec& sr)
	bool Traverse(Ray& ray);  //Traverse for shadow ray

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

	//Setup function for Grid traversal
	bool Init_Traverse(Ray& ray, int& ix, int& iy, int& iz, double& dtx, double& dty, double& dtz, double& tx_next, double& ty_next, double& tz_next,
		int& ix_step, int& iy_step, int& iz_step, int& ix_stop, int& iy_stop, int& iz_stop);

	AABB bbox;
};
#endif
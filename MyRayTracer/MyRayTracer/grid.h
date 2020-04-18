//below two statements are for the h. file and implemented class
#ifndef GRID_H
#define GRID_H

#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>
#include "scene.h"

using namespace std;

double _bin_size = 1;



class Grid
{
public:
	Grid(void);
	//~Grid(void);

	int getNumObjects();
	void addObject(Object* o);
	Object* getObject(unsigned int index) {};

	void Build(void);   // set up grid cells --> expanded in the grid.cpp file

	// bool Traverse(Ray& ray, Object** hitobject, Vector& hitpoint);  //(const Ray& ray, double& tmin, ShadeRec& sr)

	bool Traverse(Ray& ray, Object** hitobject, Vector& hitpoint) {//(const Ray& ray, double& tmin, ShadeRec& sr) 
		//Does not have to be changed I guess, just the way how we should use it in the RayTracer Code?
		// 
	};
	//shall we implement it in our RayTracer with bool Traverse(Ray& ray, Object** hitobject, Vector& hitindex) or as bool Traverse(const Ray& ray, double& tmin) without ShadeRec& sr since only in book and not in our RT 


	// Don't get it why we should use a Bool for this? In the Original Algo they also return a list;
	std::vector<Vector> Traverse(Ray& ray) {  //Traverse for shadow ray and return found intersected Object (?); Code based on https://github.com/francisengelmann/fast_voxel_traversal
		std::vector<Vector> visited_voxels;
		Vector ray_start = ray.origin;
		Vector ray_end = ray.direction;
		Vector ray_dir = (ray.direction - ray.origin).normalize();

		// This id of the first/current voxel which is hitted by the ray.
		// Using floor to round down
		// the implicit int-casting will round up for negative numbers.
		Vector current_voxel(std::floor(ray_start.x / _bin_size), std::floor(ray_start.y / _bin_size), std::floor(ray_start.z / _bin_size));
		Vector last_voxel(std::floor(ray_end.x / _bin_size), std::floor(ray_end.y / _bin_size), std::floor(ray_end.z / _bin_size));

		// In which direction are the voxel ids incremented.
		double stepX = (ray_dir.x >= 0) ? 1 : -1; // correct
		double stepY = (ray_dir.y >= 0) ? 1 : -1; // correct
		double stepZ = (ray_dir.z >= 0) ? 1 : -1; // correct

		  // Distance along the ray to the next voxel border from the current position (tMaxX, tMaxY, tMaxZ).
		double next_voxel_boundary_x = (current_voxel.x + stepX) * _bin_size; // correct
		double next_voxel_boundary_y = (current_voxel.y + stepY) * _bin_size; // correct
		double next_voxel_boundary_z = (current_voxel.z + stepZ) * _bin_size; // correct

		// tMaxX, tMaxY, tMaxZ -- distance until next intersection with voxel-border
		// the value of t at which the ray crosses the first vertical voxel boundary
		double tMaxX = (ray_dir.x != 0) ? (next_voxel_boundary_x - ray_start.x) / ray_dir.x : DBL_MAX; //
		double tMaxY = (ray_dir.y != 0) ? (next_voxel_boundary_y - ray_start.y) / ray_dir.y : DBL_MAX; //
		double tMaxZ = (ray_dir.z != 0) ? (next_voxel_boundary_z - ray_start.z) / ray_dir.z : DBL_MAX; //

		// tDeltaX, tDeltaY, tDeltaZ --
		// how far along the ray we must move for the horizontal component to equal the width of a voxel
		// the direction in which we traverse the grid
		// can only be FLT_MAX if we never go in that direction
		double tDeltaX = (ray_dir.x != 0) ? _bin_size / ray_dir.x * stepX : DBL_MAX;
		double tDeltaY = (ray_dir.y != 0) ? _bin_size / ray_dir.y * stepY : DBL_MAX;
		double tDeltaZ = (ray_dir.z != 0) ? _bin_size / ray_dir.z * stepZ : DBL_MAX;

		Vector diff(0, 0, 0);
		bool neg_ray = false;
		if (current_voxel.x != last_voxel.x && ray_dir.x < 0) { diff.x--; neg_ray = true; }
		if (current_voxel.y != last_voxel.y && ray_dir.y < 0) { diff.y--; neg_ray = true; }
		if (current_voxel.z != last_voxel.z && ray_dir.z < 0) { diff.z--; neg_ray = true; }
		visited_voxels.push_back(current_voxel);
		if (neg_ray) {
			current_voxel + diff;
			visited_voxels.push_back(current_voxel);
		};



		while !(std::equal(last_voxel, current_voxel)) {
			if (tMaxX < tMaxY) {
				if (tMaxX < tMaxZ) {
					current_voxel.x += stepX;
					tMaxX += tDeltaX;
				}
				else {
					current_voxel.z += stepZ;
					tMaxZ += tDeltaZ;
				}
			}
			else {
				if (tMaxY < tMaxZ) {
					current_voxel.y += stepY;
					tMaxY += tDeltaY;
				}
				else {
					current_voxel.z += stepZ;
					tMaxZ += tDeltaZ;
				}
			}
			visited_voxels.push_back(current_voxel);
		}
		return visited_voxels;
	};


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

	//Setup function for Grid traversal ; Initialisation phase -> identifying Voxel which includes Ray Origin
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


		};

	AABB bbox;
};
#endif
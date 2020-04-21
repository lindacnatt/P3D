#include <iostream>
#include <cfloat>
#include <vector>
#include <cmath>
#include "scene.h"
#include "maths.h"
#include "grid.h"


using namespace std;



int Grid::getNumObjects() {
	objects.size();
};
void Grid::addObject(Object* o) {
	objects.push_back(o);
};

Object* Grid::getObject(unsigned int index) {
	if (index >= 0 && index < objects.size())
		return objects[index];
	return NULL;
};

Vector Grid::find_min_bounds(void) {
	//AABB bbox; -> not needed to be initialised again after been declared on the ground
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

Vector Grid::find_max_bounds(void) {
	//AABB bbox; -> not needed to be initialised again after been declared on the ground
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


void Grid::Build(void) {
		// set up grid cells

		//find the minimum and maximum coordinates of the grid
		Vector p0 = find_min_bounds();
		Vector p1 = find_max_bounds();

		//store the grid coordinates in the bounding box
		AABB bbox;
		bbox.min.x = p0.x; bbox.min.y = p0.y; bbox.min.z = p0.z;
		bbox.max.x = p1.x; bbox.max.y = p1.y; bbox.max.z = p1.z;

		//compute the number of cells in the x-, y-, and z-directions
		float wx = p1.x - p0.x; //grid-extent in x-direction
		float wy = p1.y - p0.y; //grid-extent in y-direction
		float wz = p1.z - p0.z; //grid-extent in z-direction
		float multiplier = 2.0; // 8 times more cells than objects (m=2)
		
		float s = pow(wx * wy * wz / getNumObjects(), 0.3333333); //Pixel size
		int nx = multiplier * wx / s + 1;
		int ny = multiplier * wy / s + 1;
		int nz = multiplier * wz / s + 1;

		//set up the array of cells with null pointers

		int num_cells = nx * ny * nz;
		cells.reserve(getNumObjects());

		for (int j = 0; j < num_cells; j++) {
			cells.push_back({NULL}); // are {} needed?
		};

		//set up a temporary array to hold the number of objects stored in each cell

		vector<int> counts;
		counts.reserve(num_cells);

		for (int j = 0; j < num_cells; j++) {
			cells.push_back({0}); //are {} needed ?
		};

		// put objects into the cells 

		AABB obj_box; // object's bounding box
		int index; //cells array index

		for (int j = 0; j < getNumObjects(); j++) {
			obj_box = objects[j]->GetBoundingBox();

			//compute cell indices for the corners of the bounding box of the object
			int ixmin = clamp((obj_box.min.x - p0.x) * nx / (p1.x - p0.x), 0, nx - 1);
			int iymin = clamp((obj_box.min.y - p0.y) * ny / (p1.y - p0.y), 0, ny - 1);
			int izmin = clamp((obj_box.min.z - p0.z) * nz / (p1.z - p0.z), 0, nz - 1);

			int ixmax = clamp((obj_box.max.x - p0.x) * nx / (p1.x - p0.x), 0, nx - 1);
			int iymax = clamp((obj_box.max.y - p0.y) * ny / (p1.y - p0.y), 0, ny - 1);
			int izmax = clamp((obj_box.max.z - p0.z) * nz / (p1.z - p0.z), 0, nz - 1);

			//add the object to the cells

			for (int iz = izmin; iz <= izmax; iz++) { //cells in z direction
				for (int iy = iymin; iy <= iymax; iy++) { // cells in y direction
					for (int ix = ixmin; ix <= ixmax; ix++) { //cells in x direction
						index = ix + nx * iy + nx * ny * iz;

						/*if (counts[index] == 0){
							cells[index] = getObject(j);
							counts[index] += 1; //index = 1
						}
						else {
							if (counts[index] == 1) {
								//construct a compound object, object that stores arrays of any other geometric objects.
								Compound* compound_ptr = new Compound; //need alternative object or implement that one
								// add the object already in cell
								compound_ptr->addObject(cells[index]);
								//add the new project
								compound_ptr->addObject(getObject(j));

								// store compound in current cell
								cells[index] = compound_ptr;
								// index 2
								counts[index] += 1;
							}*/
							 //counts[index] > 1			// maybe only use this one
								// just add current object
							cells[index].push_back(getObject(j)); //changed add_object method from above to the push_back method
								

								//for statistics only
								counts[index] += 1;

							};

						};
					};
				};
				//erase Compound::Objects but don't delete them
				objects.erase(objects.begin(), objects.end());
				//erase temporary counts vector
				counts.erase(counts.begin(), counts.end());
				
			};
			
			ShadeRec Grid::Traverse(const Ray& ray, double& tmin, ShadeRec& sr) { 	//or shall we implement it in our RayTracer with bool Traverse(Ray& ray, Object** hitobject, Vector& hitindex)? alias hit_bare_bones_objects
					//ShadeRec sr(*this); -> disappears because is a Parameter
				float t;
				//double tmin = 100000; // = kHugeValue -> disappears because is a Parameter
				int num_objects = getNumObjects();

				for (int j = 0; j < num_objects; j++) {
					if (getObject(j)->intercepts(ray, t, sr) && (t < tmin)) {
						sr.hit_an_object = true;
						tmin = t;
						sr.color = objects[j]->GetMaterial()->GetDiffColor(); // eigentlich objects[j]->get_color() ; doch wie bekomme ich den Color von den Objekten
							//->get_color(); TBD -> can we take Set_Material()
					};
					return sr; // SR-Object saves Informations on how to shade ray-object at the hitting point
				}
			};

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


				while ((last_voxel != current_voxel)) {
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



	




			//Setup function for Grid traversal --> privtae in the grid.h file, doesnt have to be implemented
bool Init_Traverse(Ray & ray, int& ix, int& iy, int& iz, double& dtx, double& dty, double& dtz, double& tx_next, double& ty_next, double& tz_next,
int& ix_step, int& iy_step, int& iz_step, int& ix_stop, int& iy_stop, int& iz_stop);

		// use Amantides and Woo Algo

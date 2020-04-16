
#include <vector>
#include <cmath>
#include "scene.h"
#include "grid.h"
#include "maths.h"


using namespace std;

int Grid::getNumObjects() {
	objects.size();
};
void Grid::addObject(Object* o) {
	objects.push_back(o);
};

Object* getObject(unsigned int index) {
	//to be done
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
		
		float s = pow(wx * wy * wz / getNumObjects(), 0.3333333);
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

						if (counts[index] == 0){
							cells[index] = objects[j];
							counts[index] += 1; //index = 1
						}
						else {
							if (counts[index] == 1) {
								//construct a compound object, object that stores arrays of any other geometric objects.
								Compound* compound_ptr = new Compound; //need alternative object or implement that one
								// add the object already in cell
								compound_ptr->addObject(cells[index]);
								//add the new project
								compound_ptr->addObject(objects[j]);

								// store compound in current cell
								cells[index] = compound_ptr;
								// index 2
								counts[index] += 1;
							}
							else { //counts[index] > 1
								// just add current object
								cells[index]->addObject(objects[j]);
								

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

};


	
bool Traverse(Ray & ray, Object * *hitobject, Vector & hitpoint);  //(const Ray& ray, double& tmin, ShadeRec& sr)
bool Traverse(Ray & ray);  //Traverse for shadow ray




	//Setup function for Grid traversal
	bool Init_Traverse(Ray & ray, int& ix, int& iy, int& iz, double& dtx, double& dty, double& dtz, double& tx_next, double& ty_next, double& tz_next,
		int& ix_step, int& iy_step, int& iz_step, int& ix_stop, int& iy_stop, int& iz_stop);

	// use Amantides and Woo Algo


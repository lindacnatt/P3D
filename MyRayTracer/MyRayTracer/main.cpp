///////////////////////////////////////////////////////////////////////
//
// P3D Course
// (c) 2019 by Jo�o Madeiras Pereira
//Ray Tracing P3F scenes and drawing points with Modern OpenGL
//
///////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <chrono>
#include <conio.h>

#include <GL/glew.h>
#include <GL/freeglut.h>
#include <IL/il.h>

#include "scene.h"
#include "grid.h"
#include "maths.h"
#include "sampler.h"

#include <array>

#define CAPTION "Whitted Ray-Tracer"

#define VERTEX_COORD_ATTRIB 0
#define COLOR_ATTRIB 1

#define MAX_DEPTH 4

//Enable OpenGL drawing.  
bool drawModeEnabled = true;

//Draw Mode: 0 - point by point; 1 - line by line; 2 - full frame at once
int draw_mode = 1;

// Points defined by 2 attributes: positions which are stored in vertices array and colors which are stored in colors array
float* colors;
float* vertices;
int size_vertices;
int size_colors;

//Array of Pixels to be stored in a file by using DevIL library
uint8_t* img_Data;

GLfloat m[16];  //projection matrix initialized by ortho function

GLuint VaoId;
GLuint VboId[2];

GLuint VertexShaderId, FragmentShaderId, ProgramId;
GLint UniformId;

Scene* scene = NULL;
int RES_X, RES_Y;

int WindowHandle = 0;


// bools and global values
bool antialiasing = true;
bool softshadows = false;
float softshadowsample_1;
float softshadowsample_2;
bool depthoffield = true;
bool gridon = false;
Grid* grid = NULL;

Color rayTracing(Ray ray, int depth, float ior_1)  //index of refraction of medium 1 where the ray is travelling
{
	// Variables: Ray (includes origin and direction, depth and index of refraction
	//Calculate intersection  intercepts() functions returns true or false
	float dist;
	float tNear = INFINITY;
	int hitIndex;
	bool hit = false;
	Vector phit;
	Vector normal;
	if (!gridon) {
		for (int obj_i = 0; obj_i < scene->getNumObjects(); obj_i += 1)  //Looping through all objects to check if there is an intersection
		{
			if (scene->getObject(obj_i)->intercepts(ray, dist) == true && dist < tNear)					//check if ray is intercepting object, put obj_i in list of hitObjects; t in intercepts function checks that it is the closest t
			{
				tNear = dist;	// continues to make tNear smaller
				hitIndex = obj_i; // stores index for future uses when the object is needed
				hit = true;

			};
		};
	}
	if (gridon) {
		// TRAVERSE GRID
		// GET HITINDEX or HIToBJEKT and T
	}

	// Calculate color
	if (!hit)
	{
		return scene->GetBackgroundColor();
	}
	else
	{
		Vector phit = ray.direction * tNear + ray.origin;
		Vector normal = scene->getObject(hitIndex)->getNormal(phit);
		Color finalColor;
		Vector bias = normal * EPSILON; //To escape from self-intersections
		Vector I = (ray.direction * (-1));  // symmetric of the incident ray direction
		int shadownum = 0;
		// LOOP THROUGH LIGHTS
		for (int i = 0; i < scene->getNumLights(); i += 1)  // for every i, starting from 0, to the amount of lights (-1 for correct indexing), stepping 1 index per loop
		{
			// SHADOWING, including antialiasing and not
			Vector lightsource = scene->getLight(i)->position;
			if (softshadows && antialiasing)		// If antialiased and softshadows is activated, use the randomized sample to find a different point on the area light 
			{
				Vector lb = Vector(1, 0, 0);
				Vector la = Vector(0, 1, 0);
				lightsource = lightsource + lb * softshadowsample_1 + la * softshadowsample_2;  // r = c + eps1a + eps2b
			}
			Vector L = (lightsource - phit);		// unit light vector from hit point to light source
			float shadowDist = L.length();
			L.normalize();
			bool shadow = false;
			float ts;
			float intensity = L.normalize() * normal;  // angle between normal & light direction
			Vector H = (L + I).normalize();  //Half-Vector from Blinn model
			float NH = normal * H;
			if (NH < 0) NH = 0;
			Ray shadowRay = Ray((phit + bias), L);

			if (intensity > 0)  // intensity of light
			{
				if (gridon) {
					// TRAVERSE GRID
					//if hit, set shadow = true;
				}
				else {
					for (int sobj_i = 0; sobj_i < scene->getNumObjects(); sobj_i += 1)  //Looping through all objects to check if there is an intersection
					{
						if (scene->getObject(sobj_i)->intercepts(shadowRay, ts) == true && ts < shadowDist)  //check if ray towards source light is intercepting object
						{
							shadow = true;
							shadownum += 1;
							break;
						};
					}
				}

			}
			if (!shadow) //if no object is blocking source light
			{
				float Kd = scene->getObject(hitIndex)->GetMaterial()->GetDiffuse();
				Color diffuse_color = scene->getObject(hitIndex)->GetMaterial()->GetDiffColor() * Kd * intensity * scene->getLight(i)->color; // Calculate diffuse
				float Ks = scene->getObject(hitIndex)->GetMaterial()->GetSpecular();
				float spec = pow(NH, scene->getObject(hitIndex)->GetMaterial()->GetShine());
				Color specular_color = scene->getObject(hitIndex)->GetMaterial()->GetSpecColor() * Ks * spec * scene->getLight(i)->color;  // Calculate specular
				Color specDiffColor = diffuse_color + specular_color;  //combine colors
				finalColor += specDiffColor;
				
			}
		}

		if (depth < MAX_DEPTH) 
		{	
			// reflection
			if (scene->getObject(hitIndex)->GetMaterial()->GetReflection() > 0)   // !! If reflective component is bigger than 0, means it is reflective?
			{
				Vector V = (ray.direction) * (-1);		// math from slides
				Vector rRefl = normal * 2 * (V * normal) - V;
				Vector reflectionOrigin = phit + bias;
				Ray reflRay = Ray(reflectionOrigin, rRefl.normalize());
				Color reflColor = rayTracing(reflRay, depth + 1, ior_1); //iteration
				reflColor = reflColor * scene->getObject(hitIndex)->GetMaterial()->GetReflection() * scene->getObject(hitIndex)->GetMaterial()->GetSpecColor(); //  reduce rColor by the specular reflection coefficient	
				finalColor += reflColor;
			};

			// refraction (som p� scratch)
			if (scene->getObject(hitIndex)->GetMaterial()->GetTransmittance() > 0.0f)
			{
				float cosi = std::fmax(-1.f, std::fmin(1.f, ray.direction * normal));  // minus from other example, what does it mean?
				float etai = ior_1, etat = scene->getObject(hitIndex)->GetMaterial()->GetRefrIndex();  // comparing ior of hit object with ior of "previous" object (or air depending on what happened earlier
				
				Vector n = normal;
				bool outside;

				if (cosi < 0) {
					cosi = (-1) * cosi;
					outside = true;
				}
				else{
					outside = false;
					etat = 1.0;
					n = normal * (-1);  // if the ray is inside the object, swap the indices and invert the normal to get the correct result
				}

				//  Fresnels
				Vector v = (ray.direction) * (-1);
				Vector vt = n * (v * n) - v;
				float sini = (vt).length();
				float sint = (etai / etat) * sini;

				if (sint*sint <= 1) { //If not total reflection
				
					float cost = sqrt(1 - (sint) * (sint));
					vt.normalize();
					Vector refrDir = vt * sint - n * cost;
					Vector refrOrig = outside ? phit - bias : phit + bias;  // direction of usage of bias depends on whether or not point is outside or inside object

					Ray refractiveRay = Ray(refrOrig, refrDir);
					Color refrColor = rayTracing(refractiveRay, depth + 1, etat);
				
					cosi = fabsf(cosi);
					float Rp = pow((etai * cosi - etat * cost) / (etai * cosi + etat * cost), 2);
					float Rs = pow((etai * cost - etat * cosi) / (etai * cost + etat * cosi), 2);
					float KR = (Rs + Rp) * 0.5;

					finalColor += refrColor * (1 - KR);
				};
			}				
		};
		return (finalColor.clamp());
	}
}

/////////////////////////////////////////////////////////////////////// ERRORS

bool isOpenGLError() {
	bool isError = false;
	GLenum errCode;
	const GLubyte* errString;
	while ((errCode = glGetError()) != GL_NO_ERROR) {
		isError = true;
		errString = gluErrorString(errCode);
		std::cerr << "OpenGL ERROR [" << errString << "]." << std::endl;
	}
	return isError;
}

void checkOpenGLError(std::string error)
{
	if (isOpenGLError()) {
		std::cerr << error << std::endl;
		exit(EXIT_FAILURE);
	}
}

/////////////////////////////////////////////////////////////////////// SHADERs

const GLchar* VertexShader =
{
	"#version 430 core\n"

	"in vec2 in_Position;\n"
	"in vec3 in_Color;\n"
	"uniform mat4 Matrix;\n"
	"out vec4 color;\n"

	"void main(void)\n"
	"{\n"
	"	vec4 position = vec4(in_Position, 0.0, 1.0);\n"
	"	color = vec4(in_Color, 1.0);\n"
	"	gl_Position = Matrix * position;\n"

	"}\n"
};

const GLchar* FragmentShader =
{
	"#version 430 core\n"

	"in vec4 color;\n"
	"out vec4 out_Color;\n"

	"void main(void)\n"
	"{\n"
	"	out_Color = color;\n"
	"}\n"
};

void createShaderProgram()
{
	VertexShaderId = glCreateShader(GL_VERTEX_SHADER);
	glShaderSource(VertexShaderId, 1, &VertexShader, 0);
	glCompileShader(VertexShaderId);

	FragmentShaderId = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(FragmentShaderId, 1, &FragmentShader, 0);
	glCompileShader(FragmentShaderId);

	ProgramId = glCreateProgram();
	glAttachShader(ProgramId, VertexShaderId);
	glAttachShader(ProgramId, FragmentShaderId);

	glBindAttribLocation(ProgramId, VERTEX_COORD_ATTRIB, "in_Position");
	glBindAttribLocation(ProgramId, COLOR_ATTRIB, "in_Color");

	glLinkProgram(ProgramId);
	UniformId = glGetUniformLocation(ProgramId, "Matrix");

	checkOpenGLError("ERROR: Could not create shaders.");
}

void destroyShaderProgram()
{
	glUseProgram(0);
	glDetachShader(ProgramId, VertexShaderId);
	glDetachShader(ProgramId, FragmentShaderId);

	glDeleteShader(FragmentShaderId);
	glDeleteShader(VertexShaderId);
	glDeleteProgram(ProgramId);

	checkOpenGLError("ERROR: Could not destroy shaders.");
}

/////////////////////////////////////////////////////////////////////// VAOs & VBOs


void createBufferObjects()
{
	glGenVertexArrays(1, &VaoId);
	glBindVertexArray(VaoId);
	glGenBuffers(2, VboId);
	glBindBuffer(GL_ARRAY_BUFFER, VboId[0]);

	/* S� se faz a aloca��o dos arrays glBufferData (NULL), e o envio dos pontos para a placa gr�fica
	� feito na drawPoints com GlBufferSubData em tempo de execu��o pois os arrays s�o GL_DYNAMIC_DRAW */
	glBufferData(GL_ARRAY_BUFFER, size_vertices, NULL, GL_DYNAMIC_DRAW);
	glEnableVertexAttribArray(VERTEX_COORD_ATTRIB);
	glVertexAttribPointer(VERTEX_COORD_ATTRIB, 2, GL_FLOAT, 0, 0, 0);

	glBindBuffer(GL_ARRAY_BUFFER, VboId[1]);
	glBufferData(GL_ARRAY_BUFFER, size_colors, NULL, GL_DYNAMIC_DRAW);
	glEnableVertexAttribArray(COLOR_ATTRIB);
	glVertexAttribPointer(COLOR_ATTRIB, 3, GL_FLOAT, 0, 0, 0);

	// unbind the VAO
	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	//	glDisableVertexAttribArray(VERTEX_COORD_ATTRIB); 
	//	glDisableVertexAttribArray(COLOR_ATTRIB);
	checkOpenGLError("ERROR: Could not create VAOs and VBOs.");
}

void destroyBufferObjects()
{
	glDisableVertexAttribArray(VERTEX_COORD_ATTRIB);
	glDisableVertexAttribArray(COLOR_ATTRIB);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);

	glDeleteBuffers(1, VboId);
	glDeleteVertexArrays(1, &VaoId);
	checkOpenGLError("ERROR: Could not destroy VAOs and VBOs.");
}

void drawPoints()
{
	glBindVertexArray(VaoId);
	glUseProgram(ProgramId);

	glBindBuffer(GL_ARRAY_BUFFER, VboId[0]);
	glBufferSubData(GL_ARRAY_BUFFER, 0, size_vertices, vertices);
	glBindBuffer(GL_ARRAY_BUFFER, VboId[1]);
	glBufferSubData(GL_ARRAY_BUFFER, 0, size_colors, colors);

	glUniformMatrix4fv(UniformId, 1, GL_FALSE, m);

	if (draw_mode == 0) glDrawArrays(GL_POINTS, 0, 1);
	else if (draw_mode == 1) glDrawArrays(GL_POINTS, 0, RES_X);
	else glDrawArrays(GL_POINTS, 0, RES_X * RES_Y);
	glFinish();

	glUseProgram(0);
	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	checkOpenGLError("ERROR: Could not draw scene.");
}

ILuint saveImgFile(const char* filename) {
	ILuint ImageId;

	ilEnable(IL_FILE_OVERWRITE);
	ilGenImages(1, &ImageId);
	ilBindImage(ImageId);

	ilTexImage(RES_X, RES_Y, 1, 3, IL_RGB, IL_UNSIGNED_BYTE, img_Data /*Texture*/);
	ilSaveImage(filename);

	ilDisable(IL_FILE_OVERWRITE);
	ilDeleteImages(1, &ImageId);
	if (ilGetError() != IL_NO_ERROR)return ilGetError();

	return IL_NO_ERROR;
}

/////////////////////////////////////////////////////////////////////// CALLBACKS

// Render function by primary ray casting from the eye towards the scene's objects

void renderScene()
{
	int index_pos = 0;
	int index_col = 0;
	unsigned int counter = 0;
	srand(time(0)); //using current time as seed
	if (gridon == true) {
		//Grid grid = Grid().Build;  // initialize/build grid
	};
	if (softshadows && !antialiasing)  
	{
		int num_lights = scene->getNumLights();
		for (int li = 0; li < num_lights; li++)		// for every light in the scene...
		{
			Light* light = scene->getLight(li);
			double arealight_samples = 9;			// ...9 "new lights" are created...
			for (float alx = 0; alx < sqrt(arealight_samples); alx++) {
				for (float aly = 0; aly < sqrt(arealight_samples); aly++) {
					Vector newpos = light->position + Vector(0.1f, 0, 0) * alx + Vector(0, 0.1f, 0) * aly;  // ... with a new position...
					Color newcolor = light->color * (1 / (arealight_samples));   // and color intensity 1/9th of the original color
					if ((light->position.x != newpos.x) && (light->position.y != newpos.y)) 
					{ 
						scene->addLight(new Light(newpos, newcolor));
					}
				}
			}
			light->color = light->color * (1 / arealight_samples);  //  change color of original light
		}
	}

	for (int y = 0; y < RES_Y; y++)
	{
		for (int x = 0; x < RES_X; x++)
		{
			Ray ray = Ray(Vector(), Vector());  // defined as blank ray, gets changed later
			Color color;
			const int n = 2;	// defined by us to decide how much to split the pixel in (jittering samples)
			Vector pixel;		//viewport coordinates
			float r[n*n];		// define array s and r for soft shadows
			float s[n*n];
			//depthoffield = false;
			//antialiasing = false;
			
			// Antialiasing with jittering ... combined solution from book and https://www.scratchapixel.com/code.php?id=13&origin=/lessons/3d-basic-rendering/introduction-to-shading
			if (antialiasing) {
				Color c = Color(0.0, 0.0, 0.0);
				for (int p = 0; p <= n - 1; p++) {
					for (int q = 0; q <= n - 1; q++) {
						// Vector pixel;  //viewport coordinates
						pixel.x = x + (p + rand_float()) / float(n);
						pixel.y = y + (q + rand_float()) / float(n);

						if (softshadows)
						{
							for (int i = 0; i < n; i++) {		// create n number of randomized numbers and store in array r and s
								float rand_rf = rand_float();
								float rand_sf = rand_float();
								r[i] = rand_rf;
								s[i] = rand_sf;
							}
							//shuffling array s
							for (int i = pow(n, 2); i > 1; i--) {
								int j = rand() % i;
								std::swap(s[j], s[i]);
							}
							// set softshadowsamples (global variables) to values from r and s, to be used in raytracing for randomizing a point on the area light
							// this could also be done by including r and s as values that could be optional parameters in raytracer, but only necessary to avoid using global variables
							softshadowsample_1 = s[p];
							softshadowsample_2 = r[q];
							
						}
						if (depthoffield) {
							Vector lens_sample = sample_unit_disk() * scene->GetCamera()->GetAperture() / 2.0f;
							ray = scene->GetCamera()->PrimaryRay(lens_sample, pixel);
						}
						if (!depthoffield) { ray = scene->GetCamera()->PrimaryRay(pixel); }
						c = c + rayTracing(ray, 1, 1.0);
						
					};
				}
				float nsqr = 1 / pow(n, 2);
				color = c * nsqr;	// divide by n^2 to not overexpose each pixel due to jittering
			}
	
			else {
				//Not using antialiasing
				pixel.x = x + 0.5f;
				pixel.y = y + 0.5f;
				ray = scene->GetCamera()->PrimaryRay(pixel);
				color = rayTracing(ray, 1, 1.0);
			}
			//color = scene->GetBackgroundColor(); //just for the template

			img_Data[counter++] = u8fromfloat((float)color.r());
			img_Data[counter++] = u8fromfloat((float)color.g());
			img_Data[counter++] = u8fromfloat((float)color.b());

			if (drawModeEnabled) {
				vertices[index_pos++] = (float)x;
				vertices[index_pos++] = (float)y;
				colors[index_col++] = (float)color.r();

				colors[index_col++] = (float)color.g();

				colors[index_col++] = (float)color.b();


				if (draw_mode == 0) {  // drawing point by point
					drawPoints();
					index_pos = 0;
					index_col = 0;
				}
			}
		}
		if (draw_mode == 1 && drawModeEnabled) {  // drawing line by line
			drawPoints();
			index_pos = 0;
			index_col = 0;
		}
	}
	if (draw_mode == 2 && drawModeEnabled)        //full frame at once
		drawPoints();

	printf("Drawing finished!\n");

	if (saveImgFile("RT_Output.png") != IL_NO_ERROR) {
		printf("Error saving Image file\n");
		exit(0);
	}
	printf("Image file created\n");
	glFlush();
}


// Callback function for glutCloseFunc
void cleanup()
{
	destroyShaderProgram();
	destroyBufferObjects();
}

void ortho(float left, float right, float bottom, float top,
	float nearp, float farp)
{
	m[0 * 4 + 0] = 2 / (right - left);
	m[0 * 4 + 1] = 0.0;
	m[0 * 4 + 2] = 0.0;
	m[0 * 4 + 3] = 0.0;
	m[1 * 4 + 0] = 0.0;
	m[1 * 4 + 1] = 2 / (top - bottom);
	m[1 * 4 + 2] = 0.0;
	m[1 * 4 + 3] = 0.0;
	m[2 * 4 + 0] = 0.0;
	m[2 * 4 + 1] = 0.0;
	m[2 * 4 + 2] = -2 / (farp - nearp);
	m[2 * 4 + 3] = 0.0;
	m[3 * 4 + 0] = -(right + left) / (right - left);
	m[3 * 4 + 1] = -(top + bottom) / (top - bottom);
	m[3 * 4 + 2] = -(farp + nearp) / (farp - nearp);
	m[3 * 4 + 3] = 1.0;
}

void reshape(int w, int h)
{
	glClear(GL_COLOR_BUFFER_BIT);
	glViewport(0, 0, w, h);
	ortho(0, (float)RES_X, 0, (float)RES_Y, -1.0, 1.0);
}

void processKeys(unsigned char key, int xx, int yy)
{
	switch (key) {

	case 27:
		glutLeaveMainLoop();
		break;

	}
}

/////////////////////////////////////////////////////////////////////// SETUP

void setupCallbacks()
{
	glutKeyboardFunc(processKeys);
	glutCloseFunc(cleanup);
	glutDisplayFunc(renderScene);
	glutReshapeFunc(reshape);
}

void setupGLEW() {
	glewExperimental = GL_TRUE;
	GLenum result = glewInit();
	if (result != GLEW_OK) {
		std::cerr << "ERROR glewInit: " << glewGetString(result) << std::endl;
		exit(EXIT_FAILURE);
	}
	GLenum err_code = glGetError();
	printf("Vendor: %s\n", glGetString(GL_VENDOR));
	printf("Renderer: %s\n", glGetString(GL_RENDERER));
	printf("Version: %s\n", glGetString(GL_VERSION));
	printf("GLSL: %s\n", glGetString(GL_SHADING_LANGUAGE_VERSION));
}

void setupGLUT(int argc, char* argv[])
{
	glutInit(&argc, argv);

	glutInitContextVersion(4, 3);
	glutInitContextFlags(GLUT_FORWARD_COMPATIBLE);
	glutInitContextProfile(GLUT_CORE_PROFILE);

	glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_GLUTMAINLOOP_RETURNS);

	glutInitWindowPosition(640, 100);
	glutInitWindowSize(RES_X, RES_Y);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA);
	glDisable(GL_DEPTH_TEST);
	WindowHandle = glutCreateWindow(CAPTION);
	if (WindowHandle < 1) {
		std::cerr << "ERROR: Could not create a new rendering window." << std::endl;
		exit(EXIT_FAILURE);
	}
}


void init(int argc, char* argv[])
{
	setupGLUT(argc, argv);
	setupGLEW();
	std::cerr << "CONTEXT: OpenGL v" << glGetString(GL_VERSION) << std::endl;
	glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
	createShaderProgram();
	createBufferObjects();
	setupCallbacks();

}


void init_scene(void)
{
	char scenes_dir[70] = "P3D_Scenes/";
	char input_user[50];
	char scene_name[70];

	while (true) {
		cout << "Input the Scene Name: ";
		cin >> input_user;
		strcpy_s(scene_name, sizeof(scene_name), scenes_dir);
		strcat_s(scene_name, sizeof(scene_name), input_user);

		ifstream file(scene_name, ios::in);
		if (file.fail()) {
			printf("\nError opening P3F file.\n");
		}
		else
			break;
	}

	scene = new Scene();
	scene->load_p3f(scene_name);
	RES_X = scene->GetCamera()->GetResX();
	RES_Y = scene->GetCamera()->GetResY();
	printf("\nResolutionX = %d  ResolutionY= %d.\n", RES_X, RES_Y);

	// Pixel buffer to be used in the Save Image function
	img_Data = (uint8_t*)malloc(3 * RES_X * RES_Y * sizeof(uint8_t));
	if (img_Data == NULL) exit(1);
}

int main(int argc, char* argv[])
{
	//Initialization of DevIL 
	if (ilGetInteger(IL_VERSION_NUM) < IL_VERSION)
	{
		printf("wrong DevIL version \n");
		exit(0);
	}
	ilInit();

	int ch;
	if (!drawModeEnabled) {

		do {
			init_scene();
			auto timeStart = std::chrono::high_resolution_clock::now();
			renderScene();  //Just creating an image file
			auto timeEnd = std::chrono::high_resolution_clock::now();
			auto passedTime = std::chrono::duration<double, std::milli>(timeEnd - timeStart).count();
			printf("\nDone: %.2f (sec)\n", passedTime / 1000);

			cout << "\nPress 'y' to render another image or another key to terminate!\n";
			delete(scene);
			free(img_Data);
			ch = _getch();
		} while ((toupper(ch) == 'Y'));
	}

	else {   //Use OpenGL to draw image in the screen
		init_scene();
		if (draw_mode == 0) { // draw image point by point
			size_vertices = 2 * sizeof(float);
			size_colors = 3 * sizeof(float);
			printf("DRAWING MODE: POINT BY POINT\n\n");
		}
		else if (draw_mode == 1) { // draw image line by line
			size_vertices = 2 * RES_X * sizeof(float);
			size_colors = 3 * RES_X * sizeof(float);
			printf("DRAWING MODE: LINE BY LINE\n\n");
		}
		else if (draw_mode == 2) { // draw full frame at once
			size_vertices = 2 * RES_X * RES_Y * sizeof(float);
			size_colors = 3 * RES_X * RES_Y * sizeof(float);
			printf("DRAWING MODE: FULL IMAGE\n\n");
		}
		else {
			printf("Draw mode not valid \n");
			exit(0);
		}
		vertices = (float*)malloc(size_vertices);
		if (vertices == NULL) exit(1);

		colors = (float*)malloc(size_colors);
		if (colors == NULL) exit(1);

		/* Setup GLUT and GLEW */
		init(argc, argv);
		glutMainLoop();
	}

	free(colors);
	free(vertices);
	printf("Program ended normally\n");
	exit(EXIT_SUCCESS);
}
///////////////////////////////////////////////////////////////////////
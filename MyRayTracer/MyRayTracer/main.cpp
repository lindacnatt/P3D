///////////////////////////////////////////////////////////////////////
//
// P3D Course
// (c) 2019 by João Madeiras Pereira
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


// antialiasing bool
bool antialiasing;

Color rayTracing(Ray ray, int depth, float ior_1)  //index of refraction of medium 1 where the ray is travelling
{
	// Variables: Ray (includes origin and direction, depth and index of refraction
	//INSERT HERE YOUR CODE
	//Calculate intersection  intercepts() functions returns true or false
	float dist;
	float tNear = INFINITY;
	int hitIndex = -1;
	Vector phit;
	Vector normal;
	for (int obj_i = 0; obj_i < scene->getNumObjects() ; obj_i += 1)  //Looping through all objects to check if there is an intersection
	{
		if (scene->getObject(obj_i)->intercepts(ray, dist) == true && dist < tNear)					//check if ray is intercepting object, put obj_i in list of hitObjects; t in intercepts function checks that it is the closest t
		{
			tNear = dist;	// continues to make tNear smaller
			hitIndex = obj_i;	// stores index for future uses when the object is needed
			
		};
	};
	if (hitIndex < 0)
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

		// LOOP THROUGH LIGHTS
		for (int i = 0; i < scene->getNumLights(); i += 1)  // for every i, starting from 0, to the amount of lights (-1 for correct indexing), stepping 1 index per loop
		{
			/*if (antialiasing)
			{
				int n = 5;
				for (int p = 0; p <= n - 1; p++) {
					for (int q = 0; q <= n - 1; q++) {
						srand(time(0)); //using current time as seed

						Vector r = c + a * rand_float() + b * rand_float();


					};
				};
			}
			else
			{
				
			};*/


			// HARD SHADOWS
			
			Vector lightsource = scene->getLight(i)->position;
			Vector L = (lightsource - phit);		// unit light vector from hit point to light source
			float shadowDist = L.length();
			L.normalize();
			bool shadow = false;
			float ts;
			float intensity = L.normalize() * normal;
			Vector H = (L + I).normalize();  //Half-Vector from Blinn model
			float NH = normal * H;
			if (NH < 0) NH = 0;
			Ray shadowRay = Ray((phit + bias), L);

			if (intensity > 0)  // intensity of light
			{
				for (int sobj_i = 0; sobj_i < scene->getNumObjects(); sobj_i += 1)  //Looping through all objects to check if there is an intersection
				{
					if (scene->getObject(sobj_i)->intercepts(shadowRay, ts) == true && ts < shadowDist)  //check if ray towards source light is intercepting object
					{
						shadow = true;
						break;
					};
				}
			}
			if (!shadow) //if no object is blocking source light
			{
				// calculate inner product between light vector and normal
				float Kd = scene->getObject(hitIndex)->GetMaterial()->GetDiffuse();
				Color diffuse_color = scene->getObject(hitIndex)->GetMaterial()->GetDiffColor() * Kd * intensity; // Calculate diffuse
				float Ks = scene->getObject(hitIndex)->GetMaterial()->GetSpecular();
				float spec = pow(NH, scene->getObject(hitIndex)->GetMaterial()->GetShine());
				Color specular_color = scene->getObject(hitIndex)->GetMaterial()->GetSpecColor() * Ks * spec;  // Calculate specular
				Color specDiffColor = diffuse_color + specular_color;  //combine colors
				finalColor += specDiffColor;
			}
		}

		if (depth < MAX_DEPTH) 
		{	
			Color reflColor;
			Color refrColor;
			float KR = 1;  // defining as 1 to make sure that 100% reflectance is used in case there is no refractivity

			// reflection
			if (scene->getObject(hitIndex)->GetMaterial()->GetReflection() > 0)   // !! If reflective component is bigger than 0, means it is reflective?
			{
				Vector V = (ray.direction) * (-1);		// math from slides
				Vector rRefl = normal * 2 * (V * normal) - V;

				/*bool outside = rRefl * normal < 0;
				Vector reflectionOrigin = outside ? phit - bias : phit + bias;*/

				Vector reflectionOrigin = phit + bias;
				Ray reflRay = Ray(reflectionOrigin, rRefl.normalize());
				Color reflColor = rayTracing(reflRay, depth + 1, ior_1); //iteration
				reflColor = reflColor * scene->getObject(hitIndex)->GetMaterial()->GetReflection(); //  reduce rColor by the specular reflection coefficient	
				finalColor += reflColor;
			};

			/*float eta = ior_1;
			
			//refraction
			if (scene->getObject(hitIndex)->GetMaterial()->GetTransmittance() > 0.0f)
			{
				Vector v = (ray.direction) * (-1);
				Vector vt = normal * (v * normal) - v;
				float sinI = (vt).length();
				float ior_2 = scene->getObject(hitIndex)->GetMaterial()->GetRefrIndex();
				float sinT = (ior_1 / ior_2) * sinI;
				float cosT = sqrt(1 - (sinT) * (sinT));

				//fresnels
				if (sinT >= 1) { //If total reflection
					float KR = 1;
				}
				else {
					float cosI = sqrt(1 - (sinI) * (sinI));
					float Rp = pow((ior_1 * cosI - ior_2 * cosT) / (ior_1 * cosI + ior_2 * cosT), 2);
					float Rs = pow((ior_1 * cosT - ior_2 * cosI) / (ior_1 * cosT + ior_2 * cosI), 2);
					float KR = (Rs + Rp) * 0.5;
				};

				Vector tr = vt * (1 / vt.length());  // math from slides
				Vector rRefr = tr * sinT + normal * (-1) * cosT;
				Ray refrRay = Ray(phit + bias, rRefr);
				bool inside = rRefr * normal > 0;
				if (!inside)
				{
					float eta = 1.0f / scene->getObject(hitIndex)->GetMaterial()->GetRefrIndex();
				}
				else
				{
					float eta = scene->getObject(hitIndex)->GetMaterial()->GetRefrIndex();
				};
				Color refrColor = rayTracing(refrRay, depth + 1, eta); //tColor = trace(scene, point, tRay direction, depth + 1);   not sure about depth +1
				refrColor = refrColor * scene->getObject(hitIndex)->GetMaterial()->GetTransmittance();// reduce tColor by the transmittance coefficient
			};
			*/
			

			// refraction (som på scratch)
			if (scene->getObject(hitIndex)->GetMaterial()->GetTransmittance() > 0.0f)
			{
			
				Vector refrDir;

				float cosi = - std::fmax(-1.f, std::fmin(1.f, ray.direction * normal));  // minus from other example, what does it mean?
				float etai = 1, etat = ior_1;
				Vector n = normal;

				if (cosi < 0) { // if the ray is inside the object, swap the indices and invert the normal to get the correct result
					cosi = -cosi;
					std::swap(etai, etat); n = normal * (-1);
				}

				// Snells law
				float eta = etai / etat;
				float k = 1 - eta * eta * (1 - cosi * cosi);

				if (k < 0) {
					Vector refrDir = Vector(1, 0, 0);
				}
				else {
					Vector refrDir = ray.direction * eta + n * (eta * cosi - sqrtf(k));
				}

				//  Fresnels
				/*Vector v = (ray.direction) * (-1);
				Vector vt = normal * (v * normal) - v;
				float sini = (vt).length();
				float sint = (etai / etat) * sini;
				float cost = sqrt(1 - (sint) * (sint));

				if (sint >= 1) { //If total reflection
					float KR = 1;
				}
				else {
					float cosI = sqrt(1 - (sini) * (sini));
					float Rp = pow((etai * cosi - etat * cost) / (etai * cosi + etat * cost), 2);
					float Rs = pow((etai * cost - etat * cosi) / (etai * cost + etat * cosi), 2);
					float KR = (Rs + Rp) * 0.5;
				};
				*/

				bool outside = refrDir * normal < 0;
				Vector refrOrig = outside ? phit - bias : phit + bias;
				/*
				if (outside)
				{
					float eta = 1.0f / scene->getObject(hitIndex)->GetMaterial()->GetRefrIndex();
				}
				else
				{
					float eta = scene->getObject(hitIndex)->GetMaterial()->GetRefrIndex();
				};
				*/
				Ray refractiveRay = Ray(refrOrig, refrDir);
				refrColor = rayTracing(refractiveRay, depth + 1, scene->getObject(hitIndex)->GetMaterial()->GetRefrIndex()) * scene->getObject(hitIndex)->GetMaterial()->GetTransmittance();

				finalColor += refrColor;

				//finalColor += refrColor * (1 - KR) + reflColor * KR;
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

	/* Só se faz a alocação dos arrays glBufferData (NULL), e o envio dos pontos para a placa gráfica
	é feito na drawPoints com GlBufferSubData em tempo de execução pois os arrays são GL_DYNAMIC_DRAW */
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

	bool antialiasing;
	string antialiasing_answer;
	/*cout << "Do you want to use antialiasing (y/n): ";
	cin >> antialiasing_answer;*/
	if (antialiasing_answer == "y") antialiasing = true;
	else antialiasing = false;

	for (int y = 0; y < RES_Y; y++)
	{
		for (int x = 0; x < RES_X; x++)
		{
			Color color;
			Color cxy;
			int n = 5; // defined by us to decide how much to split the pixel in

			Vector pixel;  //viewport coordinates
			pixel.x = x + 0.5f;
			pixel.y = y + 0.5f;

			// Antialiasing with jittering ... combined solution from book and https://www.scratchapixel.com/code.php?id=13&origin=/lessons/3d-basic-rendering/introduction-to-shading
			if (antialiasing == true) {
				Color c = Color(0.0, 0.0, 0.0);
				for (int p = 0; p <= n - 1; p++) {
					for (int q = 0; q <= n - 1; q++) {
						Vector pixel;  //viewport coordinates
						srand(time(0)); //using current time as seed
						pixel.x = x + 0.5f * ((p + rand_float()) / n);
						pixel.y = y + 0.5f * ((q + rand_float()) / n);
						Ray ray = scene->GetCamera()->PrimaryRay(pixel);
						c = c + rayTracing(ray, 1, 1.0);
					};
				};
				float nsqr = 1 / pow(n, 2);
				color = c * nsqr;
			}

			else {
				//YOUR 2 FUNCTIONS:
				Ray ray = scene->GetCamera()->PrimaryRay(pixel);
				color = rayTracing(ray, 1, 1.0);

			};

			

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
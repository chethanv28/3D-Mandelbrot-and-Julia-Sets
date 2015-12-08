/* CS580 Homework 3 */

#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"
#define _USE_MATH_DEFINES
#include <cmath>

// --------- quaternion representation -----------------------------------------
//
// Each quaternion can be specified by four scalars q = A + Bi + Cj + Dk, so are
// stored as a float4. I’ve tried a struct containing a separate scalar and
// 3-vector to avoid a lot of swizzling, but the float4 representation ends up
// using fewer instructions. A matrix representation is also possible.
//

//---------- QuatMult()----------
//Returns the product of two quaternion number q1 and q2
//Need to change this function to create different Mandelbrot 3D fractals
float4 quatMult(float4 q1, float4 q2){
	float4 r;

	r.x = q1.x*q2.x - q1.y*q2.y - q1.z*q2.z - q1.w*q2.w;
	r.y = q1.x*q2.y + q1.y*q2.x + q1.z*q2.w - q1.w*q2.z;
	r.z = q1.x*q2.z - q1.y*q2.w + q1.z*q2.x + q1.w*q2.y;
	r.w = q1.x*q2.w + q1.y*q2.z - q1.z*q2.y + q1.w*q2.x;

	return r;

}

//maybe define a QuatSq if it makes sense...
//again this might need to change for different fractals
float4 quatSq(float4 q1){
	float4 r;

	r.x = q1.x*q1.x - q1.y*q1.y - q1.z*q1.z - q1.w*q1.w;
	r.y = 2.0*q1.x*q1.y;
	r.z = 2.0*q1.x*q1.z;
	r.w = 2.0*q1.x*q1.w;

	return r;
}

// dot product of two float3 vectors
float dot(float3 a, float3 b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

// dot product of two float4 vectors
float dot(float4 a, float4 b) {
	return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

void iteratePoint(float4 q, float4 qp, float4 c, int maxIterations){

	for (int i = 0; i < maxIterations; i++){
		qp = quatMult(q, qp);
		qp = qp * 2.0;
		q = quatSq(q) + c;

		if (dot(q, q) > ESCAPE_THRESHOLD){
			break;
		}

	}
}

float length(float4 q){
	float l = q.x*q.x + q.y*q.y + q.z*q.z + q.w*q.w;
	return l;
}

float3 normalizef3(float3 f){
	float norm = f.x*f.x + f.y*f.y + f.z*f.z;
	float3 normalizedf;
	normalizedf.x /= norm;
	normalizedf.y /= norm;
	normalizedf.z /= norm;

	return normalizedf;
}

float3 estimateNorm(float3 p, float4 c, int maxIterations){
	float3 N;
	float4 qP = float4(p, 0);
	float gradX, gradY, gradZ;

	float4 gx1 = qP - float4(DEL, 0, 0, 0);
	float4 gx2 = qP + float4(DEL, 0, 0, 0);
	float4 gy1 = qP - float4(0, DEL, 0, 0);
	float4 gy2 = qP + float4(0, DEL, 0, 0);
	float4 gz1 = qP - float4(0, 0, DEL, 0);
	float4 gz2 = qP + float4(0, 0, DEL, 0);

	for (int i = 0; i < maxIterations; i++){
		gx1 = quatSq(gx1) + c;
		gx2 = quatSq(gx2) + c;
		gy1 = quatSq(gy1) + c;
		gy2 = quatSq(gy2) + c;
		gz1 = quatSq(gz1) + c;
		gz2 = quatSq(gz2) + c;
	}

	gradX = length(gx2) - length(gx1);
	gradY = length(gy2) - length(gy1);
	gradZ = length(gz2) - length(gz1);

	N = normalizef3(float3(gradX, gradY, gradZ));

	return N;
}

float intersectQJulia(float3 r0, float3 rD, float4 c, int maxIterations, float epsilon){
	float dist;

	while (1){
		float4 z = float4(r0, 0);

		float4 zp = float4(1, 0, 0, 0);

		iteratePoint(z, zp, c, maxIterations);

		float normZ = length(z);
		dist = 0.5*normZ*log(normZ) / length(zp);

		r0 = r0 + rD*dist;

		//ending condition
		//We need to break if the ray distance is smaller than epsilon or length is past bound
		if (dist <epsilon || dot(r0, r0) > BOUNDING_RADIUS_2)
			break;
	}

	return dist;
}

float3 absf3(float3 n){
	if (n.x < 0)
		n.x*-1;
	if (n.y < 0)
		n.y*-1;
	if (n.z < 0)
		n.z*-1;

	return n;
}

float3 Phong(float3 light, float3 eye, float3 pt, float3 N)
{
	float3 diffuse = float3(1.00, 0.45, 0.25); // base color of shading
	const int specularExponent = 10; // shininess of shading
	const float specularity = 0.45; // amplitude of specular highlight
	float3 L = normalizef3(light - pt); // find the vector to the light
	float3 E = normalizef3(eye - pt); // find the vector to the eye
	float NdotL = dot(N, L); // find the cosine of the angle between light and normal
	float3 R = L - N *2.0 * NdotL; // find the reflected vector
	diffuse = diffuse + absf3(N)*0.3; // add some of the normal to the
	// color to make it more interesting
	// compute the illumination using the Phong equation
	return diffuse * max(NdotL, 0) + specularity*pow(max(dot(E, R), 0), specularExponent);
}

float3 intersectSphere(float3 rO, float3 rD)
{
	float B, C, d, t0, t1, t;
	B = 2 * dot(rO, rD);
	C = dot(rO, rO) - BOUNDING_RADIUS_2;
	d = sqrt(B*B - 4.0 * C);
	t0 = (-B + d) * 0.5;
	t1 = (-B - d) * 0.5;
	t = min(t0, t1);
	rO = rO + rD * t;
	return rO;
}

float4 pixelColor(float3 rO, // ray origin
	float3 rD, // ray direction (unit length)
	float4 mu, // quaternion constant specifying the particular set
	float epsilon, // specifies precision of intersection
	float3 eye, // location of the viewer
	float3 light, // location of a single point light
	bool renderShadows, // flag for turning self-shadowing on/off
	int maxIterations){

	const float4 backgroundColor = float4(0.3, 0.3, 0.3, 0); //define the background color of the image
	float4 color; // This color is the final output of our program.

	color = backgroundColor;
	rD = normalizef3(rD); //the ray direction is interpolated and may need to be normalized
	rO = intersectSphere(rO, rD);
	// Next, try to find a point along the ray which intersects the Julia set.
	float dist = intersectQJulia(rO, rD, mu, maxIterations, epsilon);

	if (dist < epsilon)
	{
		// Determine a "surface normal" which we’ll use for lighting calculations.
		float3 N = estimateNorm(rO, mu, maxIterations);
		// Compute the Phong illumination at the point of intersection.
		float3 phong = Phong(light, rD, rO, N);
		color.x = phong.x; color.y = phong.y; color.z = phong.z;
		color.w = 1; // (make this fragment opaque)
		// If the shadow flag is on, determine if this point is in shadow
		if (renderShadows == true)
		{
			float3 L = normalizef3(light - rO);
			rO = rO + N*epsilon*2.0;
			dist = intersectQJulia(rO, L, mu, maxIterations, epsilon);

			/*if (dist < epsilon){
			color.x = color.x * 0.4; // (darkening the shaded value is not really correct, but looks good)
			color.y = color.y * 0.4;
			color.z = color.z * 0.4;
			}*/
		}
	}
	// Return the final color which is still the background color if we didn’t hit anything.
	return color;
}

int GzRotXMat(float degree, GzMatrix mat)
{
	// Create rotate matrix : rotate along x axis
	// Pass back the matrix using mat value

	if (mat == NULL){
		return GZ_FAILURE;
	}

	//GzMatrix rotx;
	float radian = degree*M_PI / 180;
	//set matrix rows
	mat[0][0] = 1;
	mat[0][1] = 0;
	mat[0][2] = 0;
	mat[0][3] = 0;
	mat[1][0] = 0;
	mat[1][1] = cos(radian);
	mat[1][2] = -sin(radian);
	mat[1][3] = 0;
	mat[2][0] = 0;
	mat[2][1] = sin(radian);
	mat[2][2] = cos(radian);
	mat[2][3] = 0;
	mat[3][0] = 0;
	mat[3][1] = 0;
	mat[3][2] = 0;
	mat[3][3] = 1;

	return GZ_SUCCESS;
}

int GzRotYMat(float degree, GzMatrix mat)
{
	// Create rotate matrix : rotate along y axis
	// Pass back the matrix using mat value
	if (mat == NULL){
		return GZ_FAILURE;
	}

	float radian = degree*M_PI / 180;
	//set matrix rows
	mat[0][0] = cos(radian);
	mat[0][1] = 0;
	mat[0][2] = sin(radian);
	mat[0][3] = 0;
	mat[1][0] = 0;
	mat[1][1] = 1;
	mat[1][2] = 0;
	mat[1][3] = 0;
	mat[2][0] = -sin(radian);
	mat[2][1] = 0;
	mat[2][2] = cos(radian);
	mat[2][3] = 0;
	mat[3][0] = 0;
	mat[3][1] = 0;
	mat[3][2] = 0;
	mat[3][3] = 1;

	return GZ_SUCCESS;
}

int GzRotZMat(float degree, GzMatrix mat)
{
	// Create rotate matrix : rotate along z axis
	// Pass back the matrix using mat value
	if (mat == NULL){
		return GZ_FAILURE;
	}

	float radian = degree*M_PI / 180;
	//set matrix rows
	mat[0][0] = cos(radian);
	mat[0][1] = -sin(radian);
	mat[0][2] = 0;
	mat[0][3] = 0;
	mat[1][0] = sin(radian);
	mat[1][1] = cos(radian);
	mat[1][2] = 0;
	mat[1][3] = 0;
	mat[2][0] = 0;
	mat[2][1] = 0;
	mat[2][2] = 1;
	mat[2][3] = 0;
	mat[3][0] = 0;
	mat[3][1] = 0;
	mat[3][2] = 0;
	mat[3][3] = 1;
	return GZ_SUCCESS;
}

int GzTrxMat(GzCoord translate, GzMatrix mat)
{
	// Create translation matrix
	// Pass back the matrix using mat value
	if (mat == NULL){
		return GZ_FAILURE;
	}

	//set matrix rows
	mat[0][0] = 1;
	mat[0][1] = 0;
	mat[0][2] = 0;
	mat[0][3] = translate[X];
	mat[1][0] = 0;
	mat[1][1] = 1;
	mat[1][2] = 0;
	mat[1][3] = translate[Y];
	mat[2][0] = 0;
	mat[2][1] = 0;
	mat[2][2] = 1;
	mat[2][3] = translate[Z];
	mat[3][0] = 0;
	mat[3][1] = 0;
	mat[3][2] = 0;
	mat[3][3] = 1;

	return GZ_SUCCESS;
}

int GzScaleMat(GzCoord scale, GzMatrix mat)
{
	// Create scaling matrix
	// Pass back the matrix using mat value
	if (mat == NULL){
		return GZ_FAILURE;
	}

	//set matrix rows
	mat[0][0] = scale[X];
	mat[0][1] = 0;
	mat[0][2] = 0;
	mat[0][3] = 0;
	mat[1][0] = 0;
	mat[1][1] = scale[Y];
	mat[1][2] = 0;
	mat[1][3] = 0;
	mat[2][0] = 0;
	mat[2][1] = 0;
	mat[2][2] = scale[Z];
	mat[2][3] = 0;
	mat[3][0] = 0;
	mat[3][1] = 0;
	mat[3][2] = 0;
	mat[3][3] = 1;

	return GZ_SUCCESS;
}

/* Given vector (GzCoord) as input we find the normalization factor of it*/
int normalizeVector(GzCoord g){
	if (g == NULL)
		return GZ_FAILURE;

	float x2 = g[X] * g[X];
	float y2 = g[Y] * g[Y];
	float z2 = g[Z] * g[Z];

	float normal = sqrt(x2 + y2 + z2);

	GzCoord n;
	g[X] = g[X] / normal;
	g[Y] = g[Y] / normal;
	g[Z] = g[Z] / normal;

	return GZ_SUCCESS;
}

/* Return dot product of two vectors*/
float dotProduct(GzCoord g, GzCoord g1){
	float dotprod = g[X] * g1[X] + g[Y] * g1[Y] + g[Z] * g1[Z];

	return dotprod;
}

int normalizeMatrix(GzMatrix g){
	float K = 1 / (sqrt(g[0][0] * g[0][0] + g[0][1] * g[0][1] + g[0][2] * g[0][2] + g[0][3] * g[0][3]));

	for (int i = 0; i < 4; i++){
		for (int j = 0; j < 4; j++){
			g[i][j] *= K;
		}
	}

	return GZ_SUCCESS;
}

/* Given normal verts for shading we need to calculate color */
int GzShader(GzRender* render, GzCoord normal, GzColor color){
	if (render == NULL)
		return GZ_FAILURE;

	//normal should be normalized but just in case
	normalizeVector(normal);

	//set eye vector
	GzCoord eye;
	eye[X] = 0.0;
	eye[Y] = 0.0;
	eye[Z] = -1.0;

	normalizeVector(eye);

	//compute reflections for each light
	GzCoord* reflection = new GzCoord[render->numlights];
	int* normDirection = new int[render->numlights];

	float NL;
	float NE;

	NE = dotProduct(normal, eye);
	for (int i = 0; i < render->numlights; i++){
		NL = dotProduct(normal, render->lights[i].direction);
		//regular case
		if (NE > 0 && NL > 0){
			reflection[i][X] = 2 * NL*normal[X] - render->lights[i].direction[X];
			reflection[i][Y] = 2 * NL*normal[Y] - render->lights[i].direction[Y];
			reflection[i][Z] = 2 * NL*normal[Z] - render->lights[i].direction[Z];
			normalizeVector(reflection[i]);
			normDirection[i] = 1;
		}
		//flip normal
		else if (NE < 0 && NL < 0){
			reflection[i][X] = 2 * NL*(-normal[X]) - render->lights[i].direction[X];
			reflection[i][Y] = 2 * NL*(-normal[Y]) - render->lights[i].direction[Y];
			reflection[i][Z] = 2 * NL*(-normal[Z]) - render->lights[i].direction[Z];
			normalizeVector(reflection[i]);
			normDirection[i] = -1;
		}
		//do nothing
		else{
			normDirection[i] = 0;
			continue;
		}
	}

	//specular coefficient
	//Ks * sum[le (R dot E)^s]
	GzColor specularSum;
	GzColor specular;
	specularSum[0] = 0; specularSum[1] = 0; specularSum[2] = 0;
	for (int i = 0; i < render->numlights; i++){
		if (normDirection[i] == 0)
			continue;

		float RE = dotProduct(reflection[i], eye);
		//clip R dot E to keep in range [0,1]
		if (RE < 0){
			RE = 0;
		}
		else if (RE > 1){
			RE = 1;
		}

		specularSum[0] += render->lights[i].color[0] * pow(RE, render->spec);
		specularSum[1] += render->lights[i].color[1] * pow(RE, render->spec);
		specularSum[2] += render->lights[i].color[2] * pow(RE, render->spec);
	}
	specular[0] = render->Ks[0] * specularSum[0];
	specular[1] = render->Ks[1] * specularSum[1];
	specular[2] = render->Ks[2] * specularSum[2];

	//diffusion coefficient
	//Kd sum[le N dot L]
	GzColor diffusionSum = { 0, 0, 0 };
	GzColor diffusion = { 0, 0, 0 };
	for (int i = 0; i < render->numlights; i++){
		if (normDirection[i] == 0){
			continue;
		}
		//negate normal
		if (normDirection[i] == -1){
			GzCoord negativeNorm;
			negativeNorm[X] = -normal[X];
			negativeNorm[Y] = -normal[Y];
			negativeNorm[Z] = -normal[Z];
			float negNL = dotProduct(negativeNorm, render->lights[i].direction);
			diffusionSum[0] += render->lights[i].color[0] * negNL;
			diffusionSum[1] += render->lights[i].color[1] * negNL;
			diffusionSum[2] += render->lights[i].color[2] * negNL;
		}
		//regular normal
		else if (normDirection[i] == 1){
			float NL1 = dotProduct(normal, render->lights[i].direction);
			diffusionSum[0] += render->lights[i].color[0] * NL1;
			diffusionSum[1] += render->lights[i].color[1] * NL1;
			diffusionSum[2] += render->lights[i].color[2] * NL1;
		}

		if (true)
			int y = 10;

	}
	diffusion[0] = render->Kd[0] * diffusionSum[0];
	diffusion[1] = render->Kd[1] * diffusionSum[1];
	diffusion[2] = render->Kd[2] * diffusionSum[2];

	//compute ambient
	//Ka * la
	GzColor ambient = { 0, 0, 0 };
	ambient[0] = render->Ka[0] * render->ambientlight.color[0];
	ambient[1] = render->Ka[1] * render->ambientlight.color[1];
	ambient[2] = render->Ka[2] * render->ambientlight.color[2];

	//compute color
	color[0] = specular[0] + diffusion[0] + ambient[0];
	color[1] = specular[1] + diffusion[1] + ambient[1];
	color[2] = specular[2] + diffusion[2] + ambient[2];

	if (true)
		int k = 10;

	return GZ_SUCCESS;
}

//----------------------------------------------------------
// Begin main functions

int GzNewRender(GzRender **render, GzDisplay	*display)
{
	/*
	- malloc a renderer struct
	- setup Xsp and anything only done once
	- save the pointer to display
	- init default camera
	*/

	*render = new GzRender[1];
	render[0]->display = display;

	//intialize Xsp array
	render[0]->Xsp[0][0] = display->xres / 2;
	render[0]->Xsp[0][1] = 0; render[0]->Xsp[0][2] = 0;
	render[0]->Xsp[0][3] = display->xres / 2;
	render[0]->Xsp[1][1] = -(display->yres / 2);
	render[0]->Xsp[1][0] = 0; render[0]->Xsp[1][2] = 0;
	render[0]->Xsp[1][3] = display->yres / 2;
	render[0]->Xsp[2][2] = MAXINT;
	render[0]->Xsp[2][0] = 0; render[0]->Xsp[2][1] = 0; render[0]->Xsp[2][3] = 0;
	render[0]->Xsp[3][3] = 1;
	render[0]->Xsp[3][0] = 0; render[0]->Xsp[3][1] = 0; render[0]->Xsp[3][2] = 0;


	//set my camera stuff yo
	render[0]->camera.position[X] = DEFAULT_IM_X;
	render[0]->camera.position[Y] = DEFAULT_IM_Y;
	render[0]->camera.position[Z] = DEFAULT_IM_Z;

	render[0]->camera.lookat[X] = 0.0;
	render[0]->camera.lookat[Y] = 0.0;
	render[0]->camera.lookat[Z] = 0.0;

	render[0]->camera.worldup[X] = 0.0;
	render[0]->camera.worldup[Y] = 1.0;
	render[0]->camera.worldup[Z] = 0.0;

	render[0]->camera.FOV = DEFAULT_FOV;


	//this render is so empty man
	render[0]->matlevel = -1;
	render[0]->numlights = 0;

	return GZ_SUCCESS;

}

int GzFreeRender(GzRender *render)
{
	/*
	-free all renderer resources
	*/
	//GzFreeDisplay(render->display);
	//delete render->camera.lookat;
	//delete render->camera.position;
	//delete render->camera.worldup;
	//delete render->camera.Xiw;
	//delete render->camera.Xpi;
	//delete render->Xsp;
	//delete render->Ximage;
	delete render;


	return GZ_SUCCESS;
}

int GzBeginRender(GzRender *render)
{
	/*
	- setup for start of each frame - init frame buffer color,alpha,z
	- compute Xiw and projection xform Xpi from camera definition
	- init Ximage - put Xsp at base of stack, push on Xpi and Xiw
	- now stack contains Xsw and app can push model Xforms when needed
	*/

	GzInitDisplay(render->display);
	//set Xpi
	render->camera.Xpi[0][0] = 1;
	render->camera.Xpi[0][1] = 0; render->camera.Xpi[0][2] = 0; render->camera.Xpi[0][3] = 0;
	render->camera.Xpi[1][1] = 1;
	render->camera.Xpi[1][0] = 0; render->camera.Xpi[1][2] = 0; render->camera.Xpi[1][3] = 0;
	render->camera.Xpi[2][2] = tan((render->camera.FOV / 2.0)*M_PI / 180.0);
	render->camera.Xpi[2][0] = 0; render->camera.Xpi[2][1] = 0; render->camera.Xpi[2][3] = 0;
	render->camera.Xpi[3][0] = 0; render->camera.Xpi[3][1] = 0;
	render->camera.Xpi[3][2] = tan((render->camera.FOV / 2.0)*M_PI / 180.0);
	render->camera.Xpi[3][3] = 1;

	//set Xiw
	//need to calculate Z
	GzCoord cz;
	cz[X] = render->camera.lookat[X] - render->camera.position[X];
	cz[Y] = render->camera.lookat[Y] - render->camera.position[Y];
	cz[Z] = render->camera.lookat[Z] - render->camera.position[Z];
	normalizeVector(cz);

	//need to calculate Y
	GzCoord cy, upprime;
	float dot = dotProduct(render->camera.worldup, cz);
	upprime[X] = render->camera.worldup[X] - dot*cz[X];
	upprime[Y] = render->camera.worldup[Y] - dot*cz[Y];
	upprime[Z] = render->camera.worldup[Z] - dot*cz[Z];
	cy[X] = upprime[X];
	cy[Y] = upprime[Y];
	cy[Z] = upprime[Z];
	normalizeVector(cy);

	//need to calculate X bro
	//A = a2*b3 - a3*b2
	//B = a3*b1 - a1*b3
	//C = a1*b2 - a2*b1
	GzCoord cx;
	cx[X] = cy[Y] * cz[Z] - cy[Z] * cz[Y];
	cx[Y] = cy[Z] * cz[X] - cy[X] * cz[Z];
	cx[Z] = cy[X] * cz[Y] - cy[Y] * cz[X];
	normalizeVector(cx);


	render->camera.Xiw[0][0] = cx[X];
	render->camera.Xiw[0][1] = cx[Y];
	render->camera.Xiw[0][2] = cx[Z];
	render->camera.Xiw[0][3] = -dotProduct(cx, render->camera.position);
	render->camera.Xiw[1][0] = cy[X];
	render->camera.Xiw[1][1] = cy[Y];
	render->camera.Xiw[1][2] = cy[Z];
	render->camera.Xiw[1][3] = -dotProduct(cy, render->camera.position);
	render->camera.Xiw[2][0] = cz[X];
	render->camera.Xiw[2][1] = cz[Y];
	render->camera.Xiw[2][2] = cz[Z];
	render->camera.Xiw[2][3] = -dotProduct(cz, render->camera.position);
	render->camera.Xiw[3][0] = 0; render->camera.Xiw[3][1] = 0; render->camera.Xiw[3][2] = 0;
	render->camera.Xiw[3][3] = 1;

	GzMatrix I = { 1.0, 0.0, 0.0, 0.0,
		0.0, 1.0, 0.0, 0.0,
		0.0, 0.0, 1.0, 0.0,
		0.0, 0.0, 0.0, 1.0 };



	//Need to push Xsp => Xpi => Xiw
	GzPushMatrix(render, render->Xsp);
	GzPushNormalMatrix(render, I);

	GzPushMatrix(render, render->camera.Xpi);
	GzPushNormalMatrix(render, I);

	GzPushMatrix(render, render->camera.Xiw);
	//GzPushNormalMatrix(render, render->camera.Xiw);

	return GZ_SUCCESS;
}

int GzPutCamera(GzRender *render, GzCamera *camera)
{
	/*
	- overwrite renderer camera structure with new camera definition
	*/
	render->camera.FOV = camera->FOV;
	render->camera.lookat[X] = camera->lookat[X];
	render->camera.lookat[Y] = camera->lookat[Y];
	render->camera.lookat[Z] = camera->lookat[Z];
	render->camera.position[X] = camera->position[X];
	render->camera.position[Y] = camera->position[Y];
	render->camera.position[Z] = camera->position[Z];

	render->camera.worldup[X] = camera->worldup[X];
	render->camera.worldup[Y] = camera->worldup[Y];
	render->camera.worldup[Z] = camera->worldup[Z];

	return GZ_SUCCESS;
}

int GzPushNormalMatrix(GzRender *render, GzMatrix matrix){
	/*
	- push a matrix onto the Xnorm stack
	- check for stack overflow
	*/
	if (render == NULL)
		return GZ_FAILURE;
	if (matrix == NULL)
		return GZ_FAILURE;

	//we have overflow
	//render->matlevel = render->matlevel + 1;
	if (render->matlevel >= MATLEVELS)
		return GZ_FAILURE;

	//GzMatrix g;
	//get Rid of T, but should we save it? No. Cause...yeah
	//matrix[0][3] = 0;
	//matrix[1][3] = 0;
	//matrix[2][3] = 0;

	normalizeMatrix(matrix);

	matrix[3][0] = 0;

	//push that Matrix yo
	if (render->matlevel == 0){
		for (int i = 0; i < 4; i++){
			for (int j = 0; j < 4; j++){
				render->Xnorm[render->matlevel][i][j] = matrix[i][j];
			}
		}
	}
	else{
		//multiply the matrix and then we need to put it at the top
		for (int i = 0; i < 4; i++){
			for (int j = 0; j < 4; j++){
				render->Xnorm[render->matlevel][i][j] = render->Xnorm[render->matlevel - 1][i][0] * matrix[0][j] + render->Xnorm[render->matlevel - 1][i][1] * matrix[1][j] +
					render->Xnorm[render->matlevel - 1][i][2] * matrix[2][j] + render->Xnorm[render->matlevel - 1][i][3] * matrix[3][j];
			}
		}

	}

	if (true)
		int k = 10;

	return GZ_SUCCESS;
}

int GzPushMatrix(GzRender *render, GzMatrix	matrix)
{
	/*
	- push a matrix onto the Ximage stack
	- check for stack overflow
	*/
	if (render == NULL)
		return GZ_FAILURE;
	if (matrix == NULL)
		return GZ_FAILURE;


	//we have overflow
	render->matlevel = render->matlevel + 1;
	if (render->matlevel >= MATLEVELS)
		return GZ_FAILURE;

	GzMatrix g;
	float multvalue = 0;
	if (render->matlevel == 0){
		for (int i = 0; i < 4; i++){
			for (int j = 0; j < 4; j++){
				render->Ximage[render->matlevel][i][j] = matrix[i][j];
			}
		}
	}
	else{
		//multiply the matrix and then we need to put it at the top
		for (int i = 0; i < 4; i++){
			for (int j = 0; j < 4; j++){
				render->Ximage[render->matlevel][i][j] = render->Ximage[render->matlevel - 1][i][0] * matrix[0][j] + render->Ximage[render->matlevel - 1][i][1] * matrix[1][j] +
					render->Ximage[render->matlevel - 1][i][2] * matrix[2][j] + render->Ximage[render->matlevel - 1][i][3] * matrix[3][j];
			}
		}

	}

	GzPushNormalMatrix(render, matrix);

	return GZ_SUCCESS;
}

int GzPopMatrix(GzRender *render)
{
	/*
	- pop a matrix off the Ximage stack
	- check for stack underflow
	*/
	if (render == NULL)
		return GZ_FAILURE;

	//underflow
	if (render->Ximage[0] == NULL)
		return GZ_FAILURE;

	if (render->matlevel < 0)
		return GZ_FAILURE;

	//change top of pointer to one level lower
	render->matlevel--;


	return GZ_SUCCESS;
}

int GzPutAttribute(GzRender	*render, int numAttributes, GzToken	*nameList, GzPointer	*valueList) /* void** valuelist */
{
	/*
	- set renderer attribute states (e.g.: GZ_RGB_COLOR default color)
	- later set shaders, interpolaters, texture maps, and lights
	*/
	//error checking
	if (render == NULL)
		return GZ_FAILURE;
	if (nameList == NULL)
		return GZ_FAILURE;
	if (valueList == NULL)
		return GZ_FAILURE;

	for (int i = 0; i < numAttributes; i++){
		if (nameList[i] == GZ_RGB_COLOR){
			GzColor* a = (GzColor *)(valueList[i]);
			float r = a[0][0];
			float g = a[0][1];
			float b = a[0][2];

			render->flatcolor[0] = r;
			render->flatcolor[1] = g;
			render->flatcolor[2] = b;
		}
		else if (nameList[i] == GZ_AMBIENT_COEFFICIENT){
			GzColor* a = (GzColor *)valueList[i];

			render->Ka[0] = a[0][0];
			render->Ka[1] = a[0][1];
			render->Ka[2] = a[0][2];
		}
		else if (nameList[i] == GZ_AMBIENT_LIGHT){
			GzLight* l = (GzLight*)(valueList[i]);

			render->ambientlight = *l;
		}
		else if (nameList[i] == GZ_DIFFUSE_COEFFICIENT){
			GzColor* d = (GzColor *)valueList[i];

			render->Kd[0] = d[0][0];
			render->Kd[1] = d[0][1];
			render->Kd[2] = d[0][2];
		}
		else if (nameList[i] == GZ_DIRECTIONAL_LIGHT){
			GzLight* l = (GzLight*)(valueList[i]);

			render->lights[render->numlights] = *l;
			render->numlights += 1;
			if (true)
				int k = 10;
		}
		else if (nameList[i] == GZ_DISTRIBUTION_COEFFICIENT){
			float* n = (float*)valueList[i];

			render->spec = *n;
		}
		else if (nameList[i] == GZ_SPECULAR_COEFFICIENT){
			GzColor* s = (GzColor *)valueList[i];

			render->Ks[0] = s[0][0];
			render->Ks[1] = s[0][1];
			render->Ks[2] = s[0][2];
		}
		else if (nameList[i] == GZ_INTERPOLATE){
			int* n = (int *)valueList[i];

			render->interp_mode = *n;
		}
	}

	return GZ_SUCCESS;
}

int GzPutRays(GzRender	*render)
{
	/*
	- pass in a triangle description with tokens and values corresponding to
	GZ_POSITION:3 vert positions in model space
	- Xform positions of verts using matrix on top of stack
	- Clip - just discard any triangle with any vert(s) behind view plane
	- optional: test for triangles with all three verts off-screen (trivial frustum cull)
	- invoke triangle rasterizer
	*/

	//declare variables to store vertex information

	float xgap = (4.0) / render->display->xres;//(secondxdim - firstxdim) / render->display->xres;
	float ygap = (6.0) / render->display->yres;//(secondydim - firstydim) / render->display->yres;

	float4* colorValue = new float4[render->display->xres*render->display->yres];
	float* normRayOriginX = new float[render->display->xres*render->display->yres];
	float* normRayOriginY = new float[render->display->xres*render->display->yres];
	float* normRayOriginZ = new float[render->display->xres*render->display->yres];

	GzMatrix topXform;
	bool isBehindPlane = false;

	float4 rorigin;
	float4 modelrorigin;
	float W;
	//variables needed for coloring algorithm
	float4 mu = (0.20, 0.20, 0.0, 0.0);
	float3 normrayorigin;
	float3 raydir;
	float epsilon = 0.003;
	float3 light = ( 0.3, 0.3, 0.3 );
	float3 eye = { render->camera.lookat[X], render->camera.lookat[Y], render->camera.lookat[Z] };
	int xstep = 0;
	int ystep = 0;

	float firstxdim = -1.5;
	float secondxdim = 1.5;
	float firstydim = -2.5;
	float secondydim = 2.5;
	//float firstzdim = -;
	float secondzdim = 3.0;

	for (int m = 0; m < 4; m++){
		for (int n = 0; n < 4; n++){
			topXform[m][n] = render->Ximage[render->matlevel][m][n];
		}
	}

	//need to convert points in interval
	for (float i = firstxdim; i < secondxdim; i += xgap){
		xstep++;
		for (float j = firstydim; j < secondydim; j += ygap){
			ystep++;
			rorigin = float4(i, j, 2.0, 0.0);

			modelrorigin.x = rorigin.x * topXform[0][0] + rorigin.y * topXform[0][1] + rorigin.z * topXform[0][2] + topXform[0][3];
			modelrorigin.y = rorigin.x * topXform[1][0] + rorigin.y * topXform[1][1] + rorigin.z * topXform[1][2] + topXform[1][3];
			modelrorigin.z = rorigin.x * topXform[2][0] + rorigin.y * topXform[2][1] + rorigin.z * topXform[2][2] + topXform[2][3];
			W = rorigin.x * topXform[3][0] + rorigin.y * topXform[3][1] + rorigin.z * topXform[3][2] + topXform[3][3];


			//check to see if we are behind of the viewing plane
			if (modelrorigin.z < render->camera.position[Z]){
				isBehindPlane = true;
				break;
			}

			//normalize to 3d
			modelrorigin.x /= W;
			modelrorigin.y /= W;
			modelrorigin.z /= W;

			//normalized ray origin
			normrayorigin = (modelrorigin.x, modelrorigin.y, modelrorigin.z);
			//how do we calculate ray direction? 
			//Should it be to the camera look-at point translated to our current respective position from the camera's origin?
			raydir = (0.9, 0.5, 0.0);

			//We need to save where these rays are emanating from and their respective colors
			//Then we can use these for something but I'm not sure what
			float4 colorrayval = pixelColor(normrayorigin, raydir, mu, epsilon, eye, light, false, 100);
			colorValue[ystep + xstep*render->display->yres] = colorrayval;
			normRayOriginX[ystep + xstep*render->display->yres] = normrayorigin.x;
			normRayOriginY[ystep + xstep*render->display->yres] = normrayorigin.y;
			normRayOriginZ[ystep + xstep*render->display->yres] = normrayorigin.z;

		}
		ystep = 0;
	}

	if (true){
		int x = 10;
	}


	//go through all pixels and determine which ones should be rasterized
	//if pixel is in field of view render it
	//start with y and go across each y value
	for (int i = 0; i < render->display->xres; i++){
		for (int j = 0; j < render->display->yres; j++){

			GzIntensity r, g, b, a;
			GzDepth z = 0;

			GzGetDisplay(render->display, j, i, &r, &g, &b, &a, &z);


			//Need to put color found from ray here
			r = (GzIntensity)ctoi(colorValue[j + i*render->display->yres].x);
			g = (GzIntensity)ctoi(colorValue[j + i*render->display->yres].y);
			b = (GzIntensity)ctoi(colorValue[j + i*render->display->yres].z);
			z = normRayOriginZ[j + i*render->display->yres];
			GzPutDisplay(render->display, j, i, r, g, b, a, z);

		}
	}


	return GZ_SUCCESS;
}


/* Find min of three float values
*/
float findMin(float a, float b, float c){
	float minArray[3] = { a, b, c };
	float min;

	for (int i = 0; i < 3; i++){
		if (i == 0)
			min = a;
		else{
			if (minArray[i] < min)
				min = minArray[i];
		}
	}

	return min;
}

/* Find max of three float values
*/
float findMax(float a, float b, float c){
	float maxArray[3] = { a, b, c };
	float max;

	for (int i = 0; i < 3; i++){
		if (i == 0)
			max = a;
		else{
			if (maxArray[i] > max)
				max = maxArray[i];
		}
	}

	return max;
}

/*	Given two edges for a tri calculate the normal vector for the tri
i is the y value for a point given and j is the x value */
//float interpolateZ(Edge* e1, Edge* e2, int i, int j){
float interpolateZ(GzCoord v1, GzCoord v2, GzCoord v3, int i, int j){

	float interz = 0;

	float a1 = v2[X] - v1[X];//e1->dx; //
	float a2 = v2[Y] - v1[Y];//e1->dy; //
	float a3 = v2[Z] - v1[Z];//e1->dz; ////e1->z - e1->dz;
	float b1 = v3[X] - v2[X];//e2->dx; ////e2->x - e2->dx;
	float b2 = v3[Y] - v2[Y];//e2->dy; ////e2->y - e2->dy;
	float b3 = v3[Z] - v2[Z];//e2->dz; ////e2->z - e2->dz;
	float A = a2*b3 - a3*b2;
	float B = a3*b1 - a1*b3;
	float C = a1*b2 - a2*b1;
	//solve for D; Ax+By+Cz+D=0
	float D = -A*v1[X] - B*v1[Y] - C*v1[Z];

	//solve for the interpolated z => -(Ax+By+D)/C = z
	interz = (-A*(float)j - B*(float)i - D) / C;

	return interz;
}


/* NOT part of API - just for general assistance */

short	ctoi(float color)		/* convert float color to GzIntensity short */
{
	return(short)((int)(color * ((1 << 12) - 1)));
}


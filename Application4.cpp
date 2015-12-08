// Application4.cpp: implementation of the Application4 class.
//
//////////////////////////////////////////////////////////////////////

/*
 * application test code for homework assignment 
 */

#include "stdafx.h"
#include "CS580HW.h"
#include "Application4.h"
#include "Gz.h"
#include "disp.h"
#include "rend.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

#define INFILE  "pot4.asc"
//#define INFILE "tri.asc"
#define OUTFILE "output.ppm"


void shade(GzCoord norm, GzCoord color);

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Application4::Application4()
{

}

Application4::~Application4()
{
	Clean();
}

int Application4::Initialize()
{
/* to be filled in by the app if it sets camera params */

	GzCamera	camera;  
	int		xRes, yRes;		/* display parameters */ 

	GzToken		nameListShader[9]; 	/* shader attribute names */
	GzPointer   valueListShader[9];		/* shader attribute pointers */
	GzToken     nameListLights[10];		/* light info */
	GzPointer   valueListLights[10];
	int		shaderType, interpStyle;
	float		specpower;
	int		status; 
 
	status = 0; 

	/* 
	 * Allocate memory for user input
	 */
	m_pUserInput = new GzInput;

	/* 
	 * initialize the display and the renderer 
	 */ 

	m_nWidth = 256;	 	// frame buffer and display width
	m_nHeight = 256;	// frame buffer and display height

	status |= GzNewFrameBuffer(&m_pFrameBuffer, m_nWidth, m_nHeight);

	status |= GzNewDisplay(&m_pDisplay, m_nWidth, m_nHeight);

	status |= GzGetDisplayParams(m_pDisplay, &xRes, &yRes); 
 
	status |= GzNewRender(&m_pRender, m_pDisplay); 

/* Translation matrix */
GzMatrix	scale = 
{ 
	1.0,	0.0,	0.0,	0.0, 
	0.0,	1.0,	0.0,	0.0, 
	0.0,	0.0,	1.0,	0.0, 
	0.0,	0.0,	0.0,	1.0 
}; 
 
GzMatrix	rotateX = 
{ 
	1.0,	0.0,	0.0,	0.0, 
	0.0,	.7071,	.7071,	0.0, 
	0.0,	-.7071,	.7071,	0.0, 
	0.0,	0.0,	0.0,	1.0 
}; 
 
GzMatrix	rotateY = 
{ 
	.866,	0.0,	-0.5,	0.0, 
	0.0,	1.0,	0.0,	0.0, 
	0.5,	0.0,	.866,	0.0, 
	0.0,	0.0,	0.0,	1.0 
}; 

#if 0 	/* set up app-defined camera if desired, else use camera defaults */
	camera.position[X] = 13.2;      
  	camera.position[Y] = -8.7;
  	camera.position[Z] = -14.8;

  	camera.lookat[X] = 0.8;
  	camera.lookat[Y] = 0.7;
  	camera.lookat[Z] = 4.5;

  	camera.worldup[X] = -0.2;
  	camera.worldup[Y] = 1.0;
  	camera.worldup[Z] = 0.0;

	camera.FOV = 53.7;              /* degrees */

	status |= GzPutCamera(m_pRender, &camera); 
#endif 

	/* Start Renderer */
	status |= GzBeginRender(m_pRender);

	/* Light */
	GzLight	light1 = { {-0.7071, 0.7071, 0}, {0.5, 0.5, 0.9} };
	GzLight	light2 = { {0, -0.7071, -0.7071}, {0.9, 0.2, 0.3} };
	GzLight	light3 = { {0.7071, 0.0, -0.7071}, {0.2, 0.7, 0.3} };
	//was {.3, .3, .3}
	GzLight	ambientlight = { {0, 0, 0}, {0.4, 0.4, 0.4} };

	/* Material property */
	GzColor specularCoefficient = { 0.4, 0.4, 0.4 };
	//was {.1, .1, .1}
	GzColor ambientCoefficient = { 0.5, 0.5, 0.5 };
	GzColor diffuseCoefficient = { 0.9, 0.9, 0.9 };

/* 
  renderer is ready for frame --- define lights and shader at start of frame 
*/

        /*
         * Tokens associated with light parameters
         */
        nameListLights[0] = GZ_DIRECTIONAL_LIGHT;
        valueListLights[0] = (GzPointer)&light1;
        nameListLights[1] = GZ_DIRECTIONAL_LIGHT;
        valueListLights[1] = (GzPointer)&light2;
        nameListLights[2] = GZ_DIRECTIONAL_LIGHT;
        valueListLights[2] = (GzPointer)&light3;
        status |= GzPutAttribute(m_pRender, 3, nameListLights, valueListLights);

        nameListLights[0] = GZ_AMBIENT_LIGHT;
        valueListLights[0] = (GzPointer)&ambientlight;
        status |= GzPutAttribute(m_pRender, 1, nameListLights, valueListLights);

        /*
         * Tokens associated with shading 
         */
        nameListShader[0]  = GZ_DIFFUSE_COEFFICIENT;
        valueListShader[0] = (GzPointer)diffuseCoefficient;

	/* 
	* Select either GZ_COLOR or GZ_NORMALS as interpolation mode  
	*/
        nameListShader[1]  = GZ_INTERPOLATE;
#if 1
		//interpStyle = GZ_FLAT;
        interpStyle = GZ_COLOR;         /* Gouraud shading */
#else 
        interpStyle = GZ_NORMALS;       /* Phong shading */
#endif

        valueListShader[1] = (GzPointer)&interpStyle;
        nameListShader[2]  = GZ_AMBIENT_COEFFICIENT;
        valueListShader[2] = (GzPointer)ambientCoefficient;
        nameListShader[3]  = GZ_SPECULAR_COEFFICIENT;
        valueListShader[3] = (GzPointer)specularCoefficient;
        nameListShader[4]  = GZ_DISTRIBUTION_COEFFICIENT;
        specpower = 5;
        valueListShader[4] = (GzPointer)&specpower;

	status |= GzPutAttribute(m_pRender, 5, nameListShader, valueListShader);

	//don't need to push the matrices
	//status |= GzPushMatrix(m_pRender, scale);  
	//status |= GzPushMatrix(m_pRender, rotateY);
	//status |= GzPushMatrix(m_pRender, rotateX);

	if (status) exit(GZ_FAILURE); 

	if (status) 
		return(GZ_FAILURE); 
	else 
		return(GZ_SUCCESS); 
}

int Application4::Render() 
{
	GzToken		nameListTriangle[3]; 	/* vertex attribute names */
	GzPointer	valueListTriangle[3]; 	/* vertex attribute pointers */
	GzToken		nameListColor[3];		/* color type names */
	GzPointer	valueListColor[3];	/* color type rgb pointers */
	GzCoord		vertexList[3];		/* vertex position coordinates */ 
	GzCoord		normalList[3];		/* vertex normals */ 
	GzTextureIndex  uvList[3];		/* vertex texture map indices */ 
	GzColor		color; 
	char		dummy[256]; 
	int		status; 


	/* Initialize Display */
	status |= GzInitDisplay(m_pDisplay); 
	
	 
	// I/O File open
	/*FILE *infile;
	if( (infile  = fopen( INFILE , "r" )) == NULL )
	{
         AfxMessageBox( "The input file was not opened\n" );
		 return GZ_FAILURE;
	}*/

	FILE *outfile;
	if( (outfile  = fopen( OUTFILE , "wb" )) == NULL )
	{
         AfxMessageBox( "The output file was not opened\n" );
		 return GZ_FAILURE;
	}

	/* 
	* Walk through the list of triangles, set color 
	* and render each triangle 
	*/ 
	float4 mu = float4(-.7323, -0.2179, 0.0, 0.0);
	for (int i = 0; i < 1; i++){
		mu.x += .02;
		//colorType: 0 for z, 1 for normal
		//quatType: 0 for quaternion, 1 for hyper, 2 for complex
		//setType: 0 for Julia, 1 for Mandelbrot
		//iterations: number of iterations use for set
		GzPutRays(m_pRender, mu, 20, 0, 0, 0);

		GzFlushDisplay2File(outfile, m_pDisplay, i, mu, 0, 0, 0); 	/* write out or update display to file*/
	}
	GzFlushDisplay2FrameBuffer(m_pFrameBuffer, m_pDisplay);	// write out or update display to frame buffer

	/* 
	 * Close file
	 */ 

	//if( fclose( infile ) )
      //AfxMessageBox( "The input file was not closed\n" );

	if( fclose( outfile ) )
      AfxMessageBox( "The output file was not closed\n" );
 
	if (status) 
		return(GZ_FAILURE); 
	else 
		return(GZ_SUCCESS); 
}

int Application4::Clean()
{
	/* 
	 * Clean up and exit 
	 */ 
	int	status = 0; 

	status |= GzFreeRender(m_pRender); 
	status |= GzFreeDisplay(m_pDisplay);
	
	if (status) 
		return(GZ_FAILURE); 
	else 
		return(GZ_SUCCESS);
}

/* 
This doesn't really belong in the application program, but for this 
simplified case of a renderer that doesn't do any shading itself, this 
is the easiest place to put it.
*/

void shade(GzCoord norm, GzCoord color)
{
  GzCoord	light;
  float		coef;

  light[0] = 0.707f;
  light[1] = 0.5f;
  light[2] = 0.5f;

  coef = light[0]*norm[0] + light[1]*norm[1] + light[2]*norm[2];
  if (coef < 0) 	coef *= -1;

  if (coef > 1.0)	coef = 1.0;
  color[0] = coef*0.95f;
  color[1] = coef*0.65f;
  color[2] = coef*0.88f;
}



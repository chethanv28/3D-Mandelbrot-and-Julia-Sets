#include "disp.h" /* include your own disp.h file (e.g. hw1)*/

/* Camera defaults */
#define	DEFAULT_FOV		45.0
#define	DEFAULT_IM_Z	(0.0)  /* world coords for image plane origin */
#define	DEFAULT_IM_Y	(0.0)    /* default look-at point = 0,0,0 */
#define	DEFAULT_IM_X	(0.0)

#define	DEFAULT_AMBIENT	{0.1, 0.1, 0.1}
#define	DEFAULT_DIFFUSE	{0.7, 0.6, 0.5}
#define	DEFAULT_SPECULAR	{0.2, 0.3, 0.4}
#define	DEFAULT_SPEC		32

//fractal constants
#define BOUNDING_RADIUS_2		3.0
#define ESCAPE_THRESHOLD		1000 //481 in decimal
#define DEL						.0001 

#define	MATLEVELS	100		/* how many matrix pushes allowed */
#define	MAX_LIGHTS	10		/* how many lights allowed */

//#define size 512*512 //(render->display->xres*render->display->yres)

#ifndef GzTexture
#define GzTexture	GzPointer
#endif

#ifndef GZRENDER
#define GZRENDER
typedef struct {			/* define a renderer */
  GzDisplay		*display;
  GzCamera		camera;
  short		    matlevel;	        /* top of stack - current xform */
  GzMatrix		Ximage[MATLEVELS];	/* stack of xforms (Xsm) */
  GzMatrix		Xnorm[MATLEVELS];	/* xforms for norms (Xim) */
  GzMatrix		Xsp;		        /* NDC to screen (pers-to-screen) */
  GzMatrix		Xwi;
  GzColor		flatcolor;          /* color state for flat shaded triangles */
  int			interp_mode;
  int			numlights;
  GzLight		lights[MAX_LIGHTS];
  GzLight		ambientlight;
  GzColor		Ka, Kd, Ks;
  float		    spec;		/* specular power */
  GzTexture		tex_fun;    /* tex_fun(float u, float v, GzColor color) */
}  GzRender;
#endif

// Function declaration
// HW2
int GzNewRender(GzRender **render, GzDisplay *display);
int GzFreeRender(GzRender *render);
int GzBeginRender(GzRender	*render);
int GzPutAttribute(GzRender	*render, int numAttributes, GzToken	*nameList, 
	GzPointer *valueList);
int GzPutRays(GzRender *render, float4 mu, int iterations, int quatType, int colorType, int setType);
int normalizeVector(GzCoord g);
short	ctoi(float color);
int GzPushNormalMatrix(GzRender *render, GzMatrix matrix);
int normalizeMatrix(GzMatrix g);
int GzShader(GzRender* render, GzCoord normal, GzColor color);

// HW3
int GzPutCamera(GzRender *render, GzCamera *camera);
int GzPushMatrix(GzRender *render, GzMatrix	matrix);
int GzPopMatrix(GzRender *render);

// Object Translation
int GzRotXMat(float degree, GzMatrix mat);
int GzRotYMat(float degree, GzMatrix mat);
int GzRotZMat(float degree, GzMatrix mat);
int GzTrxMat(GzCoord translate, GzMatrix mat);
int GzScaleMat(GzCoord scale, GzMatrix mat);

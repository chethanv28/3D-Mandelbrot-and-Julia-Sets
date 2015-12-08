/* 
  disp.h -- cs580 HW1 include file for Display
*/

#include	"gz.h"

/* define general display pixel-type */
#ifndef GZ_PIXEL
typedef	struct {
  GzIntensity    red;	
  GzIntensity    green;
  GzIntensity    blue;
  GzIntensity    alpha;
  GzDepth	 z;
} GzPixel;
#define GZ_PIXEL
#endif;


#ifndef GzFloat3
#define GzFloat3
struct float3
{
	float3() {};
	float3(float s) : x(s), y(s), z(s) {}
	float3(float x, float y, float z) : x(x), y(y), z(z) {}
	float x, y, z;

	inline float3 operator*(float s) const { return float3(x*s, y*s, z*s); }
	inline float3 operator+(const float3& a) const { return float3(x + a.x, y + a.y, z + a.z); }
	inline float3 operator-(const float3& a) const { return float3(x - a.x, y - a.y, z - a.z); }

};
#endif

#ifndef GzFloat4
#define GzFloat4
struct float4
{
	float4() {};
	float4(float s) : x(s), y(s), z(s), w(s) {}
	float4(float3 f, float s) : x(f.x), y(f.y), z(f.z), w(s) {}
	float4(float x, float y, float z, float w) : x(x), y(y), z(z), w(w) {}
	float x, y, z, w;

	inline float4 operator*(float s) const { return float4(x*s, y*s, z*s, w*s); }
	inline float4 operator+(const float4& a) const { return float4(x + a.x, y + a.y, z + a.z, w + a.w); }
	inline float4 operator-(const float4& a) const { return float4(x - a.x, y - a.y, z - a.z, w - a.w); }

};
#endif 


/* define a display type */
#ifndef GZ_DISPLAY
typedef struct {
  unsigned short	xres;
  unsigned short	yres;
  GzPixel		*fbuf;		/* frame buffer array */
} GzDisplay;
#define GZ_DISPLAY
#endif;

/* define bounds on display size in case of error */
#define	MAXXRES	1024	
#define	MAXYRES	1024


// simplify display fbuf indexing 
// (optional - use this or define your own indexing method) 
#define	ARRAY(x,y)	(x+(y*display->xres))	


// Function declarations
int GzNewFrameBuffer(char** framebuffer, int width, int height);
int GzNewDisplay(GzDisplay **display, int xRes, int yRes);
int GzFreeDisplay(GzDisplay	*display);
int GzGetDisplayParams(GzDisplay *display, int *xRes, int *yRes);
int GzInitDisplay(GzDisplay	*display);
int GzPutDisplay(GzDisplay	*display, int i, int j, 
		GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z);
int GzGetDisplay(GzDisplay *display, int i, int j, 
		GzIntensity *r, GzIntensity *g, GzIntensity *b, GzIntensity *a, GzDepth	*z);
int GzFlushDisplay2File(FILE* outfile, GzDisplay *display, int number, float4 mu, int quatType, int colorType, int setType);
int GzFlushDisplay2FrameBuffer(char* framebuffer, GzDisplay* display);

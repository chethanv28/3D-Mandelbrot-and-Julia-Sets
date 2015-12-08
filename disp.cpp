/*   CS580 HW1 display functions to be completed   */

#include   "stdafx.h"  
#include	"Gz.h"
#include	"disp.h"
#include "stdio.h"
#include <iostream>
#include <string>


int GzNewFrameBuffer(char** framebuffer, int width, int height)
{
/* HW1.1 create a framebuffer for MS Windows display:
 -- allocate memory for framebuffer : 3*2 bytes(b, g, r) x width x height Now 4 bytes for each b, g, r
 -- pass back pointer 
 */
	//*framebuffer = new char[12 * width*height];//new char[3 * width*height];
	*framebuffer = new char[3 * width*height];;

	return GZ_SUCCESS;
}

int GzNewDisplay(GzDisplay	**display, int xRes, int yRes)
{
/* HW1.2 create a display:
  -- allocate memory for indicated resolution
  -- pass back pointer to GzDisplay object in display
*/
	*display = new GzDisplay[1];
	display[0]->xres = xRes;
	display[0]->yres = yRes;
	GzPixel *fbuff = new GzPixel[xRes*yRes];
	display[0]->fbuf = fbuff;
	return GZ_SUCCESS;
}


int GzFreeDisplay(GzDisplay	*display)
{
/* HW1.3 clean up, free memory */
	free(display);
	return GZ_SUCCESS;
}


int GzGetDisplayParams(GzDisplay *display, int *xRes, int *yRes)
{
/* HW1.4 pass back values for a display */
	//int *xres = new (int)(display->xres);
	xRes = new int[1];
	xRes[0] = (int)(display->xres);
	yRes = new int[1]; 
	yRes[0] = (int)(display->yres);
	return GZ_SUCCESS;
}


int GzInitDisplay(GzDisplay	*display)
{
/* HW1.5 set everything to some default values - start a new frame */
	int xres = display->xres;
	int yres = display->yres;
	int size = xres*yres;
	for (int i = 0; i < xres; i++){
		for (int j = 0; j < yres; j++){
			display->fbuf[ARRAY(i, j)].red = 0;
			display->fbuf[ARRAY(i, j)].blue = 0;
			display->fbuf[ARRAY(i, j)].green = 0;
			display->fbuf[ARRAY(i, j)].alpha = 1;
			display->fbuf[ARRAY(i, j)].z = MAXINT;
		}
	}
	return GZ_SUCCESS;
}


int GzPutDisplay(GzDisplay *display, int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z)
{
/* HW1.6 write pixel values into the display */
	int *xres = 0;
	int *yres =0;
	GzGetDisplayParams(display, xres, yres);

	if (r < 0){
		r = 0;
	}
	if (g < 0){
		g = 0;
	}
	if (b < 0){
		b = 0;
	}
	if (r>4095){
		r = 4095;
	}
	if (g>4095){
		g = 4095;
	}
	
	if (b>4095){
		b = 4095;
	}
	
	if (i >= 0 && i < display->xres && j >= 0 && j < display->yres){
		//if (display->fbuf[ARRAY(i, j)] != NULL){
		int arrayval = ARRAY(i, j);
			display->fbuf[ARRAY(i, j)].red = r;
			display->fbuf[ARRAY(i, j)].blue = b;
			display->fbuf[ARRAY(i, j)].green = g;
			display->fbuf[ARRAY(i, j)].alpha = a;
			display->fbuf[ARRAY(i, j)].z = z;
		//}
	}

	return GZ_SUCCESS;
}


int GzGetDisplay(GzDisplay *display, int i, int j, GzIntensity *r, GzIntensity *g, GzIntensity *b, GzIntensity *a, GzDepth *z)
{
	/* HW1.7 pass back a pixel value to the display */
	if (i >= 0 && i < display->xres && j >= 0 && j < display->yres){
		*r = display->fbuf[ARRAY(i, j)].red;
		*g = display->fbuf[ARRAY(i, j)].green;
		*b = display->fbuf[ARRAY(i, j)].blue;
		*a = display->fbuf[ARRAY(i, j)].alpha;
		*z = display->fbuf[ARRAY(i, j)].z;//display->fbuf[ARRAY(j, i)].z;
	}
	return GZ_SUCCESS;
}


int GzFlushDisplay2File(FILE* outfile, GzDisplay *display, int number, float4 mu, int quatType, int colorType, int setType)
{
/* HW1.8 write pixels to ppm file -- "P6 %d %d 255\r" */
	int xdim = (int)display->xres;
	int ydim = (int)display->yres;

	const char *p = "output";
	std::string s = std::to_string(number);
	const char *three = ".ppm";
	std::string total = std::string(p) + s + std::string(three);

	//total.c_str();
	outfile = fopen(total.c_str(), "wb");

	(void) fprintf(outfile, "P6 %d %d 255\r", xdim, ydim);

	for (int j = 0; j < ydim; j++){
		for (int i = 0; i < xdim; i++){
			static unsigned char color[3];
			color[0] = (byte)(display->fbuf[ARRAY(i, j)].red >> 4);
			color[1] = (byte)(display->fbuf[ARRAY(i, j)].green >> 4);
			color[2] = (byte)(display->fbuf[ARRAY(i, j)].blue >> 4);
			(void)fwrite(color, 1, 3, outfile);
			if (i == 130 && j == 230){
				int seven = 711;
			}
			if (i == 131 && j == 230){
				int seven = 711;
			}
			if (i == 130 && j == 231){
				int seven = 711;
			}
			if (i == 131 && j == 231){
				int seven = 711;
			}
		}
	}

	(void) fclose(outfile);

	return GZ_SUCCESS;
}

int GzFlushDisplay2FrameBuffer(char* framebuffer, GzDisplay *display)
{

/* HW1.9 write pixels to framebuffer: 
	- put the pixels into the frame buffer
	- CAUTION: when storing the pixels into the frame buffer, the order is blue, green, and red 
	- NOT red, green, and blue !!!
*/
	int xdim = (int)display->xres;
	int ydim = (int)display->yres;
	for (int j = 0; j < ydim; j++){
		for (int i = 0; i < xdim; i++){
			static unsigned char color[3];
			color[0] = (byte)(display->fbuf[ARRAY(i, j)].red >> 4);
			int colorr = display->fbuf[ARRAY(i, j)].red;
			int colorr1 = display->fbuf[ARRAY(i, j)].red >> 4;
			color[1] = (byte)(display->fbuf[ARRAY(i, j)].green >> 4);
			color[2] = (byte)(display->fbuf[ARRAY(i, j)].blue >> 4);
			framebuffer[3*ARRAY(i,j)] = color[2];
			framebuffer[3 * ARRAY(i,j) + 1] = color[1];
			framebuffer[3 * ARRAY(i,j) + 2] = color[0];
			if (i == 130 && j == 230){
				int seven = 711;
			}
			if (i == 131 && j == 230){
				int seven = 711;
			}
			if (i == 130 && j == 231){
				int seven = 711;
			}
			if (i == 131 && j == 231){
				int seven = 711;
			}
		}
	}
	
	return GZ_SUCCESS;
}
/** \addtogroup auxilliary
Some convenient helpers that draw various standard objects
*  @{
*/

/*!
\file Primitives.h
*/

#ifndef _PRIMITIVES_H_
#define _PRIMITIVES_H_
#include "libheaders.h"
#include "glerror.h"

//! Segment count for spheres
#define PRIM_SPHERE_SEGMENTS 16
//! Ring count for spheres
#define PRIM_SPHERE_RINGS 16

/*!
\brief Implements some helpers to draw basic primitives
*/
class Primitives
{
public:
	//! Draws a (-1, -1, -1) (1, 1, 1) cube
	static void drawNDCCube();
	//! Draws a (-1, -1) (1, 1) quad
	static void drawNDCQuad();
	//! Draws a sphere centered at (0, 0, 0) with a radius of 1
	static void drawNDCSphere();

private:
	static GLuint cubevbo;
	static GLuint cubevao;
	static GLuint quadvbo;
	static GLuint quadvao;
	static GLuint spherevbo;
	static GLuint spherevao;
};
#endif

/** @}*/
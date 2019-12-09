/** \addtogroup scene_components
Model class
*  @{
*/

/*!
\file Model.h
*/

#ifndef _MODEL_H_
#define _MODEL_H_
#include <Transform.h>
#include <vector>
#include <Mesh.h>
#include <Shader.h>

/*!
\brief Class to draw and place a set of associated meshes in the scene
*/
class Model
{
public:
	/**
	\brief Draws all meshes without setting material uniforms, using the transformation of this model.
	\param shader Shader to be drawn with
	*/
	void draw(ShaderProgram* shader, size_t rendererid = 0);
	/**
	\brief Draws all meshes with setting material uniforms, using the transformation of this model.
	\param shader Shader to be drawn with
	*/
	void drawWithMaterials(ShaderProgram* shader, size_t rendererid = 0);

	//! World transformation of this model
	Transform transform;
	//! Collection of meshes belonging to this model.
	std::vector<Mesh*> meshes;
};

#endif
/** @}*/
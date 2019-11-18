/** \addtogroup utility
Scene loading
*  @{
*/

/*!
\file SceneLoader.h
*/

#ifndef _SCENE_LOADER_H_
#define _SCENE_LOADER_H_
#include <Scene.h>
#include <string>
#include <memory>

class Scanner;

/*!
\brief Loads a scene and its referenced resources from a scene file.
*/
class SceneLoader
{
public:
	/*!
	\brief Loads a scene and its referenced resources from a scene file.
	\param scenepath	Path to the scene file
	
	All resources paths must either be absolute paths or originate from the current working directory.

	Scene file with syntax explained in the comments:
	\code
	# I'm a comment

	# Camera[1]
	# Syntax:
	# CAMERA campos(x y z) lookat(x y z) up(x y z) fovy near far
	CAMERA 0.0 0.0 2.0 0.0 0.0 0.0 0.0 1.0 0.0 90.0 0.1 100.0

	# Opaque Model[0..n]
	# Syntax:
	# Note: rotations are given as 4-component quaternions
	# MODEL_OPAQUE "path-to-some-modelfile"  position(x y z)  rotation(w x y z)  scale(x y z)
	MODEL_OPAQUE "assets/models/sponza/sponza.obj" 0.0 0.0 0.0  1.0 0.0 0.0 0.0  0.05 0.05 0.05

	# Transparent Models[0..n]
	# Syntax:
	# Note: rotations are given as 4-component quaternions
	# MODEL_TRANSPARENT "path-to-some-modelfile"  position(x y z)  rotation(w x y z)  scale(x y z)
	MODEL_TRANSPARENT "assets/models/dragon/dragon.obj" 0.0 1.5 0.0  1.0 0.0 0.0 0.0  4.0 4.0 4.0

	# Particle System[0..n]
	# Syntax:
	# maxParticles spawnrate lifetime position(x y z) direction(x y z) particleSize emitAngle velocityRange force mass {optional: color to be multiplied with the diff texture}  {optional: alphaExponent}
	PARTICLE_EMITTER 1000  200.0  3.0  -2.0 0.2 0.5  0.0 1.0 -0.2  1.5  0.2  3.0 8.0  0.0 -3.0 0.0  1.0  "assets/particles/smoke_white3.png" 1.0 0.0 0.0  1.5

	# Directional Light[0..n]
	# Syntax:
	# DIRECTIONAL_LIGHT intensity color(r g b) direction(x y z)
	DIRECTIONAL_LIGHT 10.0 1.0 1.0 1.0 -1.0 -1.0 -1.0

	# Point Light[0..n]
	# Syntax:
	# POINT_LIGHT intensity color(r g b) position(x y z)
	POINT_LIGHT 7.0 0.0 0.0 1.0 4.0 0.0 0.0

	# Ambient Light[0..n]
	# Syntax:
	# AMBIENT_LIGHT intensity color(r g b)
	AMBIENT_LIGHT 0.2 1.0 1.0 1.0

	# Tone Mapper settings (White threshold)
	# Syntax:
	# MAX_LUM maximumLuminance
	MAX_LUM 5000.0

	# Enable/Disable Normal Mapping globally
	# Syntax:
	# ENABLE_NORMAL_MAPPING (0 or 1)
	ENABLE_NORMAL_MAPPING 0

	# Initial maximum number of depth peeling passes
	# 0 means: peel as long as there are fragments left
	# Syntax:
	# DEPTH_PEEL_MAX_PASSES numberOfInitialDepthPeelPasses
	DEPTH_PEEL_MAX_PASSES 10

	# Initial bucket setup
	# Syntax:
	# NUM_BUCKETS maxBuckets numInitialBuckets firstBucketStartDepth lastBucketEndDepth
	NUM_BUCKETS 64 16 0.0 20.0

	# Offset of depth weighting function from first bucket
	# The the weighting function is shifted w.r.t. the z axis pointing towards the viewer.
	# Prevents too-large weights for the first bucket. 
	DEPTH_WEIGHT_OFFSET 1.0

	# Maximum Fragment Buffer size
	# Actual number is calculated as: screenres.x * screenres.y * FRAGLIST_SIZE_MUL
	# Syntax:
	# FRAGLIST_SIZE_MUL maximumAverageFragmentsPerPixel
	FRAGLIST_SIZE_MUL 20
	\endcode
	*/
	static std::unique_ptr<Scene> loadScene(const std::string& scenepath);
private:
	static void parseFile(Scanner* scan, Scene* scn);
	static void parseCamera(Scanner* scan, Scene* scn);
	static void parseModelOpaque(Scanner* scan, Scene* scn);
	static void parseModelTransparent(Scanner* scan, Scene* scn);
	static void parseDirLight(Scanner* scan, Scene* scn);
	static void parsePointLight(Scanner* scan, Scene* scn);
	static void parseAmbientLight(Scanner* scan, Scene* scn);
	static void parseExposure(Scanner* scan, Scene* scn);
	static void parseSkyboxRes(Scanner* scan, Scene* scn);
	static void parseClearColor(Scanner* scan, Scene* scn);
	static void parseMaterial(Scanner* scan, Scene* scn);
	static void parseEnvMap(Scanner* scan, Scene* scn);
};



#endif // !_SCENE_LOADER_H_

/** @}*/
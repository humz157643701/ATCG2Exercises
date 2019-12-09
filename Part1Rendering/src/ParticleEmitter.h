/** \addtogroup scene_components
Particle system
*  @{
*/

/*!
\file ParticleEmitter.h
*/

#ifndef _PARTICLE_EMITTER_H_
#define _PARTICLE_EMITTER_H_
#include <Transform.h>
#include <vector>
#include <Shader.h>
#include <libheaders.h>
#include <Material.h>
#include <random>

//! Size of the random texture
#define RANDOM_TEXTURE_SIZE 1024 * 1024

//! Work group size to be used for the compute shaders
#define WORK_GROUP_SIZE 1024

/**
\brief Implements a GPU-based particle system
*/
class ParticleEmitter
{
public:
	/**
	\brief Constructs a particle system.

	\param _numparticles		Maximum number of particles available
	\param _spawnrate			Number of particles to spawn per second
	\param _lifetime			Initial lifetime of all particles
	\param _position			Position of the emitter in world space
	\param _direction			Direction of the emitter in world space
	\param _particleSize		Width and height of the particle billboards
	\param _emitAngle			The maximum angle to spawn particles into in radians, centered around _direction
	\param _startVelocityRange	Range of velocities which every particle starts with. A random value is chosen in the range [_startVelocityRange.x, _startVelocityRange.y].
	\param _force				A constant force that acts on the particles
	\param _mass				Mass of each particle
	\param _material			Material of the particles. Only the diffuse texture/color is used during particle rendering.
	\param _alphaexp			Modifies how the alpha of a particle changes over time. 1.0 means linear falloff.
	*/
	ParticleEmitter(
		size_t _numparticles,
		float _spawnrate,
		float _lifetime,
		const glm::vec3& _position,
		const glm::vec3& _direction,
		float _particleSize,
		float _emitAngle,
		const glm::vec2& _startVelocityRange,
		const glm::vec3& _force,
		float _mass,
		Material* _material,
		float _alphaexp
	);

	~ParticleEmitter();

	ParticleEmitter(const ParticleEmitter& other) = delete;
	ParticleEmitter& operator=(const ParticleEmitter& other) = delete;
	ParticleEmitter(ParticleEmitter&& other);
	ParticleEmitter& operator=(ParticleEmitter&& other);

	//! Performs a simulation step for the particle system
	void update(float dt);
	//! Draws the particles without setting material unidorms
	void draw(ShaderProgram* shader);
	//! Draws the particles with setting material uniforms
	void drawWithMaterials(ShaderProgram* shader);

	ShaderProgram emitshader;
	ShaderProgram simshader;

	//! Transformation of the emitter
	Transform transform;
	glm::mat3 orientation;
	bool orientationDirty;

	GLuint vbo;
	GLuint ibo;
	GLuint vao;

	//position lifetime velocity rotation
	GLuint b_pdata;
	GLuint b_deadparticles;
	GLuint b_aliveparticles1;
	GLuint b_aliveparticles2;

	GLuint acb_counters; //dead particles, alive particles, particles to render

	//for getting some pseudo random values
	GLuint b_random;
	GLuint t_random;

	size_t numParticles;
	Material* material;
	float particleSize;
	float mass;
	float spawnrate;
	float lifetime;
	float alphaexponent;

	glm::vec3 force;
	glm::vec3 emitDirection;
	float emitAngle;
	glm::vec2 startVelocityRange;

	GLuint numAliveParticles;
	GLuint numDeadParticles;

	//time accumulator for correct spawn behaviour
	float taccum;

	//for random offsets for the noise texture
	std::default_random_engine randomEngine;
	std::uniform_int_distribution<GLuint> rintDist;

	GLuint getInvocationCount(GLuint n);
};
#endif

/** @}*/
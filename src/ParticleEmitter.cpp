#include "ParticleEmitter.h"
#include <random>
#include <vector>

ParticleEmitter::ParticleEmitter(
	size_t _numparticles,
	float _spawnrate,
	float _lifetime,
	const glm::vec3 & _position,
	const glm::vec3 & _direction,
	float _particleSize,
	float _emitAngle,
	const glm::vec2& _startVelocityRange,
	const glm::vec3 & _force,
	float _mass,
	Material * _material,
	float _alphaexp) :
	numParticles(_numparticles),
	spawnrate(_spawnrate),
	lifetime(_lifetime),
	particleSize(_particleSize),
	emitAngle(_emitAngle),
	startVelocityRange(_startVelocityRange),
	force(_force),
	mass(_mass),
	material(_material),
	emitshader(ShaderProgram::createComputeShaderProgram("assets/shaders/particlesystem/emit.comp", WORK_GROUP_SIZE, 1, 1)),
	simshader(ShaderProgram::createComputeShaderProgram("assets/shaders/particlesystem/sim.comp", WORK_GROUP_SIZE, 1, 1)),
	taccum(0.0f),
	randomEngine(),
	rintDist(0, RANDOM_TEXTURE_SIZE),
	alphaexponent(_alphaexp)
{
	// setup transform
	glm::vec3 az = glm::normalize(-_direction);
	glm::vec3 tmpup{ 0.0f, 1.0f, 0.0f };
	//if the direction is dangerously close to the canonical up axis, choose another axis for up
	if (glm::abs(glm::dot(tmpup, az)) >= 0.99f)	
		tmpup = glm::vec3{ 1.0f, 0.0f, 0.0f };
	glm::vec3 ax = glm::normalize(glm::cross(tmpup, az));
	glm::vec3 ay = glm::normalize(glm::cross(az, ax));
	glm::mat3 rotmat{ ax, ay, az };
	transform = Transform(_position, glm::quat_cast(rotmat), glm::vec3(1.0f));
	orientation = glm::mat3(
		transform.getWorldXAxis(),
		transform.getWorldYAxis(),
		transform.getWorldZAxis()
	);
	orientationDirty = false;

	//geometry data for particle
	GLfloat vertexdata[]{
		-0.5f, -0.5f, 0.0f,
		0.0f, 0.0f,
		0.0f, 0.0f, 1.0f,
		1.0f, 0.0f, 0.0f,
		0.0f, 1.0f, 0.0f,

		0.5f, -0.5f, 0.0f,
		1.0f, 0.0f,
		0.0f, 0.0f, 1.0f,
		1.0f, 0.0f, 0.0f,
		0.0f, 1.0f, 0.0f,

		0.5f, 0.5f, 0.0f,
		1.0f, 1.0f,
		0.0f, 0.0f, 1.0f,
		1.0f, 0.0f, 0.0f,
		0.0f, 1.0f, 0.0f,

		-0.5f, 0.5f, 0.0f,
		0.0f, 1.0f,
		0.0f, 0.0f, 1.0f,
		1.0f, 0.0f, 0.0f,
		0.0f, 1.0f, 0.0f
	};
	GLuint indexdata[]{
		0, 1, 2,
		0, 2, 3
	};

	//create buffers

	//vertex data
	glGenBuffers(1, &vbo);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertexdata), reinterpret_cast<const void*>(&vertexdata[0]), GL_STATIC_DRAW);
	
	//index data
	glGenBuffers(1, &ibo);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indexdata), reinterpret_cast<const void*>(&indexdata[0]), GL_STATIC_DRAW);

	//position and lifetime data
	glGenBuffers(1, &b_pdata);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, b_pdata);
	glBufferStorage(GL_SHADER_STORAGE_BUFFER, numParticles * 8 * sizeof(GLfloat), reinterpret_cast<const void*>(0), 0);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

	//dead particles list
	glGenBuffers(1, &b_deadparticles);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, b_deadparticles);
	glBufferStorage(GL_SHADER_STORAGE_BUFFER, numParticles * sizeof(GLuint), reinterpret_cast<const void*>(0), GL_MAP_WRITE_BIT);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

	//alive particles list 1
	glGenBuffers(1, &b_aliveparticles1);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, b_aliveparticles1);
	glBufferStorage(GL_SHADER_STORAGE_BUFFER, numParticles * sizeof(GLuint), reinterpret_cast<const void*>(0), 0);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

	//alive particles list 2
	glGenBuffers(1, &b_aliveparticles2);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, b_aliveparticles2);
	glBufferStorage(GL_SHADER_STORAGE_BUFFER, numParticles * sizeof(GLuint), reinterpret_cast<const void*>(0), 0);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

	//atomic counters buffer (#dead particles, #aliveparticles, #particlestorender)
	glGenBuffers(1, &acb_counters);
	glBindBuffer(GL_ATOMIC_COUNTER_BUFFER, acb_counters);
	glBufferStorage(GL_ATOMIC_COUNTER_BUFFER, 3 * sizeof(GLuint), reinterpret_cast<const void*>(0), GL_MAP_WRITE_BIT | GL_MAP_READ_BIT);
	glBindBuffer(GL_ATOMIC_COUNTER_BUFFER, 0);

	//setup vao
	//in vertex shader:
	// 0: vec3 position		|BUF1
	// 1: vec2 uv			|
	// 2: vec3 normal		|	once for all particles
	// 3: vec3 tangent		|
	// 4: vec3 bitangent	|
	const size_t vdstride = (3 + 2 + 3 + 3 + 3) * 4;

	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);

	//attributes 0 .. 4	
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, vdstride, reinterpret_cast<const void*>(0));
	glVertexAttribDivisor(0, 0);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, vdstride, reinterpret_cast<const void*>(3 * 4));
	glVertexAttribDivisor(1, 0);
	glEnableVertexAttribArray(2);
	glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, vdstride, reinterpret_cast<const void*>((3 + 2) * 4));
	glVertexAttribDivisor(2, 0);
	glEnableVertexAttribArray(3);
	glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, vdstride, reinterpret_cast<const void*>((3 + 2 + 3) * 4));
	glVertexAttribDivisor(3, 0);
	glEnableVertexAttribArray(4);
	glVertexAttribPointer(4, 3, GL_FLOAT, GL_FALSE, vdstride, reinterpret_cast<const void*>((3 + 2 + 3 + 3) * 4));
	glVertexAttribDivisor(4, 0);

	//indices
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);

	glBindVertexArray(0);

	//initialize simulation

	// initialize counter buffers
	glBindBuffer(GL_ATOMIC_COUNTER_BUFFER, acb_counters);
	GLuint* acbuf = static_cast<GLuint*>(glMapBufferRange(GL_ATOMIC_COUNTER_BUFFER, 0, 3 * sizeof(GLuint), GL_MAP_WRITE_BIT));
	acbuf[0] = static_cast<GLuint>(numParticles);		//dead particle count
	acbuf[1] = static_cast<GLuint>(0);					//alive particle count
	acbuf[2] = static_cast<GLuint>(0);					//currently unused
	glUnmapBuffer(GL_ATOMIC_COUNTER_BUFFER);

	//initialize dead particle list with indices from 0 to numparticles - 1
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, b_deadparticles);
	GLuint* dpcbuf = static_cast<GLuint*>(glMapBufferRange(GL_SHADER_STORAGE_BUFFER, 0, numParticles * sizeof(GLuint), GL_MAP_WRITE_BIT));
	for (GLuint i = 0; i < numParticles; ++i)
		dpcbuf[i] = i;
	glUnmapBuffer(GL_SHADER_STORAGE_BUFFER);

	glMemoryBarrier(GL_ALL_BARRIER_BITS);

	//create texture with random values
	std::vector<GLfloat> rdata;
	rdata.reserve(RANDOM_TEXTURE_SIZE * 4);
	std::default_random_engine eng;
	std::uniform_real_distribution<GLfloat> dist(0.0f, 1.0f);
	for (size_t i = 0; i < RANDOM_TEXTURE_SIZE * 4; ++i)
	{
		rdata.push_back(dist(eng));
	}

	glGenBuffers(1, &b_random);
	glBindBuffer(GL_TEXTURE_BUFFER, b_random);
	glBufferData(GL_TEXTURE_BUFFER, RANDOM_TEXTURE_SIZE * 4 * sizeof(GLfloat), rdata.data(), GL_STATIC_DRAW);
	glBindBuffer(GL_TEXTURE_BUFFER, 0);

	glGenTextures(1, &t_random);
	glBindTexture(GL_TEXTURE_BUFFER, t_random);
	glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA32F, b_random);
	glBindTexture(GL_TEXTURE_BUFFER, 0);

	numAliveParticles = 0;
}

ParticleEmitter::~ParticleEmitter()
{
	if (ibo)
		glDeleteBuffers(1, &ibo);
	if (vbo)
		glDeleteBuffers(1, &vbo);
	if (vao)
		glDeleteVertexArrays(1, &vao);

	if (b_pdata)
		glDeleteBuffers(1, &b_pdata);
	//if (b_velocities)
	//	glDeleteBuffers(1, &b_velocities);
	if (b_deadparticles)
		glDeleteBuffers(1, &b_deadparticles);
	if (b_aliveparticles1)
		glDeleteBuffers(1, &b_aliveparticles1);
	if (b_aliveparticles2)
		glDeleteBuffers(1, &b_aliveparticles2);
	if (acb_counters)
		glDeleteBuffers(1, &acb_counters);

	if (t_random)
		glDeleteTextures(1, &t_random);
	if (b_random)
		glDeleteBuffers(1, &b_random);
}

ParticleEmitter::ParticleEmitter(ParticleEmitter && other) :
	emitshader(std::move(other.emitshader)),
	simshader(std::move(other.simshader)),
	transform(std::move(other.transform)),
	vbo(other.vbo),
	ibo(other.ibo),
	vao(other.ibo),
	b_pdata(other.b_pdata),
	//b_velocities(other.b_velocities),
	b_deadparticles(other.b_deadparticles),
	b_aliveparticles1(other.b_aliveparticles1),
	b_aliveparticles2(other.b_aliveparticles2),
	acb_counters(other.acb_counters),
	numParticles(other.numParticles),
	material(other.material),
	particleSize(other.particleSize),
	mass(other.mass),
	spawnrate(other.spawnrate),
	lifetime(other.lifetime),
	force(other.force),
	emitDirection(other.emitDirection),
	emitAngle(other.emitAngle),
	startVelocityRange(other.startVelocityRange),
	t_random(other.t_random),
	b_random(other.b_random),
	orientation(other.orientation),
	orientationDirty(other.orientationDirty),
	numAliveParticles(0),
	taccum(0.0f),
	randomEngine(),
	rintDist(0, RANDOM_TEXTURE_SIZE),
	alphaexponent(other.alphaexponent)
{		
		other.vbo					= 0;
		other.ibo					= 0;
		other.vao					= 0;
		other.b_pdata				= 0;
		//other.b_velocities			= 0;
		other.b_deadparticles		= 0;
		other.b_aliveparticles1		= 0;
		other.b_aliveparticles2		= 0;
		other.acb_counters			= 0;
		other.t_random = 0;
		other.b_random = 0;
}

ParticleEmitter & ParticleEmitter::operator=(ParticleEmitter && other)
{
	if (this == &other)
		return *this;

	if (ibo)
		glDeleteBuffers(1, &ibo);
	if (vbo)
		glDeleteBuffers(1, &vbo);
	if (vao)
		glDeleteVertexArrays(1, &vao);

	if (b_pdata)
		glDeleteBuffers(1, &b_pdata);
	//if (b_velocities)
	//	glDeleteBuffers(1, &b_velocities);
	if (b_deadparticles)
		glDeleteBuffers(1, &b_deadparticles);
	if (b_aliveparticles1)
		glDeleteBuffers(1, &b_aliveparticles1);
	if (b_aliveparticles2)
		glDeleteBuffers(1, &b_aliveparticles2);
	if (acb_counters)
		glDeleteBuffers(1, &acb_counters);

	if (t_random)
		glDeleteTextures(1, &t_random);

	if (b_random)
		glDeleteBuffers(1, &b_random);

	emitshader = std::move(other.emitshader);
	simshader = std::move(other.simshader);
	transform = std::move(other.transform);
	orientation = other.orientation;
	orientationDirty = other.orientationDirty;

	numParticles = other.numParticles;
	material = other.material;
	particleSize = other.particleSize;
	mass = other.mass;
	spawnrate = other.spawnrate;
	lifetime = other.lifetime;
	force = other.force;
	emitDirection = other.emitDirection;
	emitAngle = other.emitAngle;
	startVelocityRange = other.startVelocityRange;
	alphaexponent = other.alphaexponent;

	vbo = other.vbo;
	ibo = other.ibo;
	vao = other.vao;
	b_pdata = other.b_pdata;
	//b_velocities = other.b_velocities;
	b_deadparticles = other.b_deadparticles;
	b_aliveparticles1 = other.b_aliveparticles1;
	b_aliveparticles2 = other.b_aliveparticles2;
	acb_counters = other.acb_counters;

	t_random = other.t_random;
	b_random = other.b_random;

	other.vbo = 0;
	other.ibo = 0;
	other.vao = 0;
	other.b_pdata = 0;
	//other.b_velocities = 0;
	other.b_deadparticles = 0;
	other.b_aliveparticles1 = 0;
	other.b_aliveparticles2 = 0;
	other.acb_counters = 0;
	other.t_random = 0;
	other.b_random = 0;

	numAliveParticles = 0;
	taccum = 0.0f;

	return *this;
}

void ParticleEmitter::update(float dt)
{
	if (orientationDirty)
	{
		orientation = glm::mat3(
			transform.getWorldXAxis(),
			transform.getWorldYAxis(),
			transform.getWorldZAxis()
		);
		orientationDirty = false;
	}

	//read dead count
	glMemoryBarrier(GL_BUFFER_UPDATE_BARRIER_BIT);
	glBindBufferBase(GL_ATOMIC_COUNTER_BUFFER, 0, acb_counters);
	GLuint* buf = static_cast<GLuint*>(glMapBufferRange(GL_ATOMIC_COUNTER_BUFFER, 0, 3 * sizeof(GLuint), GL_MAP_READ_BIT));
	numDeadParticles = buf[0];
	glUnmapBuffer(GL_ATOMIC_COUNTER_BUFFER);

	//calculate number of particles to emit
	GLuint pts = 0;
	float psi = 1.0f / spawnrate;
	taccum += dt;
	while (taccum > psi && pts < numDeadParticles)
	{
		++pts;
		taccum -= psi;
	}

	if (pts > 0)
	{
		//emit pass
		glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT | GL_ATOMIC_COUNTER_BARRIER_BIT);
		emitshader.use();
		glBindBufferBase(GL_ATOMIC_COUNTER_BUFFER, 0, acb_counters);
		glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, b_pdata);
		glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, b_deadparticles);
		glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, b_aliveparticles1);

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_BUFFER, t_random);
		emitshader.setUniform("t_random", 0);
		emitshader.setUniform("numRandomTexels", static_cast<GLuint>(RANDOM_TEXTURE_SIZE));
		emitshader.setUniform("randomOffset", rintDist(randomEngine));

		emitshader.setUniform("orientation", orientation, false);
		emitshader.setUniform("position", transform.getWorldPosition());
		emitshader.setUniform("emitAngle", emitAngle);
		emitshader.setUniform("startVelocityRange", startVelocityRange);
		emitshader.setUniform("lifetime", lifetime);
		emitshader.setUniform("numParticles", pts);
		
		glDispatchCompute(getInvocationCount(pts), 1, 1);
	}

	//read alive count and set to 0
	glMemoryBarrier(GL_BUFFER_UPDATE_BARRIER_BIT);
	glBindBufferBase(GL_ATOMIC_COUNTER_BUFFER, 0, acb_counters);
	buf = static_cast<GLuint*>(glMapBufferRange(GL_ATOMIC_COUNTER_BUFFER, 0, 3 * sizeof(GLuint), GL_MAP_READ_BIT | GL_MAP_WRITE_BIT));
	numAliveParticles = buf[1];	
	buf[1] = 0;
	glUnmapBuffer(GL_ATOMIC_COUNTER_BUFFER);

	//sim pass
	if (numAliveParticles)
	{
		glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT | GL_ATOMIC_COUNTER_BARRIER_BIT);
		simshader.use();
		glBindBufferBase(GL_ATOMIC_COUNTER_BUFFER, 0, acb_counters);
		glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, b_pdata);
		glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, b_deadparticles);
		glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, b_aliveparticles1);
		glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, b_aliveparticles2);

		simshader.setUniform("dt", dt);
		simshader.setUniform("force", force);
		simshader.setUniform("mass", mass);
		simshader.setUniform("numParticles", numAliveParticles);
		
		glDispatchCompute(getInvocationCount(numAliveParticles), 1, 1);

		//barrier
		glMemoryBarrier(GL_BUFFER_UPDATE_BARRIER_BIT);
		//read new alive count (for instancing)
		glBindBufferBase(GL_ATOMIC_COUNTER_BUFFER, 0, acb_counters);
		buf = static_cast<GLuint*>(glMapBufferRange(GL_ATOMIC_COUNTER_BUFFER, 0, 3 * sizeof(GLuint), GL_MAP_READ_BIT));
		numAliveParticles = buf[1];
		glUnmapBuffer(GL_ATOMIC_COUNTER_BUFFER);	

		//swap alive lists
		GLuint tmp = b_aliveparticles1;
		b_aliveparticles1 = b_aliveparticles2;
		b_aliveparticles2 = tmp;
	}
}

void ParticleEmitter::draw(ShaderProgram * shader)
{
	if (numAliveParticles)
	{
		glMemoryBarrier(GL_ATOMIC_COUNTER_BARRIER_BIT | GL_SHADER_STORAGE_BARRIER_BIT);
		glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, b_aliveparticles1);
		glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, b_pdata);
		shader->setUniform("model_matrix", transform.getLocalToWorldMatrix(), false);
		shader->setUniform("particleSize", particleSize);
		shader->setUniform("particleLifetime", lifetime);
		shader->setUniform("particleExponent", alphaexponent);
		shader->setUniform("isParticle", 1);
		glBindVertexArray(vao);
		glDrawElementsInstanced(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0, numAliveParticles);
		glBindVertexArray(0);
		shader->setUniform("isParticle", 0);
	}
}

void ParticleEmitter::drawWithMaterials(ShaderProgram * shader)
{
	assert(material);
	shader->saveTU();
	material->bind(shader);
	draw(shader);
	shader->restoreTU();
}

GLuint ParticleEmitter::getInvocationCount(GLuint n)
{
	return n / WORK_GROUP_SIZE + static_cast<GLuint>(n % WORK_GROUP_SIZE != 0);
}

#version 440 core
    
layout(location = 0) in vec3 position;

struct Camera
{
    mat4 viewMatrix;
    mat4 projectionMatrix;
    vec3 worldPos;
    float renderWidth;
    float renderHeight;
    float nearPlane;
    float farPlane;    
};
    
uniform mat4 model_matrix;
uniform Camera camera;

//particle stuff
layout(std430, binding = 2) buffer particleListLayout
{
  uint pidx[];
} sb_particleIndexList;

struct Particle
{
    vec3 pos;
    float lt;
    vec3 vel;
    float rot;
};

layout(binding = 3, std430) buffer pdataLayout
{
    Particle particles[];
} sb_pdata;

uniform bool isParticle = false;
uniform float particleSize = 1.0f;

out float viewSpaceDepth;

void main()
{
    mat4 mmat = model_matrix;

    if(isParticle)
    {
        Particle p = sb_pdata.particles[sb_particleIndexList.pidx[gl_InstanceID]];
        vec3 campos = camera.worldPos;

        vec3 pz = normalize(campos - p.pos);
        vec3 px = normalize(cross(vec3(0.0, 1.0, 0.0), pz));
        vec3 py = cross(pz, px);

        //rotation matrix for random rotation
        mat4 rmat = mat4(
            vec4(cos(p.rot), sin(p.rot), 0.0 , 0.0),
            vec4(-sin(p.rot), cos(p.rot), 0.0 , 0.0),
            vec4(0.0, 0.0, 1.0, 0.0),
            vec4(0.0, 0.0, 0.0, 1.0)
        );

        mat4 smat = mat4(
            vec4(particleSize, 0.0, 0.0, 0.0),
            vec4(0.0, particleSize, 0.0, 0.0),
            vec4(0.0, 0.0, particleSize ,0.0),
            vec4(0.0, 0.0, 0.0 , 1.0)
        );

        mat4 tmat = mat4(
             vec4(px, 0.0),
             vec4(py, 0.0),
             vec4(pz, 0.0),
             vec4(p.pos, 1.0)
        );

        mmat = tmat * rmat * smat;
    }
    vec4 viewPos = camera.viewMatrix * mmat * vec4(position, 1.0f);
    gl_Position = camera.projectionMatrix * viewPos;
    viewSpaceDepth = viewPos.z;
}
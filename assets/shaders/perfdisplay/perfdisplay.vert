#version 440 core
    
layout(location = 0) in vec3 position;
layout(location = 1) in vec2 uv;
layout(location = 2) in vec3 normal;
layout(location = 3) in vec3 tangent;
layout(location = 4) in vec3 bitangent;

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

struct Material
{
    bool constDiff;
    sampler2D diff;
    vec4 cdiff;
    bool constSpec;
    sampler2D spec;
    vec4 cspec;
    bool constShin;
    sampler2D shin;
    float cshin;
    bool normalMap;
    sampler2D norm;    
    bool nmisbm;
};
uniform Material material;
    
uniform mat4 model_matrix;
uniform Camera camera;
    
out struct VertexData
{
    vec3 viewspacePosition;
    vec2 uv;
    vec3 normal;
    mat3 tbn;
} vertexData;

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
   
    mat4 modelview = camera.viewMatrix * mmat;
    vec4 vpos = modelview * vec4(position, 1.0f);
    gl_Position = camera.projectionMatrix * vpos;

    mat3 nmat = inverse(transpose(mat3(modelview)));
    vec3 N = normalize(nmat * normal);

    vertexData.viewspacePosition = vpos.xyz;
    vertexData.normal = N;
    vertexData.uv = uv;

    if(material.normalMap)
    {
        vec3 T = normalize(nmat * tangent);
        //T = normalize(T - dot(T, N) * N);
        vec3 B = normalize(nmat * bitangent);//normalize(cross(N, T));       
        mat3 tbn = mat3(T, B, N);
        vertexData.tbn = tbn;
    }    
}
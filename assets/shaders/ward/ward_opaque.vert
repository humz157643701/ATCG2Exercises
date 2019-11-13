#version 440 core
    
layout(location = 0) in vec3 position;
layout(location = 1) in vec2 uv;
layout(location = 2) in vec3 normal;
layout(location = 3) in vec3 tangent;

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
    
out struct VertexData
{
    vec3 viewspacePosition;
    vec2 uv;
    vec3 normal;
    mat3 tbn;
} vertexData;
    
void main()
{
    mat4 mmat = model_matrix;   
   
    mat4 modelview = camera.viewMatrix * mmat;
    vec4 vpos = modelview * vec4(position, 1.0f);
    gl_Position = camera.projectionMatrix * vpos;

    mat3 nmat = inverse(transpose(mat3(modelview)));
    vec3 N = normalize(nmat * normal);

    vertexData.viewspacePosition = vpos.xyz;
    vertexData.normal = N;
    vertexData.uv = uv;
 
    vec3 T = normalize(nmat * tangent);
    T = normalize(T - dot(T, N) * N);
    vec3 B = normalize(cross(N, T));       
    mat3 tbn = mat3(T, B, N);
    vertexData.tbn = tbn;  
}
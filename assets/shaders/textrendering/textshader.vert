#version 440 core
    
layout(location = 0) in vec2 quadv;
layout(location = 1) in vec4 pos;
layout(location = 2) in uint charindex;

layout(binding = 1, std140) uniform UVScale
{
    vec2 uvsc[128];
};

out vec3 uvl;
    
void main()
{
   gl_Position = vec4((quadv.x * pos.z + pos.x) * 2.0f - 1.0f, (quadv.y * pos.w + pos.y) * 2.0f - 1.0f, 0.0f, 1.0f);
   uvl = vec3(quadv.x * uvsc[charindex].x, (1.0 - quadv.y) * uvsc[charindex].y, float(charindex));
}
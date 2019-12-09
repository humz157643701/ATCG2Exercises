#version 440 core
     
in vec3 uvl;

uniform vec4 fontcolor;
uniform sampler2DArray glyphtex;

out vec4 color;

void main()
{
    color = vec4(fontcolor.rgb, texture(glyphtex, uvl).r * fontcolor.a);
}
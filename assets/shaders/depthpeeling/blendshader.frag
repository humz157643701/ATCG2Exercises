#version 440 core
//fragment shader output -----------------------------------------------------------------------------------------------
out vec4 color;    
in vec2 uv;
uniform sampler2D layer;
    
//main function --------------------------------------------------------------------------------------------------------
void main()
{
    vec4 layercol = texture(layer, uv);
    color = vec4(layercol.rgb * layercol.a, layercol.a);
}
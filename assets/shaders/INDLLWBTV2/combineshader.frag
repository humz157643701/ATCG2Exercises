#version 440 core
     
//fragment shader output -----------------------------------------------------------------------------------------------
out vec4 color;
    
in vec2 uv;

uniform sampler2D tbuf;
uniform sampler2D opqbuf;
uniform float tmwhite;

vec3 TM(vec3 l);
vec3 gammaCorrect(vec3 c);
float getLuminance(vec3 linearRGB);
    
//main function --------------------------------------------------------------------------------------------------------
void main()
{
    vec4 bg = texture(opqbuf, uv);
    vec4 fg = texture(tbuf, uv);

    vec4 outcol = vec4(fg.rgb + bg.rgb * fg.a, 1.0f);
    //tone mapping
    outcol.rgb = TM(outcol.rgb);    
    //gamma correction
    outcol.rgb = gammaCorrect(outcol.rgb);    
    //output all that stuff
    color = vec4(outcol);
}
    

vec3 TM(vec3 l)
{
    float il = getLuminance(l);
    float ol = (il * (1.0f + (il / (tmwhite * tmwhite))))/(1.0f + il);
    //float ol = il / (1.0f + il);
    return l * (ol / il);
}
    
vec3 gammaCorrect(vec3 c)
{
    return pow(c, vec3(1.0f / 2.2f));
}
    
float getLuminance(vec3 linearRGB)
{
    return dot(vec3(0.2126729,  0.7151522,  0.0721750), gammaCorrect(linearRGB));
}
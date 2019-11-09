#version 440 core
     
//define some constants
#define MAX_SHININESS 4000
#define INSULATOR_BASE_REFLECTIVITY 0.04f
#define PI 3.14159265359f
#define MAX_DIR_LIGHTS 2
#define MAX_POINT_LIGHTS 8
#define MAX_AMBIENT_LIGHTS 2

#define NUM_HIST_BINS 32
#define DEPTH_PADDING 1e-3

//important: force early depth test
layout(early_fragment_tests) in;
    
    
//some structs for cleaner code ----------------------------------------------------------------------------------------
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
    
struct DirectionalLight
{
    float lumint;
    vec3 color;
    vec3 direction;
};
    
struct PointLight
{
    float lumint;
    vec3 color;
    vec3 position;
};

struct AmbientLight
{
    float lumint;
    vec3 color;
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

struct Fragment
{
    vec4 color;
    float depth;
    int next;
};

in struct VertexData
{
    vec3 viewspacePosition;
    vec2 uv;
    vec3 normal;
    mat3 tbn;
} vertexData;
    
//uniforms -------------------------------------------------------------------------------------------------------------
//camera
uniform Camera camera;
    
//lights
uniform DirectionalLight dirlights[MAX_DIR_LIGHTS];
uniform int dirlightcount;

uniform PointLight pointlights[MAX_POINT_LIGHTS];
uniform int pointlightcount;

uniform AmbientLight ambientlights[MAX_AMBIENT_LIGHTS];
uniform int ambientlightcount;
    
//material
uniform Material material;
uniform int normalsfromtex;

uniform int usedepthfunc;

//linked list stuff
layout(binding = 0) uniform atomic_uint ac_fragcount;
uniform uint max_fragments;
layout(r32i, binding = 0) uniform coherent restrict iimage2D im_offsetbuffer;
layout(std430, binding = 0) buffer layoutFragList
{
  Fragment fragments[]; //stride is 32 bytes
} sb_fraglist;

//particles special case
uniform bool isParticle = false;
uniform float particleExponent = 1.0;
in float particleLife;

// stuff for depth histogram
uniform sampler2D depth_map;
out layout(location = 0) vec4 depth_hist0;
out layout(location = 1) vec4 depth_hist1;
out layout(location = 2) vec4 depth_hist2;
out layout(location = 3) vec4 depth_hist3;
out layout(location = 4) vec4 depth_hist4;
out layout(location = 5) vec4 depth_hist5;
out layout(location = 6) vec4 depth_hist6;
out layout(location = 7) vec4 depth_hist7;
    
//forward declaration of some helper functions -------------------------------------------------------------------------

vec3 readNormal(vec2 tc, vec3 vnorm);

//retrieve light values
vec3 calculateDirectionalLightIntensity(int i);
vec3 calculatePointLightIntensity(vec3 P, int i);
vec3 calcAmbientLightIntensity(int i);
    
//fresnel approximation
vec3 fresnelSchlick(vec3 L, vec3 H, vec3 F0);
    
//functions applying the shading model
vec3 shade(vec3 N, vec3 L, vec3 V, vec3 diffc, float shn, vec3 F0, float metallic);
vec3 ambientShade(vec3 diffc, vec3 F0, float metallic);
float specToShininess(float spec);
    
//tone mapping and gamma correction
vec3 inverseGammaCorrect(vec3 c);

//depth histogram stuffs
void fillHistogram(float vsp_depth);
    
//main function --------------------------------------------------------------------------------------------------------
void main()
{    
    //prepare vectors
    vec3 P = vertexData.viewspacePosition;    
    vec3 V = normalize(-P);
    vec2 tc = vertexData.uv;

    vec3 outcol = vec3(0.0f);  
   
    //read material properties
    vec4 matdiffdata = material.constDiff ? material.cdiff : texture(material.diff, tc);
    vec3 mat_diffuse = inverseGammaCorrect(matdiffdata.rgb);
    if(isParticle)
    {
        mat_diffuse *= material.cdiff.rgb;
    }
    float mat_alpha = matdiffdata.a;
    if(isParticle)
    {
        mat_alpha = pow(particleLife, particleExponent) * mat_alpha;
    }
    float mat_shininess = material.constShin ? specToShininess(material.cshin) : specToShininess(texture(material.shin, tc).r);
    vec3 N = material.normalMap ? readNormal(tc, normalize(vertexData.normal)) : normalize(vertexData.normal); 
    //color = vec4(N, 1.0f); return;
    vec3 F0 = material.constSpec ? material.cspec.rgb : texture(material.spec, tc).rgb;
    if(isParticle)
        F0 = vec3(0.0f, 0.0f, 0.0f);
    float mat_metallic = 1.0f - (material.constSpec ? material.cspec.a : texture(material.spec, tc).a);
    if(isParticle)
        mat_metallic = 0.0f;

    //directional lights
    for(int i = 0; i < dirlightcount; ++i)
    {
        vec3 L = normalize(-dirlights[i].direction);
        if(isParticle)
            N = L;
        vec3 intensity = calculateDirectionalLightIntensity(i);
        outcol += shade(N, L, V, mat_diffuse, mat_shininess, F0, mat_metallic) * intensity;
    }

    //point lights
    for(int i = 0; i < pointlightcount; ++i)
    {
        vec3 L = normalize(pointlights[i].position - P);
        if(isParticle)
            N = L;
        vec3 intensity = calculatePointLightIntensity(P, i);
        outcol += shade(N, L, V, mat_diffuse, mat_shininess, F0, mat_metallic) * intensity;
    }

    //ambient lights
    for(int i = 0; i < ambientlightcount; ++i)
    {
        outcol += ambientShade(mat_diffuse, F0, mat_metallic) * calcAmbientLightIntensity(i);
    }

    // fill depth histogram
    fillHistogram(vertexData.viewspacePosition.z);
    
    //get next free storage position
    uint slot = atomicCounterIncrement(ac_fragcount);
    if(slot >= max_fragments) //do nothing if buffer is full (the counter actually counts all shaded fragments. keep that in mind for possible dynamic operation!)
    {
        discard;
    }
    else
    {
        sb_fraglist.fragments[slot].color = vec4(outcol, mat_alpha);
        sb_fraglist.fragments[slot].depth = vertexData.viewspacePosition.z;
        sb_fraglist.fragments[slot].next = imageAtomicExchange(im_offsetbuffer, ivec2(gl_FragCoord.xy), int(slot));
    }
}
    
//retrieve light values
vec3 calculateDirectionalLightIntensity(int i)
{
    return dirlights[i].lumint * dirlights[i].color;
}
    
vec3 calculatePointLightIntensity(vec3 fpos, int i)
{
    vec3 d = fpos - pointlights[i].position;
    return (pointlights[i].lumint * pointlights[i].color) * (1.0f / dot(d, d));
}
    
vec3 calcAmbientLightIntensity(int i)
{
    return ambientlights[i].lumint * ambientlights[i].color;
}

vec3 fresnelSchlick(vec3 L, vec3 H, vec3 F0)
{
    return F0 + (1.0 - F0) * pow(1.0 - max(dot(L, H), 0.0f), 5.0);
}
 
vec3 shade(vec3 N, vec3 L, vec3 V, vec3 diffc, float shn, vec3 F0, float metallic)
{
    vec3 H = normalize(V + L);
    float NdotL = max(dot(N, L), 0.0f);
    float NdotH = max(dot(N, H), 0.0f);
    vec3 F;
    if(!isParticle)
        F = fresnelSchlick(L, H, vec3(F0));
    else
        F = vec3(0.0f, 0.0f, 0.0f);
    return ((vec3(1.0f) - F) * (diffc / PI) * (1.0f - metallic) + F * pow(NdotH, shn) * ((shn + 8) / (8 * PI))) * NdotL;
}
 
vec3 ambientShade(vec3 diffc, vec3 F0, float metallic)
{
    vec3 ambientColor = diffc * vec3(1.0f - F0) * (1.0f - metallic);
    return ambientColor;
}
    
float specToShininess(float spec)
{
    const float f = sqrt(MAX_SHININESS - 1.0f);
    return  1.0f + (f * spec) * (f * spec);
}
    
vec3 inverseGammaCorrect(vec3 c)
{
    return pow(c, vec3(2.2f));
}

vec3 readNormal(vec2 tc, vec3 vnorm)
{
    if(material.nmisbm)
    {
        float height = texture(material.norm, tc).r;
        vec2 hgrad = vec2(dFdx(height), dFdy(height));
        vec3 normal = vec3(-hgrad, 1.0f);
        return normalize(vertexData.tbn * normal);
    }
    else
    {
        return normalize(vertexData.tbn * (texture(material.norm, tc).xyz * 2.0f - 1.0f));
    }
}

void fillHistogram(float vsp_depth)
{
    vec2 ddepth = texelFetch(depth_map, ivec2(gl_FragCoord.xy), 0).rg;
    ddepth.g *= -1.0;
    ddepth.r += DEPTH_PADDING;
    ddepth.g -= DEPTH_PADDING;
    float binwidth = (ddepth.r - ddepth.g) / NUM_HIST_BINS;
    // calculate bin index
    int bidx = clamp(int((ddepth.r - vsp_depth) / binwidth), 0, int(NUM_HIST_BINS - 1));

    // to avoid branches we set the contribution for all histogram bins
    depth_hist0.r = float(bidx == 0);
    depth_hist0.g = float(bidx == 1);
    depth_hist0.b = float(bidx == 2);
    depth_hist0.a = float(bidx == 3);
 
    depth_hist1.r = float(bidx == 4);
    depth_hist1.g = float(bidx == 5);
    depth_hist1.b = float(bidx == 6);
    depth_hist1.a = float(bidx == 7);
 
    depth_hist2.r = float(bidx == 8);
    depth_hist2.g = float(bidx == 9);
    depth_hist2.b = float(bidx == 10);
    depth_hist2.a = float(bidx == 11);

    depth_hist3.r = float(bidx == 12);
    depth_hist3.g = float(bidx == 13); 
    depth_hist3.b = float(bidx == 14);
    depth_hist3.a = float(bidx == 15);

    depth_hist4.r = float(bidx == 16);
    depth_hist4.g = float(bidx == 17);
    depth_hist4.b = float(bidx == 18);
    depth_hist4.a = float(bidx == 19);

    depth_hist5.r = float(bidx == 20);
    depth_hist5.g = float(bidx == 21);
    depth_hist5.b = float(bidx == 22);
    depth_hist5.a = float(bidx == 23);

    depth_hist6.r = float(bidx == 24);
    depth_hist6.g = float(bidx == 25);
    depth_hist6.b = float(bidx == 26);
    depth_hist6.a = float(bidx == 27);

    depth_hist7.r = float(bidx == 28);
    depth_hist7.g = float(bidx == 29);
    depth_hist7.b = float(bidx == 30);
    depth_hist7.a = float(bidx == 31);
}
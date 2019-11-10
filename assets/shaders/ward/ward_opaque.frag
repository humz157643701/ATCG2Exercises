#version 440 core
     
//define some constants
#define MAX_SHININESS 4000
#define INSULATOR_BASE_REFLECTIVITY 0.04f
#define PI 3.14159265359f
#define MAX_DIR_LIGHTS 2
#define MAX_POINT_LIGHTS 8
#define MAX_AMBIENT_LIGHTS 2
    
    
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

    vec2 anisotropicRoughness;
    float anisotropicRotation;
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

//tone mapping
uniform float tmwhite;

//fragment shader output -----------------------------------------------------------------------------------------------
out vec4 color;
    
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
vec3 TM(vec3 l);
vec3 gammaCorrect(vec3 c);
vec3 inverseGammaCorrect(vec3 c);
    
//color space conversion
vec3 RGBtoXYZ(vec3 rgb);
vec3 XYZtoRGB(vec3 xyz);
float getLuminance(vec3 linearRGB);
    
//main function --------------------------------------------------------------------------------------------------------
void main()
{
    //prepare vectors
    vec3 P = vertexData.viewspacePosition;    
    vec3 V = normalize(-P);
    
    vec3 outcol = vec3(0.0f);    
   
    //prepare texture coordinate
    vec2 tc = vertexData.uv;
    //read material properties
    vec4 matdiffdata = material.constDiff ? material.cdiff : texture(material.diff, tc);
    if(matdiffdata.a < 0.001f)
        discard;
    float mat_diffuse = 1;// TODO diffuse coefficient                        inverseGammaCorrect(matdiffdata.rgb);

    float mat_spec = 1;// TODO specular coefficient                          material.constShin ? specToShininess(material.cshin) : specToShininess(texture(material.shin, tc).r);

    vec3 N = material.normalMap ? readNormal(tc, normalize(vertexData.normal)) : normalize(vertexData.normal); 
    //color = vec4(N, 1.0f); return;
    float F0 = 1; // TODO total hemispherical reflectance in normal direction                               material.constSpec ? material.cspec.rgb : texture(material.spec, tc).rgb;


    float mat_metallic = 1.0f - (material.constSpec ? material.cspec.a : texture(material.spec, tc).a);


    //directional lights
    for(int i = 0; i < dirlightcount; ++i)
    {       
        vec3 L = normalize(-dirlights[i].direction);
        vec3 intensity = calculateDirectionalLightIntensity(i);
        outcol += shade(N, L, V, mat_diffuse, mat_spec, F0, material.anisotropicRoughness, material.anisotropicRotation) * intensity;
    }

    //point lights
    for(int i = 0; i < pointlightcount; ++i)
    {
        vec3 L = normalize(pointlights[i].position - P);
        vec3 intensity = calculatePointLightIntensity(P, i);
        outcol += shade(N, L, V, mat_diffuse, map_spec, F0, material.anisotropicRoughness, material.anisotropicRotation) * intensity;
    }

    //ambient lights
    for(int i = 0; i < ambientlightcount; ++i)
    {
        outcol += ambientShade(mat_diffuse, F0, mat_metallic) * calcAmbientLightIntensity(i); // TODO <- what about this?
    }
    
    //tone mapping
    outcol = TM(outcol);    
    //gamma correction
    outcol = gammaCorrect(outcol);    
    //output all that stuff
    color = vec4(outcol, 1.0f);
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

vec3 fresnelSchlick(vec3 V, vec3 H, vec3 F0)
{
    vec3 cos_thetha_v = dot(H, V);
    return F0 + (1.0 - F0) * pow(1.0 - max(cos_thetha_v, 0.0f), 5.0);
}

mat4 rotationZ( in float angle ) {
	return mat4(	cos(angle),		-sin(angle),	0,	0,
			 		sin(angle),		cos(angle),		0,	0,
							0,				0,		1,	0,
							0,				0,		0,	1);
}
 
vec3 shade(vec3 N, vec3 L, vec3 V, float diffc, float specc, float F0, vec2 alpha, float phi_alpha)
{
    vec3 H_hat = V + L;
    vec3 H = normalize(H_hat);
    float NdotL = max(dot(N, L), 0.0f);
    float NdotH = max(dot(N, H), 0.0f);
    vec3 F;
    if(!isParticle)
        F = fresnelSchlick(V, H, vec3(F0));
    else
        F = vec3(0.0f, 0.0f, 0.0f);
    H_prime = rotationZ(phi_alpha)*H_hat;
    k_d = diffc/PI;
    k_s = specc*(1/(PI*alpha.x*alpha.y))*(1/(4*pow(normalize(dot(L, H_prime)), 2)*pow(H_prime.z, 4)))*exp(-(pow(H_prime.x/alpha.x, 2)+pow(H_prime.y/alpha.y, 2))/pow(H_prime.z,2));
    return k_d+(k_s*F) *NdotL;
    //return ((vec3(1.0f) - F) * (diffc / PI) * (1.0f - metallic) + F * pow(NdotH, shn) * ((shn + 8) / (8 * PI))) * NdotL;
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
    
vec3 inverseGammaCorrect(vec3 c)
{
    return pow(c, vec3(2.2f));
}
    
float getLuminance(vec3 linearRGB)
{
    return dot(vec3(0.2126729,  0.7151522,  0.0721750), gammaCorrect(linearRGB));
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
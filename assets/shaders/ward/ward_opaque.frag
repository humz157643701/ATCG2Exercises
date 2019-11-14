#version 440 core
     
//define some constants
#define MAX_SHININESS 4000
#define INSULATOR_BASE_REFLECTIVITY 0.04
#define PI 3.14159265359
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
    vec3 color;
    vec3 direction;
};
    
struct PointLight
{
    vec3 color;
    vec3 position;
};

struct AmbientLight
{
    vec3 color;
};
    
struct Material
{
    sampler2D s_diffuse_albedo;
    sampler2D s_normals;
    sampler2D s_specular_albedo;
    sampler2D s_roughness;
    sampler2D s_aniso_rotation;
    sampler2D s_fresnel_f0;
    sampler2D s_displacement;
    sampler2D s_transparency;
};

in struct VertexData
{
    vec3 viewspacePosition;
    vec2 uv;
    vec3 normal;
    vec3 tangent;
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

//tone mapping
uniform float tmwhite;

//fragment shader output -----------------------------------------------------------------------------------------------
out vec4 color;
    
//forward declaration of some helper functions -------------------------------------------------------------------------

//calculate incident radiances
vec3 Li_directional_light(int i);
vec3 Li_point_light(vec3 P, int i);
vec3 Li_ambient_light(int i);
    
//fresnel approximation
float fresnelSchlick(vec3 H, vec3 V, vec3 F0);
    
//functions applying the shading model
vec3 brdf(
    vec3 tsp_L,
    vec3 tsp_V,
    vec3 diffuse_albedo,
    vec3 specular_albedo,
    vec2 roughness,
    float aniso_rot,
    float fresnel_f0,
    float displacement,
    float transparency
);

vec3 ambientShade(vec3 diffuse_albedo, float fresnel_f0);
    
//tone mapping and gamma correction
vec3 TM(vec3 l);
vec3 gammaCorrect(vec3 c);
vec3 inverseGammaCorrect(vec3 c);
    
//color space conversion
vec3 RGBtoXYZ(vec3 rgb);
vec3 XYZtoRGB(vec3 xyz);
float getLuminance(vec3 linearRGB);

mat3 rotationZ(float angle);
    
//main function --------------------------------------------------------------------------------------------------------
void main()
{
    // initialize output radiance
    vec3 outcol = vec3(0.0f); 

    // tangent space defined by mesh geometry
    vec3 vsp_vertex_normal = normalize(vertexData.normal);
    vec3 vsp_vertex_tangent = normalize(vertexData.tangent);
    // gram-schmidt orthogonalize tangent, build orthonormal tangent basis
    vsp_vertex_tangent = normalize(vsp_vertex_tangent - dot(vsp_vertex_tangent, vsp_vertex_normal) * vsp_vertex_normal);
    vec3 vsp_vertex_bitangent = normalize(cross(vsp_vertex_normal, vsp_vertex_tangent));
    vsp_vertex_bitangent *= sign(dot(vsp_vertex_bitangent, vsp_vertex_normal));
    // tangent to view space
    mat3 tbn = mat3(vsp_vertex_tangent, vsp_vertex_bitangent, vsp_vertex_normal);
    // view to tangent space
    mat3 itbn = transpose(tbn);       
   
    //prepare texture coordinate
    vec2 tc = vertexData.uv;

    //read material properties
    vec3 diffuse_albedo     = texture(material.s_diffuse_albedo, tc).rgb;
    vec3 specular_albedo    = texture(material.s_specular_albedo, tc).rgb;
    vec3 vtsp_bump_normal    = normalize(texture(material.s_normals, tc).rgb);
    vec2 roughness          = texture(material.s_roughness, tc).rg;
    float aniso_rotation    = texture(material.s_aniso_rotation, tc).r;
    float fresnel_f0        = texture(material.s_fresnel_f0, tc).r;
    float displacement      = texture(material.s_displacement, tc).r;
    float transparency      = texture(material.s_transparency, tc).r;

    // we have to perturb the tangent frame defined by the mesh geometry using the normal from the normal map
    vec3 vsp_bump_normal = tbn * vtsp_bump_normal;
    // gram-schmidt orthogonalize the geometry tangent
    vec3 vsp_bump_tangent = normalize(vsp_vertex_tangent - dot(vsp_vertex_tangent, vsp_bump_normal) * vsp_bump_normal);
    // calculate bitangent
    vec3 vsp_bump_bitangent = normalize(cross(vsp_bump_normal, vsp_bump_tangent));
    vsp_bump_bitangent *= sign(dot(vsp_bump_bitangent, vsp_bump_normal));
    // change of basis matrices for perturbed tangent space
    mat3 bump_tbn = mat3(vsp_bump_tangent, vsp_bump_bitangent, vsp_bump_normal);
    mat3 bump_itbn = transpose(bump_tbn);

    //prepare vectors
    vec3 vsp_P = vertexData.viewspacePosition;    
    vec3 vsp_V = normalize(-vsp_P);
    vec3 tsp_V = itbn * vsp_V;

    //directional lights
    for(int i = 0; i < dirlightcount; ++i)
    {   
        vec3 vsp_L = normalize(-dirlights[i].direction);
        vec3 tsp_L = bump_itbn * vsp_L;
        vec3 Li = Li_directional_light(i);
        outcol += brdf(tsp_L, tsp_V, diffuse_albedo, specular_albedo, roughness, aniso_rotation, fresnel_f0, displacement, transparency) * max(dot(vsp_bump_normal, vsp_L), 0.0) * Li;
    }

    //point lights
    for(int i = 0; i < pointlightcount; ++i)
    {
        vec3 vsp_L = normalize(pointlights[i].position - vsp_P);
        vec3 tsp_L = bump_itbn * vsp_L;
        vec3 Li = Li_point_light(vsp_P, i);
        outcol += brdf(tsp_L, tsp_V, diffuse_albedo, specular_albedo, roughness, aniso_rotation, fresnel_f0, displacement, transparency) * max(dot(vsp_bump_normal, vsp_L), 0.0) * Li;
    }

    //ambient lights
    for(int i = 0; i < ambientlightcount; ++i)
    {
        outcol += ambientShade(diffuse_albedo, fresnel_f0) * Li_ambient_light(i); // TODO <- what about this?
    }
    
    //tone mapping
    outcol = TM(outcol);    
    //gamma correction
    outcol = gammaCorrect(outcol);    
    //output all that stuff
    color = vec4(outcol, 1.0f);
}
    
//retrieve light values
vec3 Li_directional_light(int i)
{
    return dirlights[i].color;
}
    
vec3 Li_point_light(vec3 fpos, int i)
{
    vec3 d = fpos - pointlights[i].position;
    return (pointlights[i].color) * (1.0f / dot(d, d));
}
    
vec3 Li_ambient_light(int i)
{
    return ambientlights[i].color;
}

float fresnelSchlick(vec3 H, vec3 V, float F0)
{
    float cos_thetha_v = dot(H, V);
    return F0 + (1.0 - F0) * pow(1.0 - max(cos_thetha_v, 0.0f), 5.0);
}

mat3 rotationZ( float angle ) {
	return mat3(	cos(angle),		-sin(angle),	0.0,
			 		sin(angle),		cos(angle),		0.0,
							0.0,				0.0,		1.0);
}
 
vec3 brdf(
    vec3 tsp_L,
    vec3 tsp_V,
    vec3 diffuse_albedo,
    vec3 specular_albedo,
    vec2 roughness,
    float aniso_rot,
    float fresnel_f0,
    float displacement,
    float transparency
)
{
    // diffuse part
    vec3 kd = diffuse_albedo / PI;

    // specular part
    vec3 H_hat = tsp_L + tsp_V;
    vec3 H = normalize(H_hat);
    vec3 H_prime = rotationZ(aniso_rot) * H_hat;
    float LdotH_prime_2 = abs(dot(tsp_L, H_prime)) * abs(dot(tsp_L, H_prime));    
    float Hz_prime_2 = H_prime.z * H_prime.z;
    float Hz_prime_4 = Hz_prime_2 * Hz_prime_2;
    float Hx_prime_alphax_2 = (H_prime.x / roughness.x) * (H_prime.x / roughness.x);
    float Hy_prime_alphay_2 = (H_prime.y / roughness.y) * (H_prime.y / roughness.y);
    vec3 ks = specular_albedo * (1.0 / (PI * roughness.x * roughness.y)) * (1.0 / (4.0 * LdotH_prime_2 * Hz_prime_4)) * exp(-((Hx_prime_alphax_2 + Hy_prime_alphay_2) / Hz_prime_2));
    // fresnel term
    float F = fresnelSchlick(H, tsp_V, fresnel_f0);
    
    return kd + ks * F;
}
 
vec3 ambientShade(vec3 diffuse_albedo, float fresnel_f0)
{
    return diffuse_albedo;
}

vec3 TM(vec3 l)
{
    float il = getLuminance(l);
    float ol = (il * (1.0f + (il / (tmwhite * tmwhite))))/(1.0f + il);
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
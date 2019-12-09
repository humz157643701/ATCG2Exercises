#version 430 core
     
//define some constants
#define PI 3.14159265359
#define MAX_DIR_LIGHTS 2
#define MAX_POINT_LIGHTS 8
#define MAX_AMBIENT_LIGHTS 2
#define SPECULAR_IBL_SAMPLES 50
#define SPECULAR_IBL_LOD_BIAS_POW 0.7
    
    
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
    vec2 f_mscale;
};

in struct VertexData
{
    vec3 viewspacePosition;
    vec2 uv;
    vec3 normal;
    vec3 tangent;
} vertexData;

flat in mat3 modelrot;
    
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
    
uniform samplerCube skybox;
uniform samplerCube irradiance;
uniform float skybox_res;
uniform float skybox_lodlevels;
//material
uniform Material material;

//tone mapping
uniform float exposure;

//fragment shader output -----------------------------------------------------------------------------------------------
out vec4 color;
    
//forward declaration of some helper functions -------------------------------------------------------------------------

//calculate incident radiances
vec3 Li_directional_light(int i);
vec3 Li_point_light(vec3 P, int i);
vec3 Li_ambient_light(int i);
    
//fresnel approximation
float fresnelSchlick(vec3 H, vec3 V, float F0);
float fresnelSchlickRoughness(vec3 H, vec3 V, float F0, float roughness);
    
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

vec3 brdf_s(
    vec3 tsp_L,
    vec3 tsp_V,
    vec3 specular_albedo,
    vec2 roughness,
    float aniso_rot,
    float fresnel_f0,
    float displacement,
    float transparency
);

vec3 ambientShade(vec3 diffuse_albedo, float fresnel_f0, vec3 N, vec3 L, float roughness);
    
//tone mapping and gamma correction
vec3 TM_Exposure(vec3 l);
vec3 gammaCorrect(vec3 c);
vec3 inverseGammaCorrect(vec3 c);
    
float getLuminance(vec3 linearRGB);

mat3 rotationZ(float angle);

float RadicalInverse_VdC(uint bits) 
{
    bits = (bits << 16u) | (bits >> 16u);
    bits = ((bits & 0x55555555u) << 1u) | ((bits & 0xAAAAAAAAu) >> 1u);
    bits = ((bits & 0x33333333u) << 2u) | ((bits & 0xCCCCCCCCu) >> 2u);
    bits = ((bits & 0x0F0F0F0Fu) << 4u) | ((bits & 0xF0F0F0F0u) >> 4u);
    bits = ((bits & 0x00FF00FFu) << 8u) | ((bits & 0xFF00FF00u) >> 8u);
    return float(bits) * 2.3283064365386963e-10; // / 0x100000000
}

vec2 Hammersley(uint i, uint N)
{
    return vec2(float(i)/float(N), RadicalInverse_VdC(i));
}  
    
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
    vsp_vertex_bitangent *= sign(dot(cross(vsp_vertex_tangent, vsp_vertex_bitangent), vsp_vertex_normal));
    // tangent to view space
    mat3 tbn = mat3(vsp_vertex_tangent, vsp_vertex_bitangent, vsp_vertex_normal);
    // view to tangent space
    mat3 itbn = transpose(tbn);     
   
    //prepare texture coordinate
    vec2 tc = vertexData.uv * material.f_mscale;

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
    vsp_bump_bitangent *= sign(dot(cross(vsp_bump_tangent, vsp_bump_bitangent), vsp_bump_normal));
    // change of basis matrices for perturbed tangent space
    mat3 bump_tbn = mat3(vsp_bump_tangent, vsp_bump_bitangent, vsp_bump_normal);
    mat3 bump_itbn = transpose(bump_tbn);
    mat3 R = rotationZ(aniso_rotation);
    bump_itbn = R * bump_itbn;
    bump_tbn = transpose(bump_itbn);

    //prepare vectors
    vec3 vsp_P = vertexData.viewspacePosition;    
    vec3 vsp_V = normalize(-vsp_P);
    vec3 tsp_V = bump_itbn * vsp_V;

    // transformation from perturbed tangent space to world space (direction only)
    mat3 bump_tangent_to_world = transpose(mat3(camera.viewMatrix)) * bump_tbn;

    //directional lights
    for(int i = 0; i < dirlightcount; ++i)
    {   
        vec3 vsp_L = normalize(-dirlights[i].direction);
        vec3 tsp_L = bump_itbn * vsp_L;
        vec3 Li = Li_directional_light(i);
        outcol += brdf(tsp_L, tsp_V, diffuse_albedo, specular_albedo, roughness, aniso_rotation, fresnel_f0, displacement, transparency) *
                  max(dot(vsp_bump_normal, vsp_L), 0.0) *
                  Li;
    }

    //point lights
    for(int i = 0; i < pointlightcount; ++i)
    {
        vec3 vsp_L = normalize(pointlights[i].position - vsp_P);
        vec3 tsp_L = bump_itbn * vsp_L;
        vec3 Li = Li_point_light(vsp_P, i);
        outcol += brdf(tsp_L, tsp_V, diffuse_albedo, specular_albedo, roughness, aniso_rotation, fresnel_f0, displacement, transparency) *
                  max(dot(vsp_bump_normal, vsp_L), 0.0) *
                  Li;
    }

    outcol += ambientShade(diffuse_albedo,
                        fresnel_f0,
                        vsp_bump_normal,
                        vsp_V,
                        max(roughness.x, roughness.y)) * texture(irradiance, bump_tangent_to_world * vec3(0.0, 0.0, 1.0)).rgb;

    // specular IBL
    vec3 sibl = vec3(0.0);

    // from bump tangent space to world space
    // idea: choose sample count based on roughness
    for(int i = 0; i < SPECULAR_IBL_SAMPLES; ++i)
    {
        vec2 usamp = Hammersley(i, SPECULAR_IBL_SAMPLES);
        float hphi = atan(roughness.y * sin(2.0 * PI * usamp.y), roughness.x * cos(2.0 * PI * usamp.y));
        float htheta = atan(sqrt((-log(1.0 - usamp.x)) / (((cos(hphi) * cos(hphi)) / (roughness.x * roughness.x)) + ((sin(hphi) * sin(hphi)) / (roughness.y * roughness.y)))));
        vec3 h = vec3(
            sin(htheta) * cos(hphi),
            sin(htheta) * sin(hphi),
            cos(htheta)
        );

        vec3 tsp_L = 2.0 * dot(tsp_V, h) * h - tsp_V;
        vec3 wsp_L = bump_tangent_to_world * tsp_L;

        //vec3 b = brdf_s(tsp_L, tsp_V, specular_albedo, roughness, aniso_rotation, fresnel_f0, displacement, transparency);
        float w = 2.0 / (1.0 + (tsp_V.z / tsp_L.z));
        float F = fresnelSchlick(h, tsp_V, fresnel_f0);

        //vec3 pdf = (1.0 / (F * specular_albedo * w));
        vec3 Li = texture(skybox, wsp_L, pow(max(roughness.x, roughness.y), SPECULAR_IBL_LOD_BIAS_POW) * skybox_lodlevels).rgb;
        sibl += Li * max(tsp_L.z, 0.0) * F * specular_albedo * w;        
    }
    outcol += sibl / float(SPECULAR_IBL_SAMPLES);
    
    //tone mapping
    outcol = TM_Exposure(outcol);
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

float fresnelSchlickRoughness(vec3 H, vec3 V, float F0, float roughness)
{
    float cosTheta = max(dot(H, V), 0.0);
    return F0 + (max(1.0 - roughness, F0) - F0) * pow(1.0 - cosTheta, 5.0);
}

float fresnelSchlick(vec3 H, vec3 V, float F0)
{
    float cos_thetha_v = max(dot(H, V), 0.0);
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
    vec3 H = normalize(tsp_L + tsp_V);    
    float LdotH_prime_2 = max(dot(tsp_L, H), 0.0) * max(dot(tsp_L, H), 0.0);    
    float Hz_prime_2 = H.z * H.z;
    float Hz_prime_4 = Hz_prime_2 * Hz_prime_2;
    float Hx_prime_alphax_2 = (H.x / roughness.x) * (H.x / roughness.x);
    float Hy_prime_alphay_2 = (H.y / roughness.y) * (H.y / roughness.y);
    vec3 ks = specular_albedo * (1.0 / (PI * roughness.x * roughness.y)) * (1.0 / (4.0 * LdotH_prime_2 * Hz_prime_4)) * exp(-((Hx_prime_alphax_2 + Hy_prime_alphay_2) / Hz_prime_2));
    // fresnel term
    float F = fresnelSchlick(H, tsp_V, fresnel_f0);
    return kd + ks * F;//(1.0 - F) * kd + ks * F;
}

vec3 brdf_s(
    vec3 tsp_L,
    vec3 tsp_V,
    vec3 specular_albedo,
    vec2 roughness,
    float aniso_rot,
    float fresnel_f0,
    float displacement,
    float transparency
)
{
    // specular part
    vec3 H = normalize(tsp_L + tsp_V);    
    float LdotH_prime_2 = max(dot(tsp_L, H), 0.0) * max(dot(tsp_L, H), 0.0);    
    float Hz_prime_2 = H.z * H.z;
    float Hz_prime_4 = Hz_prime_2 * Hz_prime_2;
    float Hx_prime_alphax_2 = (H.x / roughness.x) * (H.x / roughness.x);
    float Hy_prime_alphay_2 = (H.y / roughness.y) * (H.y / roughness.y);
    vec3 ks = specular_albedo * (1.0 / (PI * roughness.x * roughness.y)) * (1.0 / (4.0 * LdotH_prime_2 * Hz_prime_4)) * exp(-((Hx_prime_alphax_2 + Hy_prime_alphay_2) / Hz_prime_2));
    // fresnel term
    float F = fresnelSchlick(H, tsp_V, fresnel_f0);
    
    return ks * F;
}
 
vec3 ambientShade(vec3 diffuse_albedo, float fresnel_f0, vec3 N, vec3 V, float roughness)
{
    float f0 = mix(float((1.0 - fresnel_f0) < 1e-6), 0.04, fresnel_f0);       
    return (diffuse_albedo / PI) * (1.0 - fresnelSchlickRoughness(N, V, f0, roughness));
}

vec3 TM_Exposure(vec3 l)
{
    float il = getLuminance(l);
    float ol = 1.0 - exp(- (il * exposure));
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
    return dot(vec3(0.2126729,  0.7151522,  0.0721750), gammaCorrect(abs(linearRGB)));
}
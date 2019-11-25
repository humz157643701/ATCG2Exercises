#version 430 core
     
//define some constants
#define PI 3.14159265359
#define MAX_DIR_LIGHTS 2
#define MAX_POINT_LIGHTS 8
#define MAX_AMBIENT_LIGHTS 2
#define SPECULAR_IBL_SAMPLES 32
#define SPECULAR_IBL_LOD_BIAS_POW 0.85
#define DIELECTRIC_F0 0.04
    
    
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
vec3 fresnelSchlick(vec3 H, vec3 V, vec3 F0);
vec3 fresnelSchlickRoughness(vec3 H, vec3 V, vec3 F0, float roughness);
float GGX_NDF(vec3 N, vec3 H, float roughness);
float GeometrySmith(vec3 N, vec3 V, vec3 L, float roughness);
    
//functions applying the shading model
vec3 brdf(
    vec3 vsp_L,
    vec3 vsp_V,
    vec3 vsp_N,
    vec3 diffuse_albedo,
    vec3 specular_albedo,
    float roughness
);

vec3 brdf_s(
    vec3 vsp_L,
    vec3 vsp_V,
    vec3 vsp_N,
    vec3 diffuse_albedo,
    vec3 specular_albedo,
    vec2 roughness
);

vec3 ambientShade(vec3 diffuse_albedo, vec3 fresnel_f0, vec3 N, vec3 V, float roughness);
    
//tone mapping and gamma correction
vec3 TM_Exposure(vec3 l);
vec3 gammaCorrect(vec3 c);
vec3 inverseGammaCorrect(vec3 c);
    
float getLuminance(vec3 linearRGB);

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
    vec3 diffuse_albedo     = texture(material.s_diffuse_albedo, tc).rgb;//inverseGammaCorrect(texture(material.s_diffuse_albedo, tc).rgb);
    vec3 specular_albedo    = inverseGammaCorrect(texture(material.s_specular_albedo, tc).rgb);
    vec3 vtsp_bump_normal    = normalize(texture(material.s_normals, tc).rgb * 2.0 - 1.0);
    float roughness          = texture(material.s_roughness, tc).r ;
    float displacement      = texture(material.s_displacement, tc).r;
    float transparency      = texture(material.s_transparency, tc).r;

    // we have to perturb the tangent frame defined by the mesh geometry using the normal from the normal map
    vec3 vsp_bump_normal = tbn * vtsp_bump_normal;

    //prepare vectors
    vec3 vsp_P = vertexData.viewspacePosition;    
    vec3 vsp_V = normalize(-vsp_P);


    // tangent space to world space
    mat3 view_to_world = transpose(mat3(camera.viewMatrix));
    //vec3 vsp_R = 2.0 * dot(vsp_V, vsp_bump_normal) * vsp_bump_normal - vsp_V;
    //directional lights
    for(int i = 0; i < dirlightcount; ++i)
    {   
        vec3 vsp_L = normalize(-dirlights[i].direction);
        vec3 Li = Li_directional_light(i);
        outcol += brdf(vsp_L, vsp_V, vsp_bump_normal, diffuse_albedo, specular_albedo, roughness) *
                  max(dot(vsp_bump_normal, vsp_L), 0.0) *
                  Li;
    }

    //point lights
    for(int i = 0; i < pointlightcount; ++i)
    {
        vec3 vsp_L = normalize(pointlights[i].position - vsp_P);
        vec3 Li = Li_point_light(vsp_P, i);
        outcol += brdf(vsp_L, vsp_V, vsp_bump_normal, diffuse_albedo, specular_albedo, roughness) *
                  max(dot(vsp_bump_normal, vsp_L), 0.0) *
                  Li;
    }

    // diffuse ibl
    outcol += ambientShade(diffuse_albedo, specular_albedo, vsp_bump_normal, vsp_V, roughness) * texture(irradiance, view_to_world * vsp_bump_normal).rgb;
    
    // // specular IBL
    // vec3 sibl = vec3(0.0);

    // // idea: choose sample count based on roughness
    // float r2 = roughness * roughness;

    // vec3 up = vec3(0.0, 1.0, 0.0);
    // vec3 right = normalize(cross(up, vsp_bump_normal));
    // up = normalize(cross(vsp_bump_normal, right));

    // for(int i = 0; i < SPECULAR_IBL_SAMPLES; ++i)
    // {
    //     vec2 usamp = Hammersley(i, SPECULAR_IBL_SAMPLES);
        
    //     float htheta = acos(sqrt((1.0 - usamp.x) / (usamp.y * (r2 - 1.0) + 1.0)));
    //     float hphi = 2.0 * PI * usamp.x;       

    //     vec3 h = vec3(
    //         sin(htheta) * cos(hphi),            
    //         sin(htheta) * sin(hphi),
    //         cos(htheta),
    //     );

    //     h = tbn * h;        

    //     vec3 vsp_L = 2.0 * dot(vsp_V, h) * h - vsp_V;

    //     float G = GeometrySmith(vsp_bump_normal, vsp_V, vsp_L, roughness);
    //     vec3 F = fresnelSchlick(h, vsp_V, specular_albedo);

    //     float LdotH = max(dot(vsp_L, h), 0.0); 
    //     float NdotH = max(dot(vsp_bump_normal, h), 0.0); 
    //     float VdotN = max(dot(vsp_L, vsp_bump_normal), 0.0);
    //     float LdotN = max(dot(vsp_L, vsp_bump_normal), 0.0);

    //     vec3 b = (F * G * LdotH) / (VdotN * NdotH);
       
    //     vec3 wsp_L = view_to_world * vsp_L;
    //     vec3 Li = texture(skybox, wsp_L, pow(roughness, SPECULAR_IBL_LOD_BIAS_POW) * skybox_lodlevels).rgb;
    //     sibl += b * Li * LdotN;             
    // }
    // outcol += sibl / float(SPECULAR_IBL_SAMPLES);

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

vec3 fresnelSchlick(vec3 H, vec3 V, vec3 F0)
{
    float cosTheta = max(dot(H, V), 0.0);
    return F0 + (1.0 - F0) * pow(1.0 - cosTheta, 5.0);
}

vec3 fresnelSchlickRoughness(vec3 H, vec3 V, vec3 F0, float roughness)
{
    float cosTheta = max(dot(H, V), 0.0);
    return F0 + (max(vec3(1.0 - roughness), F0) - F0) * pow(1.0 - cosTheta, 5.0);
}

float GGX_NDF(vec3 N, vec3 H, float roughness)
{
    float r2 = roughness * roughness;
    float r4 = r2 * r2;
    float NdotH  = max(dot(N, H), 0.0);
    float NdotH2 = NdotH * NdotH;
    float d = (NdotH2 * (r4 - 1.0) + 1.0);	
    return r2 / (PI * d * d);
}

float GeometrySmith(vec3 N, vec3 V, vec3 L, float roughness)
{
    float NdotV = max(dot(N, V), 0.0);
    float NdotL = max(dot(N, L), 0.0);
    float r = (roughness + 1.0);
    float k = (r * r) / 8.0;	
    float G1 = NdotV / (NdotV * (1.0 - k) + k);
    float G2 = NdotL / (NdotL * (1.0 - k) + k);
    return G1 * G2;    
}
 
vec3 brdf(
    vec3 vsp_L,
    vec3 vsp_V,
    vec3 vsp_N,
    vec3 diffuse_albedo,
    vec3 specular_albedo,
    float roughness
)
{   
    // cook torrance with GGX NDF, Smith shadowing/masking and Fresnel-Schlick reflectance
    vec3 h = normalize(vsp_V + vsp_L);
    float D = GGX_NDF(vsp_N, h, roughness);
    float G = GeometrySmith(vsp_N, vsp_V, vsp_L, roughness);
    vec3 F = fresnelSchlick(h, vsp_V, specular_albedo);
   
    vec3 n = D * F * G;
    float d = 4.0 * max(dot(vsp_N, vsp_V), 0.0) * max(dot(vsp_N, vsp_L), 0.0);
    vec3 specular = n / max(d, 1e-6);

    vec3 diffuse = ((vec3(1.0) - F) * diffuse_albedo) / PI;

    return diffuse + specular;
}

vec3 brdf_s(
    vec3 vsp_L,
    vec3 vsp_V,
    vec3 vsp_N,
    vec3 diffuse_albedo,
    vec3 specular_albedo,
    float roughness
)
{
   // cook torrance with GGX NDF, Smith shadowing/masking and Fresnel-Schlick reflectance
    vec3 h = normalize(vsp_V + vsp_L);
    float D = GGX_NDF(vsp_N, h, roughness);
    float G = GeometrySmith(vsp_N, vsp_V, vsp_L, roughness);
    vec3 F = fresnelSchlick(h, vsp_V, specular_albedo);
   
    vec3 n = D * F * G;
    float d = 4.0 * max(dot(vsp_N, vsp_V), 0.0) * max(dot(vsp_N, vsp_L), 0.0);
    vec3 specular = n / max(d, 1e-6);

    return specular;
}
 
vec3 ambientShade(vec3 diffuse_albedo, vec3 fresnel_f0, vec3 N, vec3 V, float roughness)
{      
    return (diffuse_albedo / PI) * (vec3(1.0) - fresnelSchlickRoughness(N, V, fresnel_f0, roughness));
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
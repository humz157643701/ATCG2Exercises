#version 440 core
#define MAX_FRAGMENTS 1024

#define BACK_REV_THRESHOLD 0.0000001

#define NUM_HIST_BINS 32
#define DEPTH_PADDING 1e-3

#define MAX_BINS 64

struct Fragment
{
    vec4 color;
    float depth;
    int next;
};

in vec2 uv;
out vec4 color;

// fragment list
uniform isampler2D im_offsetbuffer;
layout(std430, binding = 0) buffer layoutFragList
{
  Fragment fragments[]; //stride is 32 bytes
} sb_fraglist;

// bin count
uniform uint numbins;

// dual depth man of transparent geometry
uniform sampler2D depth_map;

// depth weight offset
uniform float dwoffset;

// fragment count
uniform sampler2D frag_count;

// depth histogram
uniform sampler2D depth_hist0;
uniform sampler2D depth_hist1;
uniform sampler2D depth_hist2;
uniform sampler2D depth_hist3;
uniform sampler2D depth_hist4;
uniform sampler2D depth_hist5;
uniform sampler2D depth_hist6;
uniform sampler2D depth_hist7;

// color and rev per bin
vec4 color_bins[MAX_BINS];
// normalizer per bin
float color_normalizer[MAX_BINS];
// bin array of the cumulative depth histogram
float cumulated_depth_hist[NUM_HIST_BINS];

// function forward declarations ---------------------------------------------------------------------------------------
// WBOIT weighting function
//float depthWeight(float z, float a, float start, float end);
float depthWeight(float depth, float startdepth, float lambda);
// bin smoothing kernel
float binWeight(float z, float start, float end);
// cumulated depth histogram
void computeCumulatedDepthHist(ivec2 fc);
// calculates final color bin for some fragment at depth d
uint calculateBinIndex(float d, vec2 ddb);

uniform uint selbin;

//main function --------------------------------------------------------------------------------------------------------
void main()
{        
    ivec2 frag_coord = ivec2(gl_FragCoord.xy);
    uint fragCount = uint(texelFetch(frag_count, frag_coord, 0).r + 0.5);
    if(fragCount == 0)
    {
        color = vec4(0.0f, 0.0f, 0.0f, 1.0f);
        return;
    }
    // generate cumulative depth histogram
    computeCumulatedDepthHist(frag_coord);

    // clear bin contents
    for(int i = 0; i < numbins; ++i)
    {
        color_bins[i] = vec4(0.0, 0.0, 0.0, 1.0);
        color_normalizer[i] = 0.0;
    }

    // sample depth values from dual depth buffer
    vec2 ddepth = texelFetch(depth_map, ivec2(gl_FragCoord.xy), 0).rg;
    ddepth.g *= -1.0;
    ddepth.r += DEPTH_PADDING;
    ddepth.g -= DEPTH_PADDING;

    // iterate over all fragments
    int fidx = texelFetch(im_offsetbuffer, frag_coord, 0).r;

    // total background visibility
    float totalrev = 1.0f;

    // fetch start index for this pixel's fragment list
    int index = texelFetch(im_offsetbuffer, frag_coord, 0).r;
    while(index > -1)
    {
        Fragment f = sb_fraglist.fragments[index];        
        index = f.next; // traverse the list

        // calculate bin index
        uint binindex = calculateBinIndex(f.depth, ddepth);

        // if(binindex != selbin)
        //     continue;

        // calculate contribution to surrounding color bins
        float w = f.color.a * depthWeight(f.depth, ddepth.r, dwoffset);//depthWeight(f.depth, f.color.a, ddepth.r + dwoffset * (ddepth.r - ddepth.g), ddepth.g); // * binweight
        color_bins[binindex].rgb += f.color.rgb * w;
        color_bins[binindex].a *= (1.0 - f.color.a); // a * binweight
        color_normalizer[binindex] += w;

        // update background visibility
        totalrev *= (1.0 - f.color.a);
    }

    // normalize wboit results per bin, blend bins together, normalize
    vec3 acc_c = vec3(0.0, 0.0, 0.0);
    float a_dst = 1.0;
    float acc_normalizer = 0.0;

    for(int k = 0; k < numbins; ++k)
    {
        acc_c += (color_bins[k].rgb / max(color_normalizer[k], 1e-10)) * a_dst * (1.0 - color_bins[k].a);
        acc_normalizer += a_dst * (1.0 - color_bins[k].a);
        a_dst *= color_bins[k].a;
    }

    color = vec4((acc_c / max(acc_normalizer, 1e-10)) * (1.0 - totalrev), totalrev);
}

// float depthWeight(float z, float a, float start, float end)
// {
//     float absz = (abs(z) - abs(start)) * (500.0 / (abs(end) - abs(start)));
//     return a * max(0.01, min(3000.0, 10.0 / (0.00001 + (absz / 10.0) * (absz / 10.0) * (absz / 10.0) + (absz / 200.0) * (absz / 200.0) * (absz / 200.0) * (absz / 200.0) * (absz / 200.0) * (absz / 200.0))));
// }

float depthWeight(float depth, float startdepth, float lambda)
{
    return exp(- lambda * (startdepth - depth)) * (1.0 - 1e-5) + 1e-5;
}

float binWeight(float z, float start, float end)
{
    // tent function ?
    // gaussian with certain standard deviation ?
    // epanechnikov kernel ?
    return 1.0;
}

uint calculateBinIndex(float d, vec2 ddb)
{
    // calculate depth histogram index
    float bin_width = (ddb.x - ddb.y) / float(NUM_HIST_BINS);    
    int hidx = clamp(int((ddb.x - d) / bin_width), int(0), int(NUM_HIST_BINS - 1));
    // calculate color bin index
    // read bin index from cumulated scaled histogram
    float hist_n = cumulated_depth_hist[clamp(hidx - 1, 0, NUM_HIST_BINS - 1)] * float(sign(hidx)); // 0 if hidx = 0
    float hist_n1 = cumulated_depth_hist[hidx];
    uint cbidx = uint(hist_n + ((hist_n1 - hist_n) / bin_width) * (ddb.x - d - float(hidx) * bin_width));
    return clamp(cbidx, 0, numbins - 1);
}

void computeCumulatedDepthHist(ivec2 fc)
{
    vec4 dh0 = floor(texelFetch(depth_hist0, fc, 0) + 0.5);
    vec4 dh1 = floor(texelFetch(depth_hist1, fc, 0) + 0.5);
    vec4 dh2 = floor(texelFetch(depth_hist2, fc, 0) + 0.5);
    vec4 dh3 = floor(texelFetch(depth_hist3, fc, 0) + 0.5);
    vec4 dh4 = floor(texelFetch(depth_hist4, fc, 0) + 0.5);
    vec4 dh5 = floor(texelFetch(depth_hist5, fc, 0) + 0.5);
    vec4 dh6 = floor(texelFetch(depth_hist6, fc, 0) + 0.5);
    vec4 dh7 = floor(texelFetch(depth_hist7, fc, 0) + 0.5);

    cumulated_depth_hist[0] =                             dh0.r;
    cumulated_depth_hist[1] =  cumulated_depth_hist[0] +  dh0.g;
    cumulated_depth_hist[2] =  cumulated_depth_hist[1] +  dh0.b;
    cumulated_depth_hist[3] =  cumulated_depth_hist[2] +  dh0.a;

    cumulated_depth_hist[4] =  cumulated_depth_hist[3] +  dh1.r;
    cumulated_depth_hist[5] =  cumulated_depth_hist[4] +  dh1.g;
    cumulated_depth_hist[6] =  cumulated_depth_hist[5] +  dh1.b;
    cumulated_depth_hist[7] =  cumulated_depth_hist[6] +  dh1.a;

    cumulated_depth_hist[8] =  cumulated_depth_hist[7] +  dh2.r;
    cumulated_depth_hist[9] =  cumulated_depth_hist[8] +  dh2.g;
    cumulated_depth_hist[10] = cumulated_depth_hist[9] +  dh2.b;
    cumulated_depth_hist[11] = cumulated_depth_hist[10] + dh2.a;

    cumulated_depth_hist[12] = cumulated_depth_hist[11] + dh3.r;
    cumulated_depth_hist[13] = cumulated_depth_hist[12] + dh3.g;
    cumulated_depth_hist[14] = cumulated_depth_hist[13] + dh3.b;
    cumulated_depth_hist[15] = cumulated_depth_hist[14] + dh3.a;

    cumulated_depth_hist[16] = cumulated_depth_hist[15] + dh4.r;
    cumulated_depth_hist[17] = cumulated_depth_hist[16] + dh4.g;
    cumulated_depth_hist[18] = cumulated_depth_hist[17] + dh4.b;
    cumulated_depth_hist[19] = cumulated_depth_hist[18] + dh4.a;

    cumulated_depth_hist[20] = cumulated_depth_hist[19] + dh5.r;
    cumulated_depth_hist[21] = cumulated_depth_hist[20] + dh5.g;
    cumulated_depth_hist[22] = cumulated_depth_hist[21] + dh5.b;
    cumulated_depth_hist[23] = cumulated_depth_hist[22] + dh5.a;

    cumulated_depth_hist[24] = cumulated_depth_hist[23] + dh6.r;
    cumulated_depth_hist[25] = cumulated_depth_hist[24] + dh6.g;
    cumulated_depth_hist[26] = cumulated_depth_hist[25] + dh6.b;
    cumulated_depth_hist[27] = cumulated_depth_hist[26] + dh6.a;

    cumulated_depth_hist[28] = cumulated_depth_hist[27] + dh7.r;
    cumulated_depth_hist[29] = cumulated_depth_hist[28] + dh7.g;
    cumulated_depth_hist[30] = cumulated_depth_hist[29] + dh7.b;
    cumulated_depth_hist[31] = cumulated_depth_hist[30] + dh7.a;

    //rescale histogram such that the last value coincides with the number of desired color bins
    float hscale = float(numbins) / cumulated_depth_hist[31];

    for(int i = 0; i < NUM_HIST_BINS; ++i)
        cumulated_depth_hist[i] *= hscale;
}




// old shader code
/*
float getViewSpaceDepth(float d)
{
    return - proj.x / (2.0 * d - 1.0 + proj.y);
}

uint fcount = uint(texelFetch(frag_count, frag_coord).r + 0.5);
    if(fcount == 0)
    {
        color = vec4(0.0f, 0.0f, 0.0f, 1.0f);
        return;
    }

    FragData flist[MAX_FRAGMENTS];
    
    int index = texelFetch(im_offsetbuffer, ivec2(gl_FragCoord.xy), 0).r;
    float totalrev = 1.0f;
    fcount = 0;
    while(index > -1 && fcount < MAX_FRAGMENTS)
    {
        Fragment f = sb_fraglist.fragments[index];
        flist[fcount++] = FragData(f.color, f.depth);
        totalrev *= (1.0f - f.color.a);
        index = f.next;        
    }    

    vec3 acc_c = vec3(0.0f, 0.0f, 0.0f);
    float acc_n = 0.0f;
    float acc_rev = 1.0f;

    vec2 bzlim = vec2(ub_bucketlist.buckets[0].x, ub_bucketlist.buckets[numbuckets - 1].y);

    //the first fragments should start at the center of the first bucket
    float basedepth = texelFetch(tdepthmap, ivec2(gl_FragCoord.xy), 0).r + abs(ub_bucketlist.buckets[0].y - ub_bucketlist.buckets[0].x) / 2.0f;
    
    for(uint b = 0; b < numbuckets; ++b)
    {
        vec2 cb = vec2(ub_bucketlist.buckets[b].x, ub_bucketlist.buckets[b].y);

        vec3 bck_c = vec3(0.0f, 0.0f, 0.0f);
        float bck_n = 0.0f;
        float bck_rev = 1.0f;
        
        for(uint f = 0; f < fcount; ++f)
        {
            float z = flist[f].depth - basedepth;
            float a = flist[f].color.a;
            float b = bucketweight(z, cb.x, cb.y);
            float w = depthFunc(z - dwoffset, a, bzlim.x, bzlim.y);
            float tw = a * w * b;
            bck_c += flist[f].color.rgb * tw;
            bck_rev *= (1.0f - a * b);
            bck_n += tw;
        }

        acc_c += (bck_c / max(bck_n, 1e-10)) * acc_rev * (1.0f - bck_rev);
        acc_n += acc_rev * (1.0f - bck_rev);
        acc_rev *= bck_rev;

        if(acc_rev <= BACK_REV_THRESHOLD)
            break;      
    }

    color = vec4((acc_c / max(acc_n, 1e-10)) * (1.0f - totalrev), totalrev);
*/
#version 440 core

in float viewSpaceDepth;
layout(location = 0) out vec2 min_max_depth;
layout(location = 1) out float fragment_count;
void main()
{
    min_max_depth.r = viewSpaceDepth;
    min_max_depth.g = -viewSpaceDepth;
    fragment_count = 1.0;
}
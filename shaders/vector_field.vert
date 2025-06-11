// Shader for instanced vector field rendering.
// Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
#version 140
#extension GL_ARB_explicit_attrib_location : enable

// Vertex attributes
layout (location = 0) in vec3 v_position; // bind v_position to attribute 0
layout (location = 1) in vec3 v_normal;   // bind v_normal   to attribute 1

layout (location = 2) in vec3 arrowPos;
layout (location = 3) in vec3 arrowVec;
layout (location = 4) in vec4 arrowColor;

uniform float u_arrowAlignment;          // offset of tail from "arrowPos" in units of arrowVec (for aligning the arrow center or tip at "arrowPos")
uniform float u_arrowRelativeScreenSize; // desired length in pixels of the longest arrow relative to the viewport size
uniform float u_targetDepth;

uniform vec3 u_lightEyePos;
uniform vec3 u_diffuseIntensity;
uniform vec3 u_ambientIntensity;

// Transformation matrices
uniform mat4 u_modelViewMatrix;
uniform mat4 u_projectionMatrix;
uniform mat3 u_normalMatrix;

out vec4 v2f_color;

void main() {
    mat3 rotmag; // rotation matrix scaled by the vector magnitude

    float len = length(arrowVec);
    rotmag[0] = arrowVec;
    // cross([0, 0, 1], dir) or cross([0, 1, 0], dir depending which entry of dir is smaller
    rotmag[1] = normalize((abs(arrowVec[2]) < abs(arrowVec[1])) ? vec3(-arrowVec[1], arrowVec[0], 0) : vec3(arrowVec[2], 0, -arrowVec[0]));
    rotmag[2] = cross(rotmag[1], arrowVec);
    rotmag[1] *= len;

    // Determine the NDC length of a unit "reference arrow" vector lying parallel to the eye space's x axis and emanating from the *view's target*
    // We use this to normalize vectors to have a user-specified pixel size.
    // First determine how this vector is stretched by the modelView matrix...
    float s = 1.0 / length(vec3(u_normalMatrix[0].x, u_normalMatrix[1].x, u_normalMatrix[2].x)); // reciprocal of norm of inverse model view matrix column 0
    // Then determine the length of the final reference vector segment [s, 0, objectOriginDepth, 1] - [0, 0, objectOriginDepth, 1] in NDC
    float referenceArrowLen = s * length(u_projectionMatrix[0].xyz) / u_targetDepth;
    float scale = 2.0 * u_arrowRelativeScreenSize / referenceArrowLen;

    // Do Gouraud shading, diffuse lighting in eye space
    vec3 pos     = arrowPos + scale * (u_arrowAlignment * arrowVec + rotmag * v_position);
    vec4 eye_pos = u_modelViewMatrix * vec4(pos, 1.0);
    vec3 l       = u_lightEyePos - eye_pos.xyz;          // unnormalized
    vec3 n       = u_normalMatrix * (rotmag * v_normal); // unnormalized
    v2f_color.rgb = (u_ambientIntensity + u_diffuseIntensity * (dot(n, l) / (length(l) * length(n)))) * arrowColor.rgb;
    v2f_color.a   = arrowColor.a;
    gl_Position = u_projectionMatrix * (eye_pos);
}

////////////////////////////////////////////////////////////////////////////////
/*! @file
//  Thick line rendering in OpenGL 3.3.
//  GL_LINE_SMOOTH + glLineWidth is no longer supported on many drivers.
//  The following is an improved version of the approach
//  suggested in https://stackoverflow.com/questions/60440682/drawing-a-line-in-modern-opengl
//
//  Author:  Julian Panetta (jpanetta), jpanetta@ucdavis.edu
//  Company:  University of California, Davis
*///////////////////////////////////////////////////////////////////////////////
#version 150
uniform mat4  u_mvp;
uniform vec2  u_resolution;
uniform float u_thickness;
uniform bool  u_lineStripMode; // GL_LINE_STRIP vs GL_LINES behavior
uniform bool  u_contiguous;

uniform samplerBuffer u_vcoord_tex;

float miter_limit = cos(7 * 3.14159265358979f / 16);

// out vec4 v2f_color; // for debugging

void main() {
    int li = gl_VertexID / 6;
    int local_vi = gl_VertexID % 6;

    // 0-2 4
    // |/ /|
    // 1 5-3
    int endpoint = int(local_vi < 2 || local_vi > 4); // 0 for left, 1 for right
    float normalOffset = local_vi % 2 == 0 ? u_thickness : -u_thickness;

    // Get stencil points and edge unit normals in half-pixel space for
    // the current corner.
    vec4 p[3];
    vec2 n[2];

    for (int i = 0; i < 3; ++i) {
        if (u_lineStripMode) p[i] = u_mvp * texelFetch(u_vcoord_tex,     li + endpoint + i);
        else                 p[i] = u_mvp * texelFetch(u_vcoord_tex, 2 * li + endpoint + i - 1); // undefined data from out-of-bounds sampling will be ignored.
        p[i].xy *= u_resolution / p[i].w;
    }

    vec2 t;
    for (int i = 0; i < 2; ++i) {
        t = normalize(p[i + 1].xy - p[i].xy);
        n[i] = vec2(-t.y, t.x);
    }

    float convexSide = dot(n[0], t) < 0 ? 1.0f : -1.0f;
    // v2f_color = (normalOffset * convexSide < 0) ? vec4(1.0, 0.0, 0.0, 1.0) : vec4(0.0, 1.0, 0.0, 1.0);
    // TODO: nicer limited-miter geometry (4 tris/line, conditional on convexSide?)

    vec2 normal = n[1 - endpoint]; // normal of current line
    if (u_contiguous && u_lineStripMode) {
        vec2 endpointNormal  = normalize(n[0] + n[1]);
        float cos_theta = dot(endpointNormal, normal);
        if (cos_theta > miter_limit) {
            normal = endpointNormal;
            normalOffset /= cos_theta;
        }
    }

    gl_Position = vec4((p[1].xy + normal * normalOffset) * (p[1].w / u_resolution), p[1].zw);
}

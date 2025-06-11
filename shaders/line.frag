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

uniform vec4 u_color;
out vec4 fragColor;

// in vec4 v2f_color; // for debugging

void main() {
    fragColor = u_color;
}

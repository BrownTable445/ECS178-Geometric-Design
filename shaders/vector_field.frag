// Shader for instanced vector field rendering.
// Author:  Julian Panetta (jpanetta), julian.panetta@gmail.com
#version 140

in vec4 v2f_color;
out vec4 fragColor;

void main() {
    // Gouraud shading
    fragColor = v2f_color;
}

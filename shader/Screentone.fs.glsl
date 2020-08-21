#version 440 core

uniform int triangleType;

uniform sampler2D tex;
in vec2 uv;

void main(void)
{
    if (triangleType <= 1)
        gl_FragColor = vec4(0.0, 0.0, 0.0, 1.0); 
    else if (triangleType <= 2)
        gl_FragColor = vec4(1.0, 1.0, 1.0, 1.0);
    else {
        gl_FragColor = texture(tex, uv).rgba;
    }
}
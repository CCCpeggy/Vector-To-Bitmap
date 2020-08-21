#version 440 core

uniform int inImgWidth;
uniform int inImgHeight;

uniform int texWidth;
uniform int texHeight;

uniform int outImgWidth;
uniform int outImgHeight;

uniform int windowWidth;
uniform int windowHeight;

layout(location = 0) in vec2 position;
out vec2 uv;

void main(void)
{
	vec2 pos = vec2(position.x / inImgWidth * 2 * outImgWidth / windowWidth - 1, position.y / inImgHeight * 2 * outImgHeight / windowHeight - 1);
	uv = vec2(position.x / inImgWidth * windowWidth / texWidth, position.y / inImgHeight * windowHeight / texHeight);
	gl_Position = vec4(pos, 0.0, 1.0);
}
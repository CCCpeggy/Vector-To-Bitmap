#version 430 core

out VertexData
{
	vec2 texcoord;
} vertexData;
const vec2[] points = vec2[4](vec2(-1, -1), vec2(1, -1), vec2(-1, 1), vec2(1, 1));
const vec2[] uv = vec2[4](vec2(0, 1), vec2(1, 1), vec2(0, 0), vec2(1, 0));

void main()
{
	gl_Position = vec4(points[gl_VertexID], 0.0, 1.0);
	vertexData.texcoord = uv[gl_VertexID];
}

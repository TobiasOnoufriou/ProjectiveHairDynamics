#version 450 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aNormal;

uniform mat4 u_ViewProjection;

out vec3 v_position;
out vec3 v_normal;

void main()
{
    v_position = aPos;
    v_normal = aNormal;
    gl_Position = u_ViewProjection * vec4(aPos, 1.0);
}

#version 450 core
in vec3 v_position;
in vec3 v_normal;

uniform vec3 u_CameraPosition;
uniform vec3 u_Color;

out vec4 pixel;

// Two-sided Lambert from a fixed key light plus a headlight, and a little
// Blinn-Phong highlight, so the deforming surface reads clearly.
void main()
{
    vec3 n = normalize(v_normal);
    vec3 view = normalize(u_CameraPosition - v_position);
    if (dot(n, view) < 0.0) n = -n;  // Show back faces (e.g. when the camera is inside).

    const vec3 base = u_Color;
    const vec3 key = normalize(vec3(0.4, 1.0, 0.3));
    float diffuse = 0.65 * max(dot(n, key), 0.0) + 0.35 * max(dot(n, view), 0.0);
    float specular = pow(max(dot(n, normalize(key + view)), 0.0), 40.0);
    pixel = vec4(base * (0.15 + diffuse) + vec3(0.25) * specular, 1.0);
}

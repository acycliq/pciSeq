const fShader = `
#define PI 3.14159265359

float lineSegment(vec2 p, vec2 a, vec2 b) {
    float thickness = 1.0/100.0;
    vec2 pa = p - a, ba = b - a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
    return step(0.05, length(pa - ba*h));
}

mat2 rot(float a){
    return mat2( cos(a), sin(a), -sin(a), cos(a));
}

void main() 
{  
    vec2 uv = (gl_PointCoord - vec2(0.5)) * rot(PI); // rotate
    uv = uv + vec2(0.5);
    float line_1 = 1.0 - lineSegment(uv, vec2(0.05), vec2(0.5, 0.95));
    float line_2 = 1.0 - lineSegment(uv, vec2(0.5, 0.95), vec2(0.95, 0.05));
    float line_3 = 1.0 - lineSegment(uv, vec2(0.95, 0.05), vec2(0.05, 0.05));
    float shaper = line_1 + line_2 + line_3;

    gl_FragColor = vec4(vec3(shaper), 1.0);
}
`
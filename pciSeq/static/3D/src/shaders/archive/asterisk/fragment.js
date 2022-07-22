const fShader = `

float lineSegment(vec2 p, vec2 a, vec2 b) {
    float thickness = 1.0/100.0;
    vec2 pa = p - a, ba = b - a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
    return step(0.05, length(pa - ba*h));
}

void main() 
{
    
    float diag_1 = 1.0 - lineSegment(gl_PointCoord, vec2(0.30, 0.30), vec2(0.70, 0.70));
    float diag_2 = 1.0 - lineSegment(gl_PointCoord, vec2(0.30, 0.70), vec2(0.70, 0.30));
    float cross = diag_1 + diag_2;
    
    float vertical = 1.0 - lineSegment(gl_PointCoord, vec2(0.5, 0.05), vec2(0.5, 0.95));
    float horizontal = 1.0 - lineSegment(gl_PointCoord, vec2(0.10, 0.5), vec2(0.90, 0.5));
    float plus = vertical + horizontal;
    
    float shaper = plus + cross;

    gl_FragColor = vec4(vec3(shaper), 1.0);
}
`
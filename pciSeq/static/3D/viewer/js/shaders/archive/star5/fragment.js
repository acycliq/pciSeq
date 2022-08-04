const fShader = `

float lineSegment(vec2 p, vec2 a, vec2 b) {
    float thickness = 1.0/100.0;
    vec2 pa = p - a, ba = b - a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
    return step(0.05, length(pa - ba*h));
}


void main() 
{

    float side = 1.0;  // fragment has side length = 1.0
    float r = side/2.0;
    vec2 p = vec2(r);  // center of the fragment
    float eps = 0.1;   // will be adjusting the tips of the star to avoid clipping
    
    // Points start at the bottom left tip of the star and move counter clockwise
    vec2 A = vec2(p.x + 0.618027443 * r, p.y + (1.0 - eps) * r); //
    vec2 B = vec2(p.x + 0.5 * r, p.y + 0.27637816 * r);
    vec2 C = vec2(p.x + (1.0 - eps) * r, p.y - 0.236080209 * r);
    vec2 D = vec2(p.x + 0.309026865 * r, p.y - 0.341634306 * r);
    vec2 E = vec2(p.x + 0.0 * r, p.y - (1.0 - eps) * r);
    vec2 F = vec2(p.x - 0.309026865 * r, p.y - 0.341634306 * r);
    vec2 G = vec2(p.x - (1.0 - eps) * r, p.y - 0.236080209 * r);
    vec2 H = vec2(p.x - 0.5 * r, p.y + 0.27637816 * r);
    vec2 I = vec2(p.x - 0.618027443 * r, p.y + (1.0 - eps) * r);
    vec2 J = vec2(p.x, p.y + 0.658351875 * r);
    vec2 K = vec2(p.x, p.y + 0.658351875 * r);
    
    // Draw now the star
    float line_1  = 1.0 - lineSegment(gl_PointCoord, A, B);
    float line_2  = 1.0 - lineSegment(gl_PointCoord, B, C);
    float line_3  = 1.0 - lineSegment(gl_PointCoord, C, D);
    float line_4  = 1.0 - lineSegment(gl_PointCoord, D, E);
    float line_5  = 1.0 - lineSegment(gl_PointCoord, E, F);
    float line_6  = 1.0 - lineSegment(gl_PointCoord, F, G);
    float line_7  = 1.0 - lineSegment(gl_PointCoord, G, H);
    float line_8  = 1.0 - lineSegment(gl_PointCoord, H, I);
    float line_9  = 1.0 - lineSegment(gl_PointCoord, I, J);
    float line_10 = 1.0 - lineSegment(gl_PointCoord, J, K);
    float line_11 = 1.0 - lineSegment(gl_PointCoord, K, A);
    
    float shaper = line_1 + line_2 + line_3 + line_4 + line_5 + line_6 + line_7 + line_8 + line_9 + line_10 + line_11;

    gl_FragColor = vec4(vec3(shaper), 1.0);
}
`

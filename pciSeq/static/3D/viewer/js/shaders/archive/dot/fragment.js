const fShader = `

float draw_circle(vec2 coord, float radius) {
    // calculate the distance from the center(which is at [0.5, 0.5])
    float d = distance(gl_PointCoord, vec2(0.5));
    
    // this will return 1.0 for all fragments inside the radius
    return step(d, radius);
}


void main() 
{
    float circle_1 = draw_circle(gl_PointCoord, 0.25);
    float circle_2 = 1.0 - draw_circle(gl_PointCoord, 0.08);
    vec3 color = vec3(circle_1 * circle_2);

    gl_FragColor = vec4(color, 1.0);
}
`


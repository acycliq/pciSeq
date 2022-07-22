const fShader = `
uniform float uSize;
uniform float r;
uniform float g;
uniform float b;
uniform float a;
void main()
{
    float strength = gl_PointCoord.y;
    gl_FragColor = vec4(vec2(gl_PointCoord.y), 1.0, 1.0);
    
}
`

const fShader_asterisk = `
varying vec3 view;
uniform float zThres;

float lineSegment(vec2 p, vec2 a, vec2 b) {
    float thickness = 1.0/100.0;
    vec2 pa = p - a, ba = b - a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
    return step(0.05, length(pa - ba*h));
}

uniform float r;
uniform float g;
uniform float b;
uniform float a;
varying vec3 vColor;

void main() {
  // gl_FragColor = vec4(1.0, 1.0, 0.0, 0.5);
  if (view.z > zThres){
    // points at the far back
    gl_FragColor = vec4(vColor.r, vColor.g, vColor.b, a);
  }
  else {
    float diag_1 = 1.0 - lineSegment(gl_PointCoord, vec2(0.30, 0.30), vec2(0.70, 0.70));
    float diag_2 = 1.0 - lineSegment(gl_PointCoord, vec2(0.30, 0.70), vec2(0.70, 0.30));
    float cross = diag_1 + diag_2;
    
    float vertical = 1.0 - lineSegment(gl_PointCoord, vec2(0.5, 0.05), vec2(0.5, 0.95));
    float horizontal = 1.0 - lineSegment(gl_PointCoord, vec2(0.10, 0.5), vec2(0.90, 0.5));
    float plus = vertical + horizontal;
    
    float shaper = plus + cross;

    gl_FragColor = vec4(vec3(shaper) * vColor, a);
    if (gl_FragColor.r + gl_FragColor.g + gl_FragColor.b == 0.0) discard;
  }
}
`;
// ---------------------------------------------------------------------------------------------------------------------

const fShader_circle = `
varying vec3 view;
uniform float zThres;

float draw_circle(vec2 coord, float radius) {
    // calculate the distance from the center(which is at [0.5, 0.5])
    float d = distance(gl_PointCoord, vec2(0.5));
    
    // this will return 1.0 for all fragments inside the radius
    return step(d, radius);
}

uniform float r;
uniform float g;
uniform float b;
uniform float a;
varying vec3 vColor;

void main() {
  // gl_FragColor = vec4(1.0, 1.0, 0.0, 0.5);
  if (view.z > zThres){
    // points at the far back
    gl_FragColor = vec4(vColor.r, vColor.g, vColor.b, a);
  }
  else {
    float circle_1 = draw_circle(gl_PointCoord, 0.5);
    float circle_2 = 1.0 - draw_circle(gl_PointCoord, 0.4);
    vec3 color = vec3(circle_1 * circle_2);

    gl_FragColor = vec4(color  * vColor, a);
    if (gl_FragColor.r + gl_FragColor.g + gl_FragColor.b == 0.0) discard;
  }  
}
`;
// ---------------------------------------------------------------------------------------------------------------------

const fShader_cross = `
varying vec3 view;
uniform float zThres;

float lineSegment(vec2 p, vec2 a, vec2 b) {
    float thickness = 1.0/100.0;
    vec2 pa = p - a, ba = b - a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
    return step(0.05, length(pa - ba*h));
}

uniform float r;
uniform float g;
uniform float b;
uniform float a;
varying vec3 vColor;

void main() {
  // gl_FragColor = vec4(1.0, 1.0, 0.0, 0.5);
  if (view.z > zThres){
    // points at the far back
    gl_FragColor = vec4(vColor.r, vColor.g, vColor.b, a);
  }
  else {
    
    float diag_1 = 1.0 - lineSegment(gl_PointCoord, vec2(0.05, 0.05), vec2(0.95, 0.95));
    float diag_2 = 1.0 - lineSegment(gl_PointCoord, vec2(0.05, 0.95), vec2(0.95, 0.05));
    
    float shaper = diag_1 + diag_2;

    gl_FragColor = vec4(vec3(shaper) * vColor, a);
    if (gl_FragColor.r + gl_FragColor.g + gl_FragColor.b == 0.0) discard;
  }  
}
`;
// ---------------------------------------------------------------------------------------------------------------------

const fShader_diamond = `
varying vec3 view;
uniform float zThres;

float lineSegment(vec2 p, vec2 a, vec2 b) {
    float thickness = 1.0/100.0;
    vec2 pa = p - a, ba = b - a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
    return step(0.05, length(pa - ba*h));
}

uniform float r;
uniform float g;
uniform float b;
uniform float a;
varying vec3 vColor;

void main() {
  // gl_FragColor = vec4(1.0, 1.0, 0.0, 0.5);
  if (view.z > zThres){
    // points at the far back
    gl_FragColor = vec4(vColor.r, vColor.g, vColor.b, a);
  }
  else {
    float side_1 = 1.0 - lineSegment(gl_PointCoord, vec2(0.5, 0.05), vec2(0.9, 0.5));
    float side_2 = 1.0 - lineSegment(gl_PointCoord, vec2(0.9, 0.5), vec2(0.5, 0.95));
    float side_3 = 1.0 - lineSegment(gl_PointCoord, vec2(0.5, 0.95), vec2(0.1, 0.5));
    float side_4 = 1.0 - lineSegment(gl_PointCoord, vec2(0.1, 0.5), vec2(0.5, 0.05));
    
    float shaper = side_1 + side_2 + side_3 + side_4;

    gl_FragColor = vec4(vec3(shaper) * vColor, a);
    if (gl_FragColor.r + gl_FragColor.g + gl_FragColor.b == 0.0) discard;
  }
}
`;
// ---------------------------------------------------------------------------------------------------------------------

const fShader_point = `
varying vec3 view;
uniform float zThres;

float draw_circle(vec2 coord, float radius) {
    // calculate the distance from the center(which is at [0.5, 0.5])
    float d = distance(gl_PointCoord, vec2(0.5));
    
    // this will return 1.0 for all fragments inside the radius
    return step(d, radius);
}

uniform float r;
uniform float g;
uniform float b;
uniform float a;
varying vec3 vColor;

void main() {
  // gl_FragColor = vec4(1.0, 1.0, 0.0, 0.5);
  if (view.z > zThres){
    // points at the far back
    gl_FragColor = vec4(vColor.r, vColor.g, vColor.b, a);
  }
  else {
    float circle_1 = draw_circle(gl_PointCoord, 0.25);
    float circle_2 = 1.0 - draw_circle(gl_PointCoord, 0.08);
    vec3 color = vec3(circle_1 * circle_2);

    gl_FragColor = vec4(color * vColor, a);
    if (gl_FragColor.r + gl_FragColor.g + gl_FragColor.b == 0.0) discard;
  }  
}
`;
// ---------------------------------------------------------------------------------------------------------------------

const fShader_plus = `
varying vec3 view;
uniform float zThres;

float lineSegment(vec2 p, vec2 a, vec2 b) {
    float thickness = 1.0/100.0;
    vec2 pa = p - a, ba = b - a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
    return step(0.05, length(pa - ba*h));
}

uniform float r;
uniform float g;
uniform float b;
uniform float a;
varying vec3 vColor;

void main() {
  // gl_FragColor = vec4(1.0, 1.0, 0.0, 0.5);
  if (view.z > zThres){
    // points at the far back
    gl_FragColor = vec4(vColor.r, vColor.g, vColor.b, a);
  }
  else {
    float vertical = 1.0 - lineSegment(gl_PointCoord, vec2(0.5, 0.05), vec2(0.5, 0.95));
    float horizontal = 1.0 - lineSegment(gl_PointCoord, vec2(0.05, 0.5), vec2(0.95, 0.5));
    float shaper = vertical + horizontal;

    gl_FragColor = vec4(vec3(shaper) * vColor, a);
    if (gl_FragColor.r + gl_FragColor.g + gl_FragColor.b == 0.0) discard;
  }  
}
`;
// ---------------------------------------------------------------------------------------------------------------------

const fShader_square = `
varying vec3 view;
uniform float zThres;

float lineSegment(vec2 p, vec2 a, vec2 b) {
    float thickness = 1.0/100.0;
    vec2 pa = p - a, ba = b - a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
    return step(0.05, length(pa - ba*h));
}

uniform float r;
uniform float g;
uniform float b;
uniform float a;
varying vec3 vColor;

void main() {
  // gl_FragColor = vec4(1.0, 1.0, 0.0, 0.5);
  if (view.z > zThres){
    // points at the far back
    gl_FragColor = vec4(vColor.r, vColor.g, vColor.b, a);
  }
  else {
    float top = 1.0 - lineSegment(gl_PointCoord, vec2(0.05, 0.05), vec2(0.95, 0.05));
    float right = 1.0 - lineSegment(gl_PointCoord, vec2(0.95, 0.05), vec2(0.95, 0.95));
    float bottom = 1.0 - lineSegment(gl_PointCoord, vec2(0.95, 0.95), vec2(0.05, 0.95));
    float left = 1.0 - lineSegment(gl_PointCoord, vec2(0.05, 0.95), vec2(0.05, 0.05));
    float shaper = top + right + bottom + left;

    gl_FragColor = vec4(vec3(shaper) * vColor, a);
    if (gl_FragColor.r + gl_FragColor.g + gl_FragColor.b == 0.0) discard;
  }
}
`;
// ---------------------------------------------------------------------------------------------------------------------

const fShader_star5 = `
varying vec3 view;
uniform float zThres;

float lineSegment(vec2 p, vec2 a, vec2 b) {
    float thickness = 1.0/100.0;
    vec2 pa = p - a, ba = b - a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
    return step(0.05, length(pa - ba*h));
}

uniform float r;
uniform float g;
uniform float b;
uniform float a;
varying vec3 vColor;

void main() {
  // gl_FragColor = vec4(1.0, 1.0, 0.0, 0.5);
  if (view.z > zThres){
    // points at the far back
    gl_FragColor = vec4(vColor.r, vColor.g, vColor.b, a);
  }
  else {
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

    gl_FragColor = vec4(vec3(shaper) * vColor, a);
    if (gl_FragColor.r + gl_FragColor.g + gl_FragColor.b == 0.0) discard;
  }
}
`;
// ---------------------------------------------------------------------------------------------------------------------

const fShader_star6 = `
varying vec3 view;
uniform float zThres;

// IT'S OK BUT NOT PERFECT. THE TIPS ARE GETTING CLIPPED AND NEEDS SOME REVISION TO FIX THIS 

float lineSegment(vec2 p, vec2 a, vec2 b) {
    float thickness = 1.0/100.0;
    vec2 pa = p - a, ba = b - a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
    return step(0.05, length(pa - ba*h));
}

uniform float r;
uniform float g;
uniform float b;
uniform float a;
varying vec3 vColor;

void main() {
  // gl_FragColor = vec4(1.0, 1.0, 0.0, 0.5);
  if (view.z > zThres){
    // points at the far back
    gl_FragColor = vec4(vColor.r, vColor.g, vColor.b, a);
  }
  else {
    float side = 1.0;  // fragment has side length = 1.0
    float r = side/2.0;
    vec2 p = vec2(r);  // center of the fragment
    float eps = 0.05;   // will be adjusting the tips of the star to avoid clipping
    
    // Points start at the bottom left tip of the star and move clockwise
    vec2 A = vec2(p.x + 0.50 * r, p.y + 0.87 * r);
    vec2 B = vec2(p.x, p.y + 0.50 * r);
    vec2 C = vec2(p.x - 0.50 * r, p.y + 0.87 * r);
    vec2 D = vec2(p.x - 0.43 * r, p.y + 0.25 * r);
    vec2 E = vec2(p.x - r, p.y);
    vec2 F = vec2(p.x - 0.43 * r, p.y - 0.25 * r);
    vec2 G = vec2(p.x - 0.50 * r, p.y - 0.87 * r);
    vec2 H = vec2(p.x, p.y - 0.50 * r);
    vec2 I = vec2(p.x + 0.50 * r, p.y - 0.87 * r);
    vec2 J = vec2(p.x + 0.43 * r, p.y - 0.25 * r);
    vec2 K = vec2(p.x + r, p.y);
    vec2 L = vec2(p.x + 0.43 * r, p.y + 0.25 * r);
    
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
    float line_11 = 1.0 - lineSegment(gl_PointCoord, K, L);
    float line_12 = 1.0 - lineSegment(gl_PointCoord, L, A);
    
    
    float shaper = line_1 + line_2 + line_3 + line_4 + line_5 + line_6 + line_7 + line_8 + line_9 + line_10 + line_11 + line_12;

    gl_FragColor = vec4(vec3(shaper) * vColor, a);
    if (gl_FragColor.r + gl_FragColor.g + gl_FragColor.b == 0.0) discard;
  }  
}
`;
// ---------------------------------------------------------------------------------------------------------------------

const fShader_triangleDown = `
varying vec3 view;
uniform float zThres;

float lineSegment(vec2 p, vec2 a, vec2 b) {
    float thickness = 1.0/100.0;
    vec2 pa = p - a, ba = b - a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
    return step(0.05, length(pa - ba*h));
}

uniform float r;
uniform float g;
uniform float b;
uniform float a;
varying vec3 vColor;

void main() {
  // gl_FragColor = vec4(1.0, 1.0, 0.0, 0.5);
  if (view.z > zThres){
    // points at the far back
    gl_FragColor = vec4(vColor.r, vColor.g, vColor.b, a);
  }
  else {
    float line_1 = 1.0 - lineSegment(gl_PointCoord, vec2(0.05), vec2(0.5, 0.95));
    float line_2 = 1.0 - lineSegment(gl_PointCoord, vec2(0.5, 0.95), vec2(0.95, 0.05));
    float line_3 = 1.0 - lineSegment(gl_PointCoord, vec2(0.95, 0.05), vec2(0.05, 0.05));
    float shaper = line_1 + line_2 + line_3;

    gl_FragColor = vec4(vec3(shaper) * vColor, a);
    if (gl_FragColor.r + gl_FragColor.g + gl_FragColor.b == 0.0) discard;
  }
}
`;
// ---------------------------------------------------------------------------------------------------------------------

const fShader_triangleLeft = `
#define PI 3.14159265359
varying vec3 view;
uniform float zThres;

float lineSegment(vec2 p, vec2 a, vec2 b) {
    float thickness = 1.0/100.0;
    vec2 pa = p - a, ba = b - a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
    return step(0.05, length(pa - ba*h));
}

mat2 rot(float a){
    return mat2( cos(a), sin(a), -sin(a), cos(a));
}

uniform float r;
uniform float g;
uniform float b;
uniform float a;
varying vec3 vColor;

void main() {
  // gl_FragColor = vec4(1.0, 1.0, 0.0, 0.5);
  if (view.z > zThres){
    // points at the far back
    gl_FragColor = vec4(vColor.r, vColor.g, vColor.b, a);
  }
  else {
    vec2 uv = (gl_PointCoord - vec2(0.5)) * rot(PI * 0.5); // rotate
    uv = uv + vec2(0.5);
    float line_1 = 1.0 - lineSegment(uv, vec2(0.05), vec2(0.5, 0.95));
    float line_2 = 1.0 - lineSegment(uv, vec2(0.5, 0.95), vec2(0.95, 0.05));
    float line_3 = 1.0 - lineSegment(uv, vec2(0.95, 0.05), vec2(0.05, 0.05));
    float shaper = line_1 + line_2 + line_3;

    gl_FragColor = vec4(vec3(shaper) * vColor, a);
    if (gl_FragColor.r + gl_FragColor.g + gl_FragColor.b == 0.0) discard;
  }  
}
`;
// ---------------------------------------------------------------------------------------------------------------------

const fShader_triangleRight = `
#define PI 3.14159265359
varying vec3 view;
uniform float zThres;

float lineSegment(vec2 p, vec2 a, vec2 b) {
    float thickness = 1.0/100.0;
    vec2 pa = p - a, ba = b - a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
    return step(0.05, length(pa - ba*h));
}

mat2 rot(float a){
    return mat2( cos(a), sin(a), -sin(a), cos(a));
}

uniform float r;
uniform float g;
uniform float b;
uniform float a;
varying vec3 vColor;

void main() {
  // gl_FragColor = vec4(1.0, 1.0, 0.0, 0.5);
  if (view.z > zThres){
    // points at the far back
    gl_FragColor = vec4(vColor.r, vColor.g, vColor.b, a);
  }
  else {
    vec2 uv = (gl_PointCoord - vec2(0.5)) * rot(1.5 * PI); // rotate
    uv = uv + vec2(0.5);
    float line_1 = 1.0 - lineSegment(uv, vec2(0.05), vec2(0.5, 0.95));
    float line_2 = 1.0 - lineSegment(uv, vec2(0.5, 0.95), vec2(0.95, 0.05));
    float line_3 = 1.0 - lineSegment(uv, vec2(0.95, 0.05), vec2(0.05, 0.05));
    float shaper = line_1 + line_2 + line_3;

    gl_FragColor = vec4(vec3(shaper) * vColor, a);
    if (gl_FragColor.r + gl_FragColor.g + gl_FragColor.b == 0.0) discard;
  }  
}
`;
// ---------------------------------------------------------------------------------------------------------------------

const fShader_triangleUp = `
#define PI 3.14159265359
varying vec3 view;
uniform float zThres;

float lineSegment(vec2 p, vec2 a, vec2 b) {
    float thickness = 1.0/100.0;
    vec2 pa = p - a, ba = b - a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
    return step(0.05, length(pa - ba*h));
}

mat2 rot(float a){
    return mat2( cos(a), sin(a), -sin(a), cos(a));
}

uniform float r;
uniform float g;
uniform float b;
uniform float a;
varying vec3 vColor;

void main() {
  // gl_FragColor = vec4(1.0, 1.0, 0.0, 0.5);
  if (view.z > zThres){
    gl_FragColor = vec4(vColor.r, vColor.g, vColor.b, a);
  }
  else {
    vec2 uv = (gl_PointCoord - vec2(0.5)) * rot(PI); // rotate
    uv = uv + vec2(0.5);
    float line_1 = 1.0 - lineSegment(uv, vec2(0.05), vec2(0.5, 0.95));
    float line_2 = 1.0 - lineSegment(uv, vec2(0.5, 0.95), vec2(0.95, 0.05));
    float line_3 = 1.0 - lineSegment(uv, vec2(0.95, 0.05), vec2(0.05, 0.05));
    float shaper = line_1 + line_2 + line_3;

    gl_FragColor = vec4(vec3(shaper) * vColor, a);
    if (gl_FragColor.r + gl_FragColor.g + gl_FragColor.b == 0.0) discard;
  }
}
`;

const fShader_solidSquare = `
uniform vec2 u_resolution;
uniform vec2 u_mouse;
uniform float u_time;
uniform float zThres;
uniform float r;
uniform float g;
uniform float b;
uniform float a;
varying vec3 view;


void main(){
      if (view.z > zThres){
      // points at the far back
        gl_FragColor = vec4(r, g, b, a);
      }
      else {
      // points close to the camera
        vec2 st = gl_FragCoord.xy/u_resolution;
        float pct = 0.0;
    
        // a. The DISTANCE from the pixel to the center
        pct = distance(st, vec2(0.5));
    
        // b. The LENGTH of the vector
        //    from the pixel to the center
        // vec2 toCenter = vec2(0.5)-st;
        // pct = length(toCenter);
    
        // c. The SQUARE ROOT of the vector
        //    from the pixel to the center
        // vec2 tC = vec2(0.5)-st;
        // pct = sqrt(tC.x*tC.x+tC.y*tC.y);
    
        vec3 color = vec3(pct);
        gl_FragColor = vec4(1.0, 0.0, 0.0, a);
      }
}
`;
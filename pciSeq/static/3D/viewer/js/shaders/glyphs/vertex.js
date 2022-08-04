const vShader_glyphs = `
uniform float glyphSize;
uniform float dotSize;
uniform float zThres;
uniform bool attenuate;
varying vec3 view;
attribute float scale;
attribute float size;
attribute vec3 mycolor;
varying vec3 vColor;
void main()
{
    vColor = mycolor;
    
    // Position
    vec4 modelPosition = modelMatrix * vec4(position, 1.0);
    vec4 viewPosition = viewMatrix * modelPosition;
    vec4 projectedPosition = projectionMatrix * viewPosition;
    gl_Position = projectedPosition;
    
    view = -1.0 * viewPosition.xyz;  // NDC in [-1, 1] (by perspective divide)
    
    // Size
    if (view.z > zThres) {
        // points at the far back
        gl_PointSize = dotSize;
    }
    else {
        // points close to the camera
        gl_PointSize = 100.0 * scale * glyphSize; //1500.0;
        if (attenuate){
            gl_PointSize = size * (1. / - viewPosition.z);
        }
    }
    
    // Controls the attenuation
    //  gl_PointSize *= (1.0 / - viewPosition.z);
    
}
`;


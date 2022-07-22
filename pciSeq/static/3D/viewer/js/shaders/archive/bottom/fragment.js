const fShader = `

void main() {
                vec3 color = vec3(0.0);
                vec3 shape = vec3(0.0);

                // True (ie 1.0) if x between 0.45 and 0.55  
                float x_1 = 1.0 - step(0.55, gl_PointCoord.x);
                float x_2 = step(0.45, gl_PointCoord.x);
                
                // True (ie 1.0) if y between 0.45 and 0.55  
                float y_1 = 1.0 - step(0.55, gl_PointCoord.y);
                float y_2 = step(0.45, gl_PointCoord.y);
                
                // cross
                shape = vec3(x_1 * x_2 + y_1 * y_2);
                color = vec3(1.0, 0.0, 0.5);
                
                gl_FragColor = vec4(shape * color, 1.0);
}
`
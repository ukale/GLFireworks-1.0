//
//  ParticleShader.fsh
//  SteerMe
//
//  Created by Admin Dev on 2/26/15.
//  Copyright (c) 2015 Kinetic Bytes. All rights reserved.
//

// Fragment Shader

static const char* ParticleFS = STRINGIFY
(
 precision highp float;
 
 // Varying
 varying vec3  v_colorOffset;
 varying float v_alpha;
 varying vec2  v_texOffset;
 
 // Uniforms
 uniform vec3 u_color;
 uniform sampler2D u_sampler;
 uniform float u_texScale;
 
 void main(void)
 {
    vec2 pointCoord = gl_PointCoord*u_texScale + v_texOffset;
    vec4 texture = texture2D(u_sampler, pointCoord);

    vec4 color = vec4(clamp(u_color + v_colorOffset, vec3(0.0), vec3(1.0)), v_alpha);
    gl_FragColor = texture * color;
 }
);
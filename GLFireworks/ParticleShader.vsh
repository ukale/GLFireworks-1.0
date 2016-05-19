//
//  ParticleShader.vsh
//  SteerMe
//
//  Created by Admin Dev on 2/26/15.
//  Copyright (c) 2015 Kinetic Bytes. All rights reserved.
//

// Vertex Shader

static const char* ParticleVS = STRINGIFY
(
 
 // Attributes
 attribute float  a_ID;
 attribute float  a_direction;
 attribute float  a_speedOffset;
 attribute float  a_sizeOffset;
 attribute vec3   a_colorOffset;
 attribute vec2   a_texOffset;
 
 // Uniforms
 uniform mat4     u_mvpMatrix;
 uniform float    u_speed;
 uniform float    u_initialSize;
 uniform vec2     u_position;
 uniform float    u_updateFrame;
 uniform float    u_minAlpha;
 
 // Varying
 varying vec3     v_colorOffset;
 varying float    v_alpha;
 varying vec2     v_texOffset;
 
 void main(void)
 {
    float activeFrame = max(0.0, u_updateFrame-a_ID);
    float active = clamp(u_updateFrame-a_ID+1.0, 0.0, 1.0);
    
    const float gravity = 1.0;
    const float drag = 0.984;

    //size decreases as particle flies outward
    float attenuation = pow(drag, activeFrame);
    float s = mix(0.0, u_initialSize, attenuation);

    //direction of outward motion for this particle
    float displacement = (u_speed + a_speedOffset)*(1.0 - attenuation)/(1.0 - drag);    //using summation of geometric series
    float x = cos(a_direction)*displacement;
    float y = sin(a_direction)*displacement + gravity*activeFrame;
    
    //outputs for fragment shader
    v_colorOffset = a_colorOffset;
    v_alpha = max(u_minAlpha, attenuation);  //as size decreases, make it fade away
    v_texOffset = a_texOffset;
    
    //Required OpenGL ES 2.0 outputs
    gl_Position = u_mvpMatrix * vec4(u_position + vec2(x,y), 0.0, 1.0);
    gl_PointSize = active*max(0.0, (s + a_sizeOffset));    //because we are using GL_POINTS for rendering
 }
);
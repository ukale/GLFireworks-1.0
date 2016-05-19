//
//  ParticleShader.h
//  SteerMe
//
//  Created by Admin Dev on 2/26/15.
//  Copyright (c) 2015 Kinetic Bytes. All rights reserved.
//

#ifndef __SteerMe__ParticleShader__
#define __SteerMe__ParticleShader__

#include <OpenGLES/ES2/glext.h>

class ParticleShader
{
public:
    ParticleShader() {}
    virtual ~ParticleShader();
    virtual void Init();
    virtual void Uninit();

protected:
    GLuint BuildProgram(const char* vertexShaderSource, const char* fragmentShaderSource);
    GLuint BuildShader(const char* source, GLenum shaderType);
    
public:
    GLuint program;
    
    //attribute handles
    GLint a_ID;
    GLint a_direction;
    GLint a_speedOffset;
    GLint a_sizeOffset;
    GLint a_colorOffset;
    GLint a_texOffset;
    
    //uniform handles
    GLuint u_mvpMatrix;
    GLint u_speed;
    GLint u_initialSize;
    GLint u_color;
    GLint u_sampler;
    GLint u_position;
    GLint u_updateFrame;
    GLint u_texScale;
    GLint u_minAlpha;           //minimum opacity to maintain for particles as they decay
};

#endif /* defined(__SteerMe__ParticleShader__) */

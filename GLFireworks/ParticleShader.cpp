//
//  ParticleShader.cpp
//  SteerMe
//
//  Created by Admin Dev on 2/26/15.
//  Copyright (c) 2015 Kinetic Bytes. All rights reserved.
//

#include <iostream>
#include "ParticleShader.h"

#define STRINGIFY(A) #A
#include "ParticleShader.vsh"
#include "ParticleShader.fsh"

ParticleShader::~ParticleShader()
{
    Uninit();
}

void ParticleShader::Init()
{
    //Program
    program = BuildProgram(ParticleVS, ParticleFS);
    
    //Attributes
    //get bindings allocated by opengl
    //(or we could use glBindAttributeLocation to explicitly specify bindings if the same VAO needs to be used in multiple shaders)
    a_ID = glGetAttribLocation(program, "a_ID");
    a_direction = glGetAttribLocation(program, "a_direction");
    a_speedOffset = glGetAttribLocation(program, "a_speedOffset");
    a_sizeOffset = glGetAttribLocation(program, "a_sizeOffset");
    a_colorOffset = glGetAttribLocation(program, "a_colorOffset");
    a_texOffset = glGetAttribLocation(program, "a_texOffset");
    
    //Uniforms
    u_mvpMatrix = glGetUniformLocation(program, "u_mvpMatrix");
    u_speed = glGetUniformLocation(program, "u_speed");
    u_initialSize = glGetUniformLocation(program, "u_initialSize");
    u_color = glGetUniformLocation(program, "u_color");
    u_sampler = glGetUniformLocation(program, "u_sampler");
    u_position = glGetUniformLocation(program, "u_position");
    u_updateFrame = glGetUniformLocation(program, "u_updateFrame");
    u_texScale = glGetUniformLocation(program, "u_texScale");
    u_minAlpha = glGetUniformLocation(program, "u_minAlpha");
}

void ParticleShader::Uninit()
{
    glDeleteProgram(program);
}

GLuint ParticleShader::BuildProgram(const char* vertexShaderSource, const char* fragmentShaderSource)
{
    // Build shaders
    GLuint vertexShader = BuildShader(vertexShaderSource, GL_VERTEX_SHADER);
    GLuint fragmentShader = BuildShader(fragmentShaderSource, GL_FRAGMENT_SHADER);
    
    // Create program
    GLuint programHandle = glCreateProgram();
    
    // Attach shaders
    glAttachShader(programHandle, vertexShader);
    glAttachShader(programHandle, fragmentShader);
    
    // Link program
    glLinkProgram(programHandle);
    
    // Check for errors
    GLint linkSuccess;
    glGetProgramiv(programHandle, GL_LINK_STATUS, &linkSuccess);
    if (linkSuccess == GL_FALSE)
    {
        //NSLog(@"GLSL Program Error");
        GLchar messages[1024];
        glGetProgramInfoLog(programHandle, sizeof(messages), 0, &messages[0]);
        std::cout << messages;
        exit(1);
    }
    
    // Delete shaders
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);
    
    return programHandle;
}

GLuint ParticleShader::BuildShader(const char* source, GLenum shaderType)
{
    // Create the shader object
    GLuint shaderHandle = glCreateShader(shaderType);
    
    // Load the shader source
    glShaderSource(shaderHandle, 1, &source, 0);
    
    // Compile the shader
    glCompileShader(shaderHandle);
    
    // Check for errors
    GLint compileSuccess;
    glGetShaderiv(shaderHandle, GL_COMPILE_STATUS, &compileSuccess);
    if (compileSuccess == GL_FALSE)
    {
        //NSLog(@"GLSL Shader Error");
        GLchar messages[1024];
        glGetShaderInfoLog(shaderHandle, sizeof(messages), 0, &messages[0]);
        std::cout << messages;
        exit(1);
    }
    return shaderHandle;
}

//
//  ParticleEmitter.cpp
//  SteerMe
//
//  Created by Admin Dev on 2/26/15.
//  Copyright (c) 2015 Kinetic Bytes. All rights reserved.
//

#include "ParticleEmitter.h"

ParticleEmitter::ParticleEmitter(int numParticles, int particlesPerFrame, Vector2f pos, ParticleShader* shader, float texturePacking, float minAlpha, bool colorize, float initialSize, float sizeOffset, float speed, float speedOffset, float decayPeriod)
{
    this->numParticles = numParticles;
    this->particlesPerFrame = particlesPerFrame;
    this->pos.CopyOf(pos);
    this->texturePacking = texturePacking;
    this->minAlpha = minAlpha;
    this->colorize = colorize;
    this->initialSize = initialSize;        //120.0f
    this->maxSizeOffset = sizeOffset;       //2.8f
    this->speed = speed;                    //2.8f
    this->maxSpeedOffset = speedOffset;     //30.0f
    this->decayPeriod = decayPeriod;
    
    particles = new Particle[numParticles];
    color[0] = colorize ? fRandom() : 1.0f;
    color[1] = colorize ? fRandom() : 1.0f;
    color[2] = colorize ? fRandom() : 1.0f;
    updateFrame = 0.0f;
    this->shader = shader;
    
    Init();
}

ParticleEmitter::~ParticleEmitter()
{
    Uninit();
}

void ParticleEmitter::Render(MatrixF& mvpMatrix)
{
    // Uniforms
    glUniformMatrix4fv(shader->u_mvpMatrix, 1, 0, mvpMatrix.f);
    glUniform1f(shader->u_speed, speed);
    glUniform1f(shader->u_initialSize, initialSize);
    glUniform3f(shader->u_color, color[0], color[1], color[2]);
    glUniform1i(shader->u_sampler, 0);
    glUniform2f(shader->u_position, pos.x, pos.y);
    glUniform1f(shader->u_updateFrame, updateFrame);
    glUniform1f(shader->u_texScale, 1.0f/texturePacking);
    glUniform1f(shader->u_minAlpha, minAlpha);

    pos.x += 0.0f;
    pos.y += 0.0f;

    // Draw particles
    glBindVertexArrayOES(VAO);
    glDrawArrays(GL_POINTS, 0, numParticles);
    glBindVertexArrayOES(0);
}

bool ParticleEmitter::IsActive()
{
    updateFrame += 1.0f;
    if (updateFrame < ((float)numParticles/particlesPerFrame + decayPeriod)) //typically 3 second (180.0f) particle decay looks good
        return true;
    else
        return false;
}

void ParticleEmitter::Init()
{
    float maxColorOffset = colorize ? 0.25f : 0.0f;
    
    //create particles
    for (int i=0; i<numParticles; i++)
    {
        particles[i].ID = (i/particlesPerFrame)*1.0f;
        particles[i].direction = DegToRad(fRandom()*359.95f);
        
        //set up some variations in particle parameters
        particles[i].speedOffset = (fRandom()*2.0f-1.0f)*maxSpeedOffset;
        particles[i].sizeOffset = (fRandom()*2.0f-1.0f)*maxSizeOffset;
        particles[i].colorOffset[0] = (fRandom()*2.0f-1.0f)*maxColorOffset;
        particles[i].colorOffset[1] = (fRandom()*2.0f-1.0f)*maxColorOffset;
        particles[i].colorOffset[2] = (fRandom()*2.0f-1.0f)*maxColorOffset;
        particles[i].texOffset[0] = (int)(fRandom()*texturePacking)/texturePacking;
        particles[i].texOffset[1] = (int)(fRandom()*texturePacking)/texturePacking;
    }
    //set up the VAO
    glGenVertexArraysOES(1, &VAO);
    glBindVertexArrayOES(VAO);
    
    glGenBuffers(1, &vertexBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
    glBufferData(GL_ARRAY_BUFFER, numParticles*sizeof(Particle), particles, GL_STATIC_DRAW);
    
    glEnableVertexAttribArray(shader->a_ID);
    glEnableVertexAttribArray(shader->a_direction);
    glEnableVertexAttribArray(shader->a_speedOffset);
    glEnableVertexAttribArray(shader->a_sizeOffset);
    glEnableVertexAttribArray(shader->a_colorOffset);
    glEnableVertexAttribArray(shader->a_texOffset);
    
    glVertexAttribPointer(shader->a_ID, 1, GL_FLOAT, GL_FALSE, sizeof(Particle), (void*)(offsetof(Particle, ID)));
    glVertexAttribPointer(shader->a_direction, 1, GL_FLOAT, GL_FALSE, sizeof(Particle), (void*)(offsetof(Particle, direction)));
    glVertexAttribPointer(shader->a_speedOffset, 1, GL_FLOAT, GL_FALSE, sizeof(Particle), (void*)(offsetof(Particle, speedOffset)));
    glVertexAttribPointer(shader->a_sizeOffset, 1, GL_FLOAT, GL_FALSE, sizeof(Particle), (void*)(offsetof(Particle, sizeOffset)));
    glVertexAttribPointer(shader->a_colorOffset, 3, GL_FLOAT, GL_FALSE, sizeof(Particle), (void*)(offsetof(Particle, colorOffset)));
    glVertexAttribPointer(shader->a_texOffset, 2, GL_FLOAT, GL_FALSE, sizeof(Particle), (void*)(offsetof(Particle, texOffset)));
    
    glBindVertexArrayOES(0);
}

void ParticleEmitter::Uninit()
{
    glDeleteBuffers(1, &vertexBuffer);
    glDeleteVertexArraysOES(1, &VAO);
}

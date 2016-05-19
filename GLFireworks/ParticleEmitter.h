//
//  ParticleEmitter.h
//  SteerMe
//
//  Created by Admin Dev on 2/26/15.
//  Copyright (c) 2015 Kinetic Bytes. All rights reserved.
//

#ifndef __SteerMe__ParticleEmitter__
#define __SteerMe__ParticleEmitter__

#include <stdlib.h>
#include "ParticleShader.h"
#include "MatrixF.h"
#include "Vector2f.h"

typedef struct Particle
{
    float ID;
    float direction;
    float speedOffset;
    float sizeOffset;
    float colorOffset[3];
    float texOffset[2];
}
Particle;

/////////////////////////////////

class ParticleEmitter
{
public:
    ParticleEmitter(int numParticles, int particlesPerFrame, Vector2f pos, ParticleShader* shader, float texturePacking, float minAlpha, bool colorize, float initialSize, float sizeOffset, float speed, float speedOffset, float decayPeriod);
    virtual ~ParticleEmitter();
    virtual void Render(MatrixF& mvpMatrix);
    virtual bool IsActive();
    
private:
    ParticleEmitter();                         // Don't implement
    ParticleEmitter(ParticleEmitter const&);   // Don't implement
    void operator=(ParticleEmitter const&);    // Don't implement

    virtual void Init();
    virtual void Uninit();
    float fRandom() { return rand()*1.0f/RAND_MAX; }
    
    //data
public:
    
private:
    GLuint VAO;
    Vector2f pos;
    GLuint vertexBuffer;
    
    int numParticles;
    int particlesPerFrame;
    
    Particle* particles;
    float speed;
    float initialSize;
    float color[3];
    float updateFrame;
    float texturePacking; //sub-textures per side of texture atlas (assuming square atlas)
    float minAlpha; //minimum opacity to maintain for the particles as they decay
    bool colorize;
    float maxSizeOffset;
    float maxSpeedOffset;
    float decayPeriod;
    
    ParticleShader* shader;
};

#endif /* defined(__SteerMe__ParticleEmitter__) */

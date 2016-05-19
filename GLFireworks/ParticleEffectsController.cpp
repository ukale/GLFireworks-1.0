//
//  ParticleEffectsController.cpp
//  SteerMe
//
//  Created by Admin Dev on 2/26/15.
//  Copyright (c) 2015 Kinetic Bytes. All rights reserved.
//

#include "ParticleEffectsController.h"

ParticleEffectsController::ParticleEffectsController()
{
    shader = NULL;
}

ParticleEffectsController::~ParticleEffectsController()
{
    Cleanup();
}

void ParticleEffectsController::Cleanup()
{
    //clean up emitters
    int sz = (int)emitters.size();
    for (int i=0; i<sz; i++)
    {
        ParticleEmitter* emitter = emitters[i];
        if (emitter != NULL) delete emitter;
    }
    emitters.clear();
    
    //clean up shader
    if (shader != NULL)
        delete shader;
}

void ParticleEffectsController::AddEmitter(int numParticles, int particlesPerFrame, Vector2f pos, float texturePacking, float minAlpha, bool colorize, float initialSize, float sizeOffset, float speed, float speedOffset, float decayPeriod)
{
    if (shader == NULL)
    {
        //create a single shader for all emitters controlled by this controller
        shader = new ParticleShader();
        shader->Init();
    }
    ParticleEmitter* emitter = new ParticleEmitter(numParticles, particlesPerFrame, pos, shader, texturePacking, minAlpha, colorize, initialSize, sizeOffset, speed, speedOffset, decayPeriod);
    emitters.push_back(emitter);
}

void ParticleEffectsController::RenderEmitters(MatrixF& mvpMatrix)
{
    if (shader == NULL) return;
    glUseProgram(shader->program);
    for (std::vector<ParticleEmitter*>::iterator it = emitters.begin(); it != emitters.end(); ++it)
    {
        ParticleEmitter* emitter = (ParticleEmitter*)(*it);
        if (emitter->IsActive())
            emitter->Render(mvpMatrix);   //only render active emitters
    }
}

void ParticleEffectsController::UpdateLifecycle()
{
    for (std::vector<ParticleEmitter*>::iterator it = emitters.begin(); it != emitters.end();)
    {
        ParticleEmitter* emitter = (ParticleEmitter*)(*it);
        if (!emitter->IsActive())   //check if emitter has expired
        {
            it = emitters.erase(it);
            delete emitter;
        }
        else
            ++it;
    }
}

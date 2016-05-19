//
//  ParticleEffectsController.h
//  SteerMe
//
//  Created by Admin Dev on 2/26/15.
//  Copyright (c) 2015 Kinetic Bytes. All rights reserved.
//

#ifndef __SteerMe__ParticleEffectsController__
#define __SteerMe__ParticleEffectsController__

#include <vector>
#include "ParticleEmitter.h"
#include "ParticleShader.h"

class ParticleEffectsController
{
public:
    ParticleEffectsController();
    virtual ~ParticleEffectsController();
    void AddEmitter(int numParticles, int particlesPerFrame, Vector2f pos, float texturePacking, float minAlpha, bool colorize, float initialSize, float sizeOffset, float speed, float speedOffset, float decayPeriod);
    void Cleanup();
    void UpdateLifecycle();
    void RenderEmitters(MatrixF& mvpMatrix);
    bool HasEmitters() { return emitters.size() > 0; }
    
    //data
protected:
    ParticleShader* shader;
    std::vector<ParticleEmitter*> emitters;
};

#endif /* defined(__SteerMe__ParticleEffectsController__) */

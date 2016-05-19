//
//  MainViewController.m
//  SteerMe
//
//  Created by Admin Dev on 2/26/15.
//  Copyright (c) 2015 Kinetic Bytes. All rights reserved.
//

#import "MainViewController.h"
#import "ParticleEffectsController.h"
#import "MatrixF.h"

/////////////////

@interface MainViewController ()

@property (nonatomic, strong) GLKTextureInfo *particlesTexture;

@end


////////////////////////

@implementation MainViewController


ParticleEffectsController effectsController;


- (void)viewDidLoad
{
    [super viewDidLoad];
    
    // Set up context
    EAGLContext* context = [[EAGLContext alloc] initWithAPI:kEAGLRenderingAPIOpenGLES2];
    [EAGLContext setCurrentContext:context];
    
    self.preferredFramesPerSecond = 60; // 4
    
    // Set up GLKview
    GLKView* view = (GLKView*) self.view;
    view.context = context;
    
    [self createParticlesTexture];
}

- (void)createParticlesTexture
{
    NSError *error;
    self.particlesTexture = [GLKTextureLoader
                             textureWithContentsOfFile:[[NSBundle mainBundle] pathForResource:@"swirl14" ofType:@"png"]
                             options:nil
                             error:&error];
    
    glBindTexture(GL_TEXTURE_2D, self.particlesTexture.name);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
}


#pragma mark - GLKViewDelegate

void ortho(float left, float right, float bottom, float top, float zNear, float zFar, float* f)
{
    float h = 1.0f / (right - left);
    float i = 1.0f / (top - bottom);
    float j = 1.0f / (zFar - zNear);
    
    memset(f, 0, 16*sizeof(float));
    f[0] = 2.0f*h;
    f[5] = 2.0f*i;
    f[10] = -2.0f*j;
    f[12] = -(right + left)*h;
    f[13] = -(top + bottom)*i;
    f[14] = -(zFar + zNear)*j;
    f[15] = 1.0f;
};

- (void)glkView:(GLKView *)view drawInRect:(CGRect)rect
{
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);
    //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    MatrixF m;
    ortho(0.0f, self.view.bounds.size.width, self.view.bounds.size.height, 0.0f, -1.0, 1.0, m.f);
    
    glBindTexture(GL_TEXTURE_2D, self.particlesTexture.name);
    effectsController.RenderEmitters(m);
}

- (void)update
{
    effectsController.UpdateLifecycle();
}

- (void)touchesBegan:(NSSet *)touches withEvent:(UIEvent *)event
{
    // Get touch point and screen information
    CGPoint touchPoint = [touches.anyObject locationInView:self.view];
    //effectsController.AddEmitter(1200, 60, Vector2f(touchPoint.x, touchPoint.y));
    effectsController.AddEmitter(2048, 512, Vector2f(touchPoint.x, touchPoint.y), 32.0f, 1.0f, true, 20.0f, 1.0f, 1.8f, 1.8f, 240.0f);
}
/*
- (void)touchesMoved:(NSSet *)touches withEvent:(UIEvent *)event
{
    // Get touch point and screen information
    CGPoint touchPoint = [touches.anyObject locationInView:self.view];
    effectsController.AddEmitter(600, 60, Vector2f(touchPoint.x, touchPoint.y));
}
*/
@end

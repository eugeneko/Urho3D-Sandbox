Vector3 scale(0.01, 0.01, 0.01);

OverrideModelScale(cache.resourceDirs[0] +
    "Default_Character/Models/Female.mdl",
    "Default_Character/Models/Female.mdl", scale);

OverrideModelScale(cache.resourceDirs[0] +
    "Default_Character/Models/Male.mdl",
    "Default_Character/Models/Male.mdl", scale);

OverrideAnimationScale(cache.resourceDirs[0] +
    "Default_Character/Animations/idle.ani",
    "Default_Character/Animations/idle.ani",
    "Default_Character/Models/Female.mdl", scale);

OverrideAnimationScale(cache.resourceDirs[0] +
    "Default_Character/Animations/walking.ani",
    "Default_Character/Animations/walking.ani",
    "Default_Character/Models/Female.mdl", scale);

OverrideAnimationScale(cache.resourceDirs[0] +
    "Default_Character/Animations/running.ani",
    "Default_Character/Animations/running.ani",
    "Default_Character/Models/Female.mdl", scale);

OverrideAnimationScale(cache.resourceDirs[0] +
    "Default_Character/Animations/walking_backward.ani",
    "Default_Character/Animations/walking_backward.ani",
    "Default_Character/Models/Female.mdl", scale);

OverrideAnimationScale(cache.resourceDirs[0] +
    "Default_Character/Animations/jump.full.ani",
    "Default_Character/Animations/jump.full.ani",
    "Default_Character/Models/Female.mdl", scale);

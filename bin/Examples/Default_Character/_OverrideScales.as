Vector3 scale(0.01, 0.01, 0.01);

OverrideModelScale(cache.resourceDirs[0] +
    "Default_Character/Models/Model.mdl",
    "Default_Character/Models/Model.mdl", scale);

OverrideAnimationScale(cache.resourceDirs[0] +
    "Default_Character/Animations/idle.ani",
    "Default_Character/Animations/idle.ani",
    "Default_Character/Models/Model.mdl", scale);

OverrideAnimationScale(cache.resourceDirs[0] +
    "Default_Character/Animations/walking.ani",
    "Default_Character/Animations/walking.ani",
    "Default_Character/Models/Model.mdl", scale);

OverrideAnimationScale(cache.resourceDirs[0] +
    "Default_Character/Animations/walking_backward.ani",
    "Default_Character/Animations/walking_backward.ani",
    "Default_Character/Models/Model.mdl", scale);

OverrideAnimationScale(cache.resourceDirs[0] +
    "Default_Character/Animations/jump.ani",
    "Default_Character/Animations/jump.ani",
    "Default_Character/Models/Model.mdl", scale);

Matrix3x4 rotation(Vector3(0,0,0), Quaternion(180, Vector3(0,1,0)), Vector3(1,1,1));

ImportCharacterAnimation(cache.resourceDirs[0] +
    "Default_Character/Animations/idle.char.xml",
    "Default_Character/Animations/idle.ani",
    "Default_Character/Models/Model_Skeleton.xml",
    "Default_Character/Models/Model.mdl",
    rotation);

ImportCharacterAnimation(cache.resourceDirs[0] +
    "Default_Character/Animations/walking.char.xml",
    "Default_Character/Animations/walking.ani",
    "Default_Character/Models/Model_Skeleton.xml",
    "Default_Character/Models/Model.mdl",
    rotation);

ImportCharacterAnimation(cache.resourceDirs[0] +
    "Default_Character/Animations/running.char.xml",
    "Default_Character/Animations/running.ani",
    "Default_Character/Models/Model_Skeleton.xml",
    "Default_Character/Models/Model.mdl",
    rotation);

ImportCharacterAnimation(cache.resourceDirs[0] +
    "Default_Character/Animations/walking_backward.char.xml",
    "Default_Character/Animations/walking_backward.ani",
    "Default_Character/Models/Model_Skeleton.xml",
    "Default_Character/Models/Model.mdl",
    rotation);

ImportCharacterAnimation(cache.resourceDirs[0] +
    "Default_Character/Animations/jump.up.char.xml",
    "Default_Character/Animations/jump.up.ani",
    "Default_Character/Models/Model_Skeleton.xml",
    "Default_Character/Models/Model.mdl",
    rotation);

ImportCharacterAnimation(cache.resourceDirs[0] +
    "Default_Character/Animations/jump.down.char.xml",
    "Default_Character/Animations/jump.down.ani",
    "Default_Character/Models/Model_Skeleton.xml",
    "Default_Character/Models/Model.mdl",
    rotation);

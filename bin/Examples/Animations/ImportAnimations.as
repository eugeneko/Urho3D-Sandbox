ImportCharacterAnimation(cache.resourceDirs[0] +
    "Animations/Kachujin_Walk.xml",
    "Animations/Kachujin_Walk.ani",
    "Models/Kachujin_Skeleton.xml",
    "Models/Kachujin/Kachujin.mdl",
    Matrix3x4());

ImportCharacterAnimation(cache.resourceDirs[0] +
    "Animations/Swat_WalkFwd.xml",
    "Animations/Swat_WalkFwd.ani",
    "Models/Swat_Skeleton.xml",
    "Models/Swat.mdl",
    Matrix3x4(Vector3(0, 0, 0), Quaternion(180, Vector3(0, 1, 0)), 1));
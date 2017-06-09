float femaleLeft  = ComputeSegmentLength("Default_Character/Models/Female.mdl", "Default_Character/Models/Female_Skeleton.xml", "LeftFoot");
float femaleRight = ComputeSegmentLength("Default_Character/Models/Female.mdl", "Default_Character/Models/Female_Skeleton.xml", "RightFoot");
Print("Female scale: " + (femaleLeft + femaleRight) / 2);

float maleLeft  = ComputeSegmentLength("Default_Character/Models/Male.mdl", "Default_Character/Models/Male_Skeleton.xml", "LeftFoot");
float maleRight = ComputeSegmentLength("Default_Character/Models/Male.mdl", "Default_Character/Models/Male_Skeleton.xml", "RightFoot");
Print("Male scale: " + (maleLeft + maleRight) / 2);

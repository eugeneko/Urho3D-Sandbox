class Animator : ScriptObject
{
    float footRotationAmount = 0;
    float footRotationGlobal = 0.5;
    float footOffset = 0;
    float rotateAnimation = 0;
    
    void DelayedStart()
    {
        AnimationController@ animController = node.GetComponent("CharacterAnimationController");
        animController.Play("CharacterAnimator/Swat_WalkFwd.ani", 0, true);
    }
    
    void Update(float timeStep)
    {
        AnimationController@ animController = node.GetComponent("CharacterAnimationController");
        CharacterAnimationController@ characterController = node.GetComponent("CharacterAnimationController");
        Node@ groundControl = node.GetChild("control:Ground");
        
        characterController.SetAnimationTransform(Matrix3x4(Quaternion(rotateAnimation, Vector3(0, 1, 0)).rotationMatrix));
        
        characterController.SetTargetTransform("LeftFoot", groundControl.transform);
        characterController.SetTargetTransform("RightFoot", groundControl.transform);

        characterController.SetTargetRotationAmount("LeftFoot", footRotationAmount);
        characterController.SetTargetRotationAmount("RightFoot", footRotationAmount);
        characterController.SetTargetRotationBalance("LeftFoot", footRotationGlobal);
        characterController.SetTargetRotationBalance("RightFoot", footRotationGlobal);
    }
}

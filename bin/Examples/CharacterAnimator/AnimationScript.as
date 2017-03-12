class Animator : ScriptObject
{
    float footRotationAmount = 0;
    float footRotationGlobal = 0.5;
    float footOffset = 0;
    
    void DelayedStart()
    {
        AnimationController@ animController = node.GetComponent("CharacterAnimationController");
        animController.Play("CharacterAnimator/Swat_WalkFwd.ani", 0, true);
    }
    
    void Update(float timeStep)
    {
        CharacterAnimationController@ animController = node.GetComponent("CharacterAnimationController");
        Node@ groundControl = node.GetChild("control:Ground");
        if (groundControl !is null)
        {
            animController.SetTargetTransform("LeftFoot", groundControl.transform);
            animController.SetTargetTransform("RightFoot", groundControl.transform);
        }

        animController.SetTargetRotationAmount("LeftFoot", footRotationAmount);
        animController.SetTargetRotationAmount("RightFoot", footRotationAmount);
        animController.SetTargetRotationBalance("LeftFoot", footRotationGlobal);
        animController.SetTargetRotationBalance("RightFoot", footRotationGlobal);
    }
}

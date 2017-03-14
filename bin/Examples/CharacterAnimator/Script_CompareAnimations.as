class Animator : ScriptObject
{    
    void DelayedStart()
    {
        AnimationController@ animController = node.GetComponent("AnimationController");
        if (animController is null)
            animController = node.GetComponent("CharacterAnimationController");
        animController.Play("CharacterAnimator/Swat_WalkFwd.ani", 0, true);
    }
}

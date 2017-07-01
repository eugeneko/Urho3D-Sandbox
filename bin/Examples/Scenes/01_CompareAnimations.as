class Animator : ScriptObject
{
    float rotation = 0;
    void DelayedStart()
    {
        CharacterAnimationController@ characterController = node.GetComponent("CharacterAnimationController");
        AnimationController@ animController = node.GetComponent("AnimationController");
        if (animController is null)
            animController = characterController;

        animController.Play("Default_Character/Animations/walking.ani", 0, true);
    }
}

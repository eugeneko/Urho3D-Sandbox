class Animator : ScriptObject
{
    float rotation = 0;
    void DelayedStart()
    {
        AnimationController@ animController = node.GetComponent("AnimationController");
        if (animController is null)
            animController = node.GetComponent("CharacterAnimationController");
        animController.Play("CharacterAnimator/Kachujin_Walk.ani", 0, true);

        CharacterAnimationController@ characterController = node.GetComponent("CharacterAnimationController");
        if (characterController !is null)
        {
            characterController.SetAnimationTransform(Matrix3x4(Quaternion(rotation, Vector3(0, 1, 0)).rotationMatrix));
        }
    }
}

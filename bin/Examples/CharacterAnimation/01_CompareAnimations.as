class Animator : ScriptObject
{
    float rotation = 0;
    void DelayedStart()
    {
        CharacterAnimationController@ characterController = node.GetComponent("CharacterAnimationController");
        AnimationController@ animController = node.GetComponent("AnimationController");
        if (animController is null)
            animController = characterController;

        animController.Play("Animations/Swat_WalkFwd.ani", 0, true);
        //if (characterController !is null)
            //characterController.SetAnimationTransform(Matrix3x4(Quaternion(rotation, Vector3(0, 1, 0)).rotationMatrix));
    }
}

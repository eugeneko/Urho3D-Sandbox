const int CTRL_FORWARD = 1;
const int CTRL_BACK = 2;
const int CTRL_LEFT = 4;
const int CTRL_RIGHT = 8;
const int CTRL_JUMP = 16;

enum State
{
    Initial,
    Idle,
    Walk,
    Jump,
    Fall
}

class Controller
{
    // Configuration
    float _thresholdX = 0.1;
    float _idleThreshold = 0.1;
    float _switchDuration = 0.2;
    float _animationRotationY = 180;
    String _walkAnimation;
    float _walkBaseVelocity = 1.5;
    String _idleAnimation;

    // External state
    float _velocityX = 0;
    float _velocityY = 0;

    // Internal state
    State _state = Initial;
    float _idleTimer = 0;
    bool _aboutToJump = false;
    
    float velocityX
    {
        get const { return _velocityX; }
        set { _velocityX = value; }
    }
    String idleAnimation
    {
        get const { return _idleAnimation; }
        set { _idleAnimation = value; }
    }
    String walkAnimation
    {
        get const { return _walkAnimation; }
        set { _walkAnimation = value; }
    }
    
    void Jump()
    {
        _aboutToJump = true;
    }
    void Update(CharacterAnimationController@ characterController, float timeStep)
    {
        AnimationController@ animController = characterController;
        bool movingX = _velocityX > _thresholdX;

        if (movingX)
            _idleTimer = 0;

        if (_state == Initial)
        {
            _state = Idle;
            animController.PlayExclusive(_idleAnimation, 0, true);
        }
        else if (_state == Idle)
        {
            if (movingX)
            {
                _state = Walk;
                animController.PlayExclusive(_walkAnimation, 0, true, _switchDuration);
            }
        }
        else if (_state == Walk)
        {
            animController.SetSpeed(_walkAnimation, Abs(_velocityX / _walkBaseVelocity));
            if (!movingX)
            {
                _idleTimer += timeStep;
                if (_idleTimer >= _idleThreshold)
                {
                    _idleTimer = 0;
                    _state = Idle;
                    animController.PlayExclusive(_idleAnimation, 0, true, _switchDuration);
                }
            }
        }
    }
}

class Main : ScriptObject
{
    Controls _controls;
    Vector3 _previousPosition;
    Vector3 _velocity;
    Controller@ _controller;
    void DelayedStart()
    {
        @_controller = Controller();
        _controller.walkAnimation = "Animations/Kachujin_Walk.ani";
        _controller.idleAnimation = "Animations/Kachujin_Idle.ani";
        _previousPosition = node.worldPosition;
    }
    void FixedUpdate(float timeStep)
    {
        const float walkSpeed = 1.8;
        float walkDirection = 0;
        if (_controls.IsDown(CTRL_LEFT))
            walkDirection += -1;
        if (_controls.IsDown(CTRL_RIGHT))
            walkDirection += 1;
        _previousPosition = node.worldPosition;
        node.worldPosition = node.worldPosition + Vector3(1, 0, 0) * timeStep * walkSpeed * walkDirection;
        _velocity = (node.worldPosition - _previousPosition) / timeStep;
        _controller.velocityX = _velocity.x;
    }
    void Update(float timeStep)
    {
        _controls.Set(CTRL_FORWARD | CTRL_BACK | CTRL_LEFT | CTRL_RIGHT | CTRL_JUMP, false);
        if (ui.focusElement is null)
        {
            _controls.Set(CTRL_FORWARD, input.keyDown[KEY_I]);
            _controls.Set(CTRL_BACK, input.keyDown[KEY_K]);
            _controls.Set(CTRL_LEFT, input.keyDown[KEY_J]);
            _controls.Set(CTRL_RIGHT, input.keyDown[KEY_L]);
            _controls.Set(CTRL_JUMP, input.keyDown[KEY_SPACE]);

            _controls.yaw = input.mousePosition.x;
            _controls.pitch = input.mousePosition.y;
        }

        CharacterAnimationController@ characterController = node.GetComponent("CharacterAnimationController");
        _controller.Update(characterController, timeStep);
    }
}

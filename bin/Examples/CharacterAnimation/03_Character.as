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
    void Update(CharacterAnimationController@ characterController, float deltaTime)
    {
        characterController.SetAnimationTransform(Matrix3x4(Quaternion(_animationRotationY, Vector3(0, 1, 0)).rotationMatrix));
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
                _idleTimer += deltaTime;
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

class MockPlayer
{
    float _time = 0.0;
    Controller@ _controller;
    
    MockPlayer(Controller@ controller, Node@ node)
    {
        @_controller = @controller;        
        _controller.walkAnimation = "Animations/Kachujin_Walk.ani";
        _controller.idleAnimation = "Animations/Kachujin_Idle.ani";
    }
    void Update(float deltaTime, Node@ node)
    {
        _time += deltaTime;
        if (_time < 1.0)
        {
            _controller.velocityX = 0;
        }
        else if (_time < 5.0)
        {
            _controller.velocityX = 1;
            node.position = node.position + Vector3(0, 0, 1) * deltaTime;
        }
        else if (_time < 6.0)
        {
            _controller.velocityX = 0;
        }
        else
        {
            _time = 0;
            _controller.velocityX = 0;
            node.position = Vector3(0, 0, 0);
        }

        CharacterAnimationController@ characterController = node.GetComponent("CharacterAnimationController");
        _controller.Update(characterController, deltaTime);
    }
}

class Main : ScriptObject
{
    Controller@ _controller;
    MockPlayer@ _player;
    void DelayedStart()
    {
        @_controller = Controller();
        _player = MockPlayer(_controller, node);
    }
    void Stop()
    {
        @_controller = null;
        @_player = null;
    }
    void Update(float deltaTime)
    {
        _player.Update(deltaTime, node);
    }
}

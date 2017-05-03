enum State
{
    Idle,
    Walk,
    Jump,
    Fall
}

class Controller
{
    float _thresholdX = 0.1;
    State _state = Idle;
    float _velocityX = 0;
    float _velocityY = 0;
    bool _aboutToJump = false;
    
    float velocityX
    {
        get const { return _velocityX; }
        set { _velocityX = value; }
    }
    void Jump()
    {
        _aboutToJump = true;
    }
    void UpdateState(CharacterAnimationController@ characterController)
    {
        AnimationController@ animController = characterController;
        if (_state == Idle)
        {
            if (_velocityX > _thresholdX)
            {
                _state = Walk;
            }
        }
        else if (_state == Walk)
        {
            if (_velocityX <= _thresholdX)
            {
                _state = Idle;
            }
        }
    }
}

class MockPlayer
{
    float _time = 0.0;
    Controller@ _controller;
    Node@ _node;
    CharacterAnimationController@ _characterController;
    AnimationController@ _animController;
    
    MockPlayer(Controller@ controller, Node@ node)
    {
        @_controller = @controller;
        _node = node;
        _characterController = node.GetComponent("CharacterAnimationController");
        _animController = _characterController;
    }
    void Update(float deltaTime)
    {
        _time += deltaTime;
        if (_time < 1.0)
        {
            _controller.velocityX = 0;
        }
        else if (_time < 10.0)
        {
            _controller.velocityX = 1;
            _node.position = _node.position + Vector3(0, 0, 1) * deltaTime;
        }
        else if (_time < 11.0)
        {
            _controller.velocityX = 0;
        }
        else
        {
            _time = 0;
            _controller.velocityX = 0;
            _node.position = Vector3(0, 0, 0);
        }
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
    void Update(float deltaTime)
    {
        _player.Update(deltaTime);
    }
}

const int CTRL_UP = 1;
const int CTRL_DOWN = 2;
const int CTRL_LEFT = 4;
const int CTRL_RIGHT = 8;
const int CTRL_SLOW = 16;

enum State
{
    Idle,
    Walk,
    WalkBackward,
    Jump,
    Fall,
    State_COUNT
}

class WalkAnimationDesc
{
    String animationName;
    Animation@ animation;
    float minVelocity;
    float baseVelocity;
    float footScale;
    float length;
}

class Animator
{
    // Configuration
    private Array<String> _animations; ///< Animations for each state.
    private Array<WalkAnimationDesc> _walkAnimations; ///< Walk animations.
    
    float modelFootScale = 1;
    float idleThreshold = 0.1;
    float walkSwitchDuration = 0.2;
    float jumpSwitchDuration = 0.1;
    float flySwitchDuration = 0.4;
    private float _switchDuration = 0.2;
    
    /// Add walk animation. Shall be ordered from negative to positive speeds.
    void AddWalkAnimation(String &in animation, float minVelocity)
    {
        WalkAnimationDesc desc;
        desc.animationName = animation;
        desc.animation = cache.GetResource("Animation", animation);
        desc.minVelocity = minVelocity;
        desc.baseVelocity = desc.animation.metadata["speed"].GetFloat();
        desc.footScale = desc.animation.metadata["footscale"].GetFloat();
        desc.length = desc.animation.length;
        _walkAnimations.Push(desc);
    }
    /// Add animation. Don't use this for walk animations.
    void AddAnimation(State state, String &in animation)
    {
        _animations[state] = animation;
    }

    // Control variables
    bool grounded = false; ///< Whether the character is on the ground.
    bool jump = false; ///< Whether the character is about to jump.
    float movementSpeed = 0; ///< Horizontal movement speed.
    float verticalSpeed = 0; ///< Vertical movement speed.

    // Internal state
    State _state = Idle;
    float _idleTimer = 0;
    
    /// Construct.
    Animator()
    {
        _animations.Resize(State_COUNT);
        _switchDuration = walkSwitchDuration;
    }
    private float GetWalkPhase(AnimationController@ animController) const
    {
        for (uint i = 0; i < _walkAnimations.length; ++i)
        {
            if (animController.IsPlaying(_walkAnimations[i].animationName))
                return animController.GetTime(_walkAnimations[i].animationName) / _walkAnimations[i].length;
        }
        return -1;
    }
    private uint SelectWalkAnimation(float velocity) const
    {
        for (uint i = 0; i < _walkAnimations.length; ++i)
        {
            float edge = _walkAnimations[i].minVelocity;
            if (edge > 0)
                break;
            if (velocity <= edge)
                return i;
        }
        for (uint i = _walkAnimations.length; i > 0; --i)
        {
            float edge = _walkAnimations[i - 1].minVelocity;
            if (edge < 0)
                break;
            if (velocity >= edge)
                return i - 1;
        }
        return 0;
    }
    private void UpdateWalkAnimation(AnimationController@ animController, float velocity, float switchDuration)
    {
        uint idx = SelectWalkAnimation(velocity);
        const String animationName = _walkAnimations[idx].animationName;
        bool isPlaying = animController.IsPlaying(animationName);
        float phase = GetWalkPhase(animController);
        animController.PlayExclusive(animationName, 0, true, switchDuration);
        if (!isPlaying)
            animController.SetTime(animationName, phase * _walkAnimations[idx].length);

        float baseVelocity = _walkAnimations[idx].baseVelocity * modelFootScale / _walkAnimations[idx].footScale;
        animController.SetSpeed(animationName, Abs(velocity / baseVelocity));
    }
    void Update(CharacterAnimationController@ characterController, float timeStep)
    {
        AnimationController@ animController = characterController;
        
        // Compute parameters
        bool movingX = Abs(movementSpeed) > 0;
        if (movingX)
            _idleTimer = 0;

        // Update state
        if (grounded)
        {
            // Jump if asked for
            if (jump)
            {
                _state = Jump;
                jump = false;
            }
            // Start movement
            else if (movingX)
            {
                _state = Walk;
            }
            // Go idle if stay for a while
            else
            {
                _idleTimer += timeStep;
                if (_idleTimer >= idleThreshold)
                {
                    _idleTimer = 0;
                    _state = Idle;
                }
            }
        }
        else
        {
            // Convert jump to fall
            if (_state == Jump && verticalSpeed <= 0)
            {
                _state = Fall;
            }
        }
        
        // Apply animation
        switch (_state)
        {
        case Idle:
            animController.PlayExclusive(_animations[Idle], 0, true, _switchDuration);
            _switchDuration = walkSwitchDuration;
            break;
        case Walk:
            UpdateWalkAnimation(animController, movementSpeed, _switchDuration);
            _switchDuration = walkSwitchDuration;
            break;
        case Jump:
            animController.PlayExclusive(_animations[Jump], 0, false, jumpSwitchDuration);
            animController.SetSpeed(_animations[Jump], 0.5);
            _switchDuration = jumpSwitchDuration;
            break;
        case Fall:
            animController.PlayExclusive(_animations[Fall], 0, false, flySwitchDuration);
            break;
        }
    }
}

class Controller
{
    // Grounding configuration
    float flyThreshold = 0.1;
    float jumpCooldown = 0.5;
    
    // Movement parameters
    float moveVelocity = 1;
    float jumpVelocity = 6;

    // Rotation parameters
    float flipDuration = 0.3;       ///< Flip from -1 to 1 duration.
    float rotationNegative = 90;    ///< Rotation angle for negative direction.
    float rotationNeutral = 0;      ///< Rotation angle for zero direction.
    float rotationPositive = -90;   ///< Rotation angle for positive direction.
    
    // Control variables
    int moveDirection = 0;          ///< Movement direction.
    bool slow = false;              ///< Whether the character is in slow movement mode.
    bool grounded = false;          ///< Whether the character is on the ground.
    bool jump = false;              ///< Whether the character is about to jump.
    int lookDirection = 1;          ///< Look direction.
    
    ///  90   -90
    ///  <  0  >
    ///     v
    /// -1  0  1
    private int _direction = 1;
    private float _currentDirection = 1;

    private float _inAirTimer = 0;
    private float _jumpTimer = 0;
    private bool _softGrounded = false;
    
    private void UpdateParameters(float timeStep)
    {
        // Update timers
        if (!grounded)
            _inAirTimer += timeStep;
        else
            _inAirTimer = 0.0f;
        _jumpTimer += timeStep;
        _softGrounded = _jumpTimer > jumpCooldown && _inAirTimer < flyThreshold;

        // Update jump state
        if (jump)
        {
            if (!_softGrounded)
                jump = false;
            else
                _jumpTimer = 0;
        }
    }
    void Update(Node@ node, RigidBody@ rigidBody, Animator@ animator, float timeStep)
    {
        // Update parameters
        UpdateParameters(timeStep);
        animator.grounded = _softGrounded;

        // Update direction and speed
        bool isMoving = moveDirection != 0;
        if (slow)
            _direction = lookDirection;
        else if (isMoving)
            _direction = moveDirection >= 0 ? 1 : -1;

        // Update jump
        if (jump)
        {
            rigidBody.ApplyImpulse(Vector3(0, 1, 0) * rigidBody.mass * jumpVelocity);
            animator.jump = true;
            jump = false;
        }

        // Apply movement
        float speed = moveVelocity * moveDirection;
        animator.movementSpeed = speed * _direction;
        Vector3 linearVelocity = rigidBody.linearVelocity;
        linearVelocity.x = speed;
        rigidBody.linearVelocity = linearVelocity;
        animator.verticalSpeed = linearVelocity.y;
        
        // Interpolate direction and apply node rotation
        float flipSpeed = 1 / (0.5 * flipDuration);
        _currentDirection = Clamp(_currentDirection + Sign(_direction - _currentDirection) * flipSpeed * timeStep, -1.0f, 1.0f);
        float angle = _currentDirection < 0
            ? Lerp(rotationNeutral, rotationNegative, -_currentDirection)
            : Lerp(rotationNeutral, rotationPositive, _currentDirection);
        node.worldRotation = Quaternion(angle, Vector3(0, 1, 0));
    }
}

class Main : ScriptObject
{
    float walkVelocity = 1.5;
    float runVelocity = 4;
    bool male = false;

    private float _modelFootScale = 1;
    private Controls _controls;
    private Animator@ _animator;
    private Controller@ _controller;
    float _velocity;

    void DelayedStart()
    {
        SubscribeToEvent(node, "NodeCollision", "HandleNodeCollision");
        SubscribeToEvent("KeyDown", "HandleKeyDown");

        @_animator = Animator();
        @_controller = Controller();
        _animator.AddAnimation(Idle, "Default_Character/Animations/idle.ani");
        _animator.AddAnimation(Jump, "Default_Character/Animations/jump.up.ani");
        _animator.AddAnimation(Fall, "Default_Character/Animations/jump.down.ani");
        _animator.AddWalkAnimation("Default_Character/Animations/walking_backward.ani", -1);
        _animator.AddWalkAnimation("Default_Character/Animations/walking.ani", 1);
        _animator.AddWalkAnimation("Default_Character/Animations/running.ani", 3);
    }
    void HandleNodeCollision(StringHash eventType, VariantMap& eventData)
    {
        VectorBuffer contacts = eventData["Contacts"].GetBuffer();

        while (!contacts.eof)
        {
            Vector3 contactPosition = contacts.ReadVector3();
            Vector3 contactNormal = contacts.ReadVector3();
            float contactDistance = contacts.ReadFloat();
            float contactImpulse = contacts.ReadFloat();

            // If contact is below node center and pointing up, assume it's a ground contact
            if (contactPosition.y < (node.position.y + 1.0f))
            {
                float level = contactNormal.y;
                if (level > 0.75)
                    _controller.grounded = true;
            }
        }
    }
    void HandleKeyDown(StringHash eventType, VariantMap& eventData)
    {
        int key = eventData["Key"].GetInt();
        if (key == KEY_H)
            _controls.Set(CTRL_SLOW, !_controls.IsDown(CTRL_SLOW));
    }
    void FixedUpdate(float timeStep)
    {
        int moveDirection = 0;
        if (_controls.IsDown(CTRL_LEFT))
            moveDirection += -1;
        if (_controls.IsDown(CTRL_RIGHT))
            moveDirection += 1;

        CharacterAnimationController@ characterController = node.GetComponent("CharacterAnimationController");
        RigidBody@ rigidBody = node.GetComponent("RigidBody");
        _controller.moveDirection = moveDirection;
        _controller.slow = _controls.IsDown(CTRL_SLOW);
        _controller.jump = _controls.IsDown(CTRL_UP);
        _controller.lookDirection = _controls.yaw >= graphics.width / 2 ? 1 : -1;
        _controller.moveVelocity = _controls.IsDown(CTRL_SLOW) ? walkVelocity : runVelocity;
        //_controller.aim = Vector2(_controls.yaw, _controls.pitch);
        _controller.Update(node, rigidBody, _animator, timeStep);
        _animator.modelFootScale = _modelFootScale;
        _animator.Update(characterController, timeStep);
        
        _controller.grounded = false;
    }
    void ApplyAttributes()
    {
        String modelName = male ? "Default_Character/Models/Male.mdl" : "Default_Character/Models/Female.mdl";
        Model@ model = cache.GetResource("Model", modelName);

        AnimatedModel@ animatedModel = node.GetComponent("AnimatedModel");
        animatedModel.model = model;
        animatedModel.ApplyMaterialList();
        _modelFootScale = model.metadata["footscale"].GetFloat();
    }
    void Update(float timeStep)
    {
        _controls.Set(CTRL_UP | CTRL_DOWN | CTRL_LEFT | CTRL_RIGHT, false);
        if (ui.focusElement is null)
        {
            _controls.Set(CTRL_UP, input.keyDown[KEY_I]);
            _controls.Set(CTRL_DOWN, input.keyDown[KEY_K]);
            _controls.Set(CTRL_LEFT, input.keyDown[KEY_J]);
            _controls.Set(CTRL_RIGHT, input.keyDown[KEY_L]);

            _controls.yaw = input.mousePosition.x;
            _controls.pitch = input.mousePosition.y;
        }
    }
}

const int CTRL_UP = 1;
const int CTRL_DOWN = 2;
const int CTRL_LEFT = 4;
const int CTRL_RIGHT = 8;
const int CTRL_JUMP = 16;
const int CTRL_SLOW = 32;

enum State
{
    Idle,
    Walk,
    Jump,
    Fall
}

class Animator
{
    // Configuration
    float _idleThreshold = 0.1;
    float _switchDuration = 0.2;
    float _animationRotationY = 180;
    String _walkAnimation;
    float _walkBaseVelocity = 1.5;
    String _idleAnimation;
    String _jumpAnimation;

    String idleAnimation    { get const { return _idleAnimation; } set { _idleAnimation = value; } }
    String walkAnimation    { get const { return _walkAnimation; } set { _walkAnimation = value; } }
    String jumpAnimation    { get const { return _jumpAnimation; } set { _jumpAnimation = value; } }

    // Control variables
    bool _grounded = false; ///< Whether the character is on the ground.
    bool _jump = false; ///< Whether the character is about to jump.
    float _movementSpeed = 0; ///< Horizontal movement speed.

    float movementSpeed { get const { return _movementSpeed; }  set { _movementSpeed = value; } }
    bool jump           { get const { return _jump; }           set { _jump = value; } }
    bool grounded       { get const { return _grounded; }       set { _grounded = value; } }

    // Internal state
    State _state = Idle;
    float _idleTimer = 0;
    
    void Update(CharacterAnimationController@ characterController, float timeStep)
    {
        AnimationController@ animController = characterController;
        bool movingX = _movementSpeed > 0;

        if (movingX)
            _idleTimer = 0;

        // Jump
        if (jump && grounded)
        {
            _state = Jump;
            jump = false;
            animController.PlayExclusive(_jumpAnimation, 0, false, _switchDuration);
            animController.SetSpeed(_jumpAnimation, 0.25);
        }
        else // Update movement
        if ((_state == Idle || _state == Jump) && grounded)
        {
            //Print("Go Walk");
            if (movingX)
            {
                _state = Walk;
                animController.PlayExclusive(_walkAnimation, 0, true, _switchDuration);
            }
            else
                animController.PlayExclusive(_idleAnimation, 0, true, _switchDuration);
        }
        else if ((_state == Walk || _state == Jump) && grounded)
        {
            //Print("Go Idle");
            animController.SetSpeed(_walkAnimation, Abs(_movementSpeed / _walkBaseVelocity));
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
            else
                animController.PlayExclusive(_walkAnimation, 0, true, _switchDuration);
        }
    }
}

class Controller
{
    float _moveThreshold = 0.01;    ///< Movement threshold.

    // Rotation parameters
    float _flipDuration = 0.3;      ///< Flip from -1 to 1 duration.
    float _rotationNegative = 90;   ///< Rotation angle for negative direction.
    float _rotationNeutral = 0;     ///< Rotation angle for zero direction.
    float _rotationPositive = -90;  ///< Rotation angle for positive direction.
    
    float _speed = 0;
    bool _slow = false;
    bool _grounded = false;
    bool _jump = false;
    
    float _inAirTimer = 0;
    float _jumpCooldown = 0;

    /// Horizontal character speed.
    float speed     { get const { return _speed;    }  set { _speed = value;    } }
    /// Whether the character is in slow movement mode.
    bool slow       { get const { return _slow;     }  set { _slow = value;     } }
    /// Whether the character is on the ground.
    bool grounded   { get const { return _grounded; }  set { _grounded = value; } }
    /// Whether the character is about to jump.
    bool jump       { get const { return _jump;     }  set { _jump = value;     } }

    ///  90   -90
    ///  <  0  >
    ///     v
    /// -1  0  1
    int _direction = 1;
    float _currentDirection = 1;
    
    void Update(Node@ node, Animator@ animator, float timeStep)
    {
        RigidBody@ rigidBody = node.GetComponent("RigidBody");
        
        // Update the in air timer. Reset if grounded
        if (!grounded)
            _inAirTimer += timeStep;
        else
            _inAirTimer = 0.0f;
        _jumpCooldown = Max(_jumpCooldown - timeStep, 0.0);
        // When character has been in air less than 1/10 second, it's still interpreted as being on ground
        bool isGrounded = _inAirTimer < 0.1;
        bool canJump = _jumpCooldown == 0 && isGrounded;
        animator.grounded = canJump;

        // Update movement speed
        bool isMoving = Abs(speed) > _moveThreshold;
        if (isMoving)
        {
            if (slow)
            {
                animator.movementSpeed = speed * _direction;
            }
            else
            {
                _direction = speed >= 0 ? 1 : -1;
                animator.movementSpeed = Abs(speed * _direction);
            }
        }
        else
        {
            animator.movementSpeed = 0.0;
        }
        
        // Update current direction and node rotation
        float flipSpeed = 1 / (0.5 * _flipDuration);
        _currentDirection = Clamp(_currentDirection + Sign(_direction - _currentDirection) * flipSpeed * timeStep, -1.0f, 1.0f);
        float angle = _currentDirection < 0
            ? Lerp(_rotationNeutral, _rotationNegative, -_currentDirection)
            : Lerp(_rotationNeutral, _rotationPositive, _currentDirection);
        node.worldRotation = Quaternion(angle, Vector3(0, 1, 0));
        
        // Apply physics
        if (jump && canJump)
        {
            rigidBody.ApplyImpulse(Vector3(0, 1, 0) * rigidBody.mass * 6);
            animator.jump = true;
            jump = false;
            _jumpCooldown = 0.5;
        }
        
        Vector3 linearVelocity = rigidBody.linearVelocity;
        linearVelocity.x = speed;
        rigidBody.linearVelocity = linearVelocity;
    }
}

class Main : ScriptObject
{
    Controls _controls;
    Animator@ _animator;
    Controller@ _controller;
    void DelayedStart()
    {
        SubscribeToEvent(node, "NodeCollision", "HandleNodeCollision");

        @_animator = Animator();
        @_controller = Controller();
        //_animator.walkAnimation = "Doll_Female/Animations/walking_backward.ani";
        //_animator.walkAnimation = "Doll_Female/Animations/running.ani";
        //_animator.walkAnimation = "Doll_Female/Animations/jogging.ani";
        _animator.walkAnimation = "Default_Character/Animations/walking.ani";
        //_animator.walkAnimation = "Models/Mutant/Mutant_Run.ani";
        //_animator.walkAnimation = "Doll_Female/Animations/idle.ani";
        //_animator.idleAnimation = "Doll_Female/Animations/walking.ani";
        _animator.idleAnimation = "Default_Character/Animations/idle.ani";
        _animator.jumpAnimation = "Default_Character/Animations/jump.ani";
        //_animator.idleAnimation = "Models/Mutant/Mutant_Run.ani";
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
    void FixedUpdate(float timeStep)
    {
        const float walkSpeed = 1.5;
        //const float walkSpeed = 2;
        float walkDirection = 0;
        if (_controls.IsDown(CTRL_LEFT))
            walkDirection += -1;
        if (_controls.IsDown(CTRL_RIGHT))
            walkDirection += 1;

        CharacterAnimationController@ characterController = node.GetComponent("CharacterAnimationController");
        RigidBody@ rigidBody = node.GetComponent("RigidBody");
        _controller.speed = walkSpeed * walkDirection;
        _controller.slow = _controls.IsDown(CTRL_SLOW);
        _controller.jump = _controls.IsDown(CTRL_UP);
        _controller.Update(node, _animator, timeStep);
        _animator.Update(characterController, timeStep);
        
        _controller.grounded = false;
    }
    void Update(float timeStep)
    {
        _controls.Set(CTRL_UP | CTRL_DOWN | CTRL_LEFT | CTRL_RIGHT | CTRL_JUMP | CTRL_SLOW, false);
        if (ui.focusElement is null)
        {
            _controls.Set(CTRL_UP, input.keyDown[KEY_I]);
            _controls.Set(CTRL_DOWN, input.keyDown[KEY_K]);
            _controls.Set(CTRL_LEFT, input.keyDown[KEY_J]);
            _controls.Set(CTRL_RIGHT, input.keyDown[KEY_L]);
            _controls.Set(CTRL_SLOW, input.keyDown[KEY_U]);
            _controls.Set(CTRL_JUMP, input.keyDown[KEY_SPACE]);

            _controls.yaw = input.mousePosition.x;
            _controls.pitch = input.mousePosition.y;
        }

        if (node.worldPosition.x > 5)
            node.worldPosition = Vector3(0, 0, 0);
    }
}

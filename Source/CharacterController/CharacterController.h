#pragma once

// #include <Urho3D/Scene/Component.h>
// #include <Urho3D/Container/ArrayPtr.h>
// #include <Urho3D/Graphics/AnimationController.h>
// #include <Urho3D/Resource/Resource.h>
// #include <Urho3D/Resource/XMLFile.h>
#include <Urho3D/Scene/LogicComponent.h>
// #include <Urho3D/Scene/Node.h>
// #include <memory>

class asIScriptEngine;

namespace Urho3D
{

class RigidBody;

// class Animation;
// class AnimatedModel;
// class Model;

}

namespace Urho3D
{

/// Character Controller.
class CharacterController : public LogicComponent
{
    URHO3D_OBJECT(CharacterController, LogicComponent);

public:
    using Component::SetNode;
    /// Construct.
    CharacterController(Context* context);
    /// Register object factory.
    static void RegisterObject(Context* context);

    /// @see LogicComponent::FixedUpdate
    virtual void FixedUpdate(float timeStep) override;

    /// Set jump velocity.
    void SetJumpVelocity(float velocity) { jumpVelocity_ = velocity; }
    /// Return jump velocity.
    float GetJumpVelocity() const { return jumpVelocity_; }

    /// Set whether the jump is requested.
    void SetJump(bool wantToJump) { wantToJump_ = wantToJump; }
    /// Set controller velocity in XZ plane.
    void SetVelocity(const Vector3& velocity) { controllerVelocity_ = velocity; }

    /// Return whether the character has jumped.
    bool HasJumped() const { return jumped_; }
    /// Return grounded flag.
    bool IsGrounded() const { return softGrounded_; }
    /// Return actual movement velocity.
    const Vector3& GetVelocity() const { return linearVelocity_; }

private:
    /// Handle node set.
    virtual void OnNodeSet(Node* node) override;
    /// Handle collision.
    void HandleCollsion(StringHash eventType, VariantMap& eventData);
    /// Check integrity of the component.
    bool CheckIntegrity();

private:
    /// Rigid body of character.
    WeakPtr<RigidBody> body_;

    /// Grounded threshold.
    float groundedThreshold_ = 0.1f;
    /// Jump cool down.
    float jumpCooldown_ = 0.5f;
    /// Jump velocity.
    float jumpVelocity_ = 6.0f;
    /// Max slope.
    float maxSlope_ = 70.0f;

    /// Whether the body is grounded.
    bool grounded_ = false;
    /// Whether the jump is requested.
    bool wantToJump_ = false;
    /// Movement velocity.
    Vector3 controllerVelocity_;

    /// Fly duration.
    float flyDuration_ = 0.0f;
    /// Jump duration.
    float jumpDuration_ = 0.0f;
    /// Contact normal.
    Vector3 contactNormal_;

    /// Whether the jump was just performed.
    bool jumped_ = false;
    /// Whether the character is grounded considering the threshold.
    bool softGrounded_ = false;
    /// Linear velocity of the body.
    Vector3 linearVelocity_;
};

/// Register classes.
void RegisterCharacterController(Context* context);
/// Register script API.
void RegisterCharacterControllerScriptAPI(asIScriptEngine* engine);

}

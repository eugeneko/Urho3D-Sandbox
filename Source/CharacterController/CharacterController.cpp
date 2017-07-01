#include "CharacterController.h"

#include <Urho3D/AngelScript/APITemplates.h>
#include <Urho3D/IO/MemoryBuffer.h>
#include <Urho3D/Physics/PhysicsEvents.h>
#include <Urho3D/Physics/RigidBody.h>
// #include <Urho3D/Container/Ptr.h>
// #include <Urho3D/Core/Context.h>
// #include <Urho3D/Graphics/AnimatedModel.h>
// #include <Urho3D/Graphics/Animation.h>
// #include <Urho3D/Graphics/AnimationController.h>
// #include <Urho3D/Graphics/AnimationState.h>
// #include <Urho3D/Graphics/DebugRenderer.h>
// #include <Urho3D/IO/FileSystem.h>
// #include <Urho3D/IO/Log.h>
// #include <Urho3D/Math/MathDefs.h>
// #include <Urho3D/Math/Sphere.h>
// #include <Urho3D/Resource/ResourceCache.h>
// #include <Urho3D/Scene/Serializable.h>
//
// #include <algorithm>

namespace Urho3D
{

extern const char* PHYSICS_CATEGORY;

CharacterController::CharacterController(Context* context)
    : LogicComponent(context)
{

}

void CharacterController::RegisterObject(Context* context)
{
    context->RegisterFactory<CharacterController>();
}

void CharacterController::FixedUpdate(float timeStep)
{
    if (!CheckIntegrity())
        return;

    linearVelocity_ = body_->GetLinearVelocity();
    currentJumpVelocity_ = Vector3(controllerVelocity_.x_, jumpVelocity_, controllerVelocity_.z_);

    // Use default contact normal if cannot compute real, normalize otherwise
    const bool contactNormalValid = contactNormal_.LengthSquared() > M_EPSILON;
    if (!contactNormalValid)
        contactNormal_ = Vector3::UP;
    else
        contactNormal_.Normalize();

    // Compute grounded state according to threshold
    grounded_ = contactNormalValid && contactNormal_.y_ > Cos(maxSlope_);
    if (!grounded_)
        flyDuration_ += timeStep;
    else
        flyDuration_ = 0.0f;
    softGrounded_ = jumpDuration_ > jumpCooldown_ && flyDuration_ < groundedThreshold_;

    // Update jump
    jumpDuration_ += timeStep;
    jumped_ = wantToJump_ && softGrounded_;
    wantToJump_ = false;
    if (jumped_)
    {
        jumpDuration_ = 0;
        softGrounded_ = false;
    }

    // Update body
    const bool shallMove = controllerVelocity_.Length() > M_EPSILON;
    if (jumped_)
    {
        // Apply jump to velocity
        /*Vector3 linearVelocity =*/ body_->GetLinearVelocity();
        body_->SetLinearVelocity(currentJumpVelocity_);
    }
    else if (softGrounded_)
    {
        // Grounded movement, clip velocity if actually moving
        const Vector3 moveDirection = controllerVelocity_.Orthogonalize(contactNormal_);
        body_->ApplyImpulse(moveDirection * body_->GetMass() * groundAcceleration_);
        Vector3 linearVelocity = body_->GetLinearVelocity();
        if (shallMove)
            linearVelocity = linearVelocity.Normalized() * Min(linearVelocity.Length(), controllerVelocity_.Length());
        body_->SetLinearVelocity(linearVelocity);
    }
    else
    {
        // Fly movement, accelerate if target velocity isn't already reached
        const Vector3 accelDirection = ((controllerVelocity_ - body_->GetLinearVelocity()) * Vector3(1, 0, 1)).Normalized();
        body_->ApplyForce(body_->GetMass() * flyAcceleration_ * accelDirection);
    }

    // Update friction
    const bool shallStand = !shallMove && softGrounded_;
    body_->SetFriction(shallStand ? staticFriction_ : dynamicFriction_);

    contactNormal_ = Vector3::ZERO;
}

void CharacterController::FixedPostUpdate(float timeStep)
{
    // Re-apply jump to ensure correct velocity
    if (jumped_)
    {
        /*Vector3 linearVelocity =*/ body_->GetLinearVelocity();
        body_->SetLinearVelocity(currentJumpVelocity_);
    }
}

void CharacterController::OnNodeSet(Node* node)
{
    if (node)
        SubscribeToEvent(node, E_NODECOLLISION, URHO3D_HANDLER(CharacterController, HandleCollsion));
    else
        UnsubscribeFromEvent(E_NODECOLLISION);
}

void CharacterController::HandleCollsion(StringHash eventType, VariantMap& eventData)
{
    MemoryBuffer contacts(eventData["Contacts"].GetBuffer());

    while (!contacts.IsEof())
    {
        /*const Vector3 position =*/ contacts.ReadVector3();
        const Vector3 normal = contacts.ReadVector3();
        /*const float distance =*/ contacts.ReadFloat();
        /*const float impulse =*/ contacts.ReadFloat();

        if (normal.y_ > 0 && (!strictSlopeLimit_ || normal.y_ > Cos(maxSlope_)))
            contactNormal_ += normal;
    }
}

bool CharacterController::CheckIntegrity()
{
    if (!node_)
        return false;
    if (!body_)
        body_ = node_->GetComponent<RigidBody>();
    return !!body_;
}

void RegisterCharacterController(Context* context)
{
    CharacterController::RegisterObject(context);
}

void RegisterCharacterControllerScriptAPI(asIScriptEngine* engine)
{
    RegisterComponent<CharacterController>(engine, "CharacterController");
    engine->RegisterObjectMethod("CharacterController", "void set_node(Node@+)", asMETHOD(CharacterController, SetNode), asCALL_THISCALL);

    engine->RegisterObjectMethod("CharacterController", "void set_jumpVelocity(float)", asMETHOD(CharacterController, SetJumpVelocity), asCALL_THISCALL);
    engine->RegisterObjectMethod("CharacterController", "float get_jumpVelocity() const", asMETHOD(CharacterController, GetJumpVelocity), asCALL_THISCALL);
    engine->RegisterObjectMethod("CharacterController", "void set_flyAcceleration(float)", asMETHOD(CharacterController, SetFlyAcceleration), asCALL_THISCALL);
    engine->RegisterObjectMethod("CharacterController", "float get_flyAcceleration() const", asMETHOD(CharacterController, GetFlyAcceleration), asCALL_THISCALL);

    engine->RegisterObjectMethod("CharacterController", "void set_jump(bool)", asMETHOD(CharacterController, SetJump), asCALL_THISCALL);
    engine->RegisterObjectMethod("CharacterController", "void set_velocity(const Vector3&in)", asMETHOD(CharacterController, SetVelocity), asCALL_THISCALL);

    engine->RegisterObjectMethod("CharacterController", "bool get_hasJumped() const", asMETHOD(CharacterController, HasJumped), asCALL_THISCALL);
    engine->RegisterObjectMethod("CharacterController", "bool get_grounded() const", asMETHOD(CharacterController, IsGrounded), asCALL_THISCALL);
    engine->RegisterObjectMethod("CharacterController", "Vector3& get_velocity() const", asMETHOD(CharacterController, GetVelocity), asCALL_THISCALL);

    engine->RegisterObjectMethod("CharacterController", "void FixedUpdate(float)", asMETHOD(CharacterController, FixedUpdate), asCALL_THISCALL);
    engine->RegisterObjectMethod("CharacterController", "void FixedPostUpdate(float)", asMETHOD(CharacterController, FixedPostUpdate), asCALL_THISCALL);
}

}

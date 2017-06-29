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
    if (jumped_)
    {
        Vector3 linearVelocity = body_->GetLinearVelocity();
        linearVelocity.y_ = jumpVelocity_;
        body_->SetLinearVelocity(linearVelocity);
    }
    else if (softGrounded_)
    {
        const Vector3 accelerationDirection = controllerVelocity_.Orthogonalize(contactNormal_);
        body_->ApplyImpulse(accelerationDirection * body_->GetMass() * acceleration_);
        Vector3 linearVelocity = body_->GetLinearVelocity();
        linearVelocity = linearVelocity.Normalized() * Min(linearVelocity.Length(), controllerVelocity_.Length());
        body_->SetLinearVelocity(linearVelocity);
    }

    // Update friction
    const bool wantToMove = controllerVelocity_.Length() > M_EPSILON;
    body_->SetFriction(!wantToMove && softGrounded_ ? staticFriction_ : dynamicFriction_);

    contactNormal_ = Vector3::ZERO;

    linearVelocity_ = body_->GetLinearVelocity();
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

    engine->RegisterObjectMethod("CharacterController", "void set_jump(bool)", asMETHOD(CharacterController, SetJump), asCALL_THISCALL);
    engine->RegisterObjectMethod("CharacterController", "void set_velocity(const Vector3&in)", asMETHOD(CharacterController, SetVelocity), asCALL_THISCALL);

    engine->RegisterObjectMethod("CharacterController", "bool get_hasJumped() const", asMETHOD(CharacterController, HasJumped), asCALL_THISCALL);
    engine->RegisterObjectMethod("CharacterController", "bool get_grounded() const", asMETHOD(CharacterController, IsGrounded), asCALL_THISCALL);
    engine->RegisterObjectMethod("CharacterController", "Vector3& get_velocity() const", asMETHOD(CharacterController, GetVelocity), asCALL_THISCALL);

    engine->RegisterObjectMethod("CharacterController", "void FixedUpdate(float)", asMETHOD(CharacterController, FixedUpdate), asCALL_THISCALL);
}

}

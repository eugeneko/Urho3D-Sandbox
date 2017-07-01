#include "Urho3DPlayer.h"
#include <CharacterAnimator/CharacterAnimator.h>
#include <CharacterController/CharacterController.h>

using namespace Urho3D;

/// FlexEnginePlayer application runs a script specified on the command line.
class FlexEnginePlayer : public Urho3DPlayer
{
    URHO3D_OBJECT(FlexEnginePlayer, Urho3DPlayer);

public:
    /// Construct.
    FlexEnginePlayer(Context* context);

    /// Setup after engine initialization. Load the script and execute its start function.
    virtual void Start();
};

#include <Urho3D/AngelScript/Script.h>
#include <Urho3D/Graphics/Renderer.h>
// #include <Urho3D/IO/FileSystem.h>
#include <Urho3D/Resource/ResourceCache.h>
//
// #include <Urho3D/DebugNew.h>

URHO3D_DEFINE_APPLICATION_MAIN(FlexEnginePlayer);

FlexEnginePlayer::FlexEnginePlayer(Context* context) :
    Urho3DPlayer(context)
{
}

void FlexEnginePlayer::Start()
{
    ResourceCache* resourceCache = GetSubsystem<ResourceCache>();
    GetSubsystem<Renderer>()->SetMinInstances(1);
    GetSubsystem<Renderer>()->SetNumExtraInstancingBufferElements(1);

    RegisterCharacterAnimator(context_);
    RegisterCharacterController(context_);

    Urho3DPlayer::Start();

    Script* scriptSubsystem = GetSubsystem<Script>();
    asIScriptEngine* scriptEngine = scriptSubsystem->GetScriptEngine();
    RegisterCharacterAnimatorScriptAPI(scriptEngine);
    RegisterCharacterControllerScriptAPI(scriptEngine);
}

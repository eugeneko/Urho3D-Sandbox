#include <Urho3D/Core/Main.h>
#include <Urho3D/Engine/Application.h>
#include <Urho3D/Engine/EngineDefs.h>
#include <Urho3D/Graphics/Model.h>
#include <Urho3D/Graphics/Animation.h>
#include <Urho3D/IO/File.h>
#include <Urho3D/IO/Log.h>
#include <Urho3D/Input/Input.h>
#include <Urho3D/UI/UI.h>
#include <CharacterAnimator/CharacterAnimator.h>

using namespace Urho3D;

template <class T>
SharedPtr<T> LoadResource(Context* context, const String& fileName)
{
    File file(context, fileName);
    SharedPtr<T> resource(new T(context));
    if (resource->Load(file))
        return resource;
    return nullptr;
}

int Main()
{
    // Check number of arguments
    const StringVector& arguments = GetArguments();
    if (arguments.Size() < 4)
    {
        Log::WriteRaw("Usage: CharacterAnimationTool model.mdl skeleton.xml animation.ani output.xml");
        return 1;
    }

    // Initialize engine
    SharedPtr<Context> context(new Context);
    SharedPtr<Engine> engine(new Engine(context));
    VariantMap parameters;
    parameters[EP_HEADLESS] = true;
    parameters[EP_RESOURCE_PATHS] = ".";
    parameters[EP_LOG_NAME] = "";
    if (!engine->Initialize(parameters))
    {
        Log::WriteRaw("Failed to initialize engine");
        return 2;
    }

    const String& modelFileName = arguments[0];
    const String& skeletonFileName = arguments[1];
    const String& animationFileName = arguments[2];
    const String& outputFileName = arguments[3];

    SharedPtr<Model> model(LoadResource<Model>(context, modelFileName));
    SharedPtr<Animation> animation(LoadResource<Animation>(context, animationFileName));
    SharedPtr<CharacterSkeleton> characterSkeleton(LoadResource<CharacterSkeleton>(context, skeletonFileName));
    if (!model || !animation || !characterSkeleton)
    {
        Log::WriteRaw("Failed to open input resources");
        return 3;
    }

    CharacterAnimation characterAnimation(context);
    characterAnimation.SetName(outputFileName);
    if (!characterAnimation.Import(*animation, *model, *characterSkeleton))
    {
        Log::WriteRaw("Failed to import animation");
        return 4;
    }

    if (!characterAnimation.SaveFile(outputFileName))
    {
        Log::WriteRaw("Failed to save destination file");
        return 5;
    }
    return 0;
}

URHO3D_DEFINE_MAIN(Main());

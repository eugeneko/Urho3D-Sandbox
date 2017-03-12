#include <Urho3D/Core/Main.h>
#include <Urho3D/Engine/Application.h>
#include <Urho3D/Engine/EngineDefs.h>
#include <Urho3D/Graphics/Model.h>
#include <Urho3D/Graphics/Animation.h>
#include <Urho3D/IO/File.h>
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
        return 2;

    // Initialize engine
    SharedPtr<Context> context(new Context);
    SharedPtr<Engine> engine(new Engine(context));

    const String& modelFileName = arguments[0];
    const String& animationFileName = arguments[1];
    const String& skeletonFileName = arguments[2];
    const String& outputFileName = arguments[3];

    Model* model = LoadResource<Model>(context, modelFileName);
    Animation* animation = LoadResource<Animation>(context, animationFileName);
    CharacterSkeleton* characterSkeleton = LoadResource<CharacterSkeleton>(context, skeletonFileName);
    if (!model || !animation || !characterSkeleton)
        return 3;

    CharacterAnimation characterAnimation(context);
    characterAnimation.SetName(outputFileName);
    if (characterAnimation.ImportAnimation(*characterSkeleton, *model, *animation))
        characterAnimation.SaveFile(outputFileName);
    else
        return 4;
    
    return 0;
}

URHO3D_DEFINE_MAIN(Main());

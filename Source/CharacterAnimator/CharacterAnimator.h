#pragma once

#include <Urho3D/Graphics/AnimationController.h>
#include <Urho3D/Resource/Resource.h>
#include <Urho3D/Resource/XMLFile.h>
#include <Urho3D/Scene/LogicComponent.h>

class asIScriptEngine;

namespace Urho3D
{

class Animation;
class AnimatedModel;
class Model;

}

namespace Urho3D
{

/// Description of character skeleton 2-segment.
struct CharacterSkeletonSegment2
{
    /// Name of the segment.
    String name_;
    /// Name of the bone that is related to the root position of the segment.
    String rootBone_;
    /// Name of the bone that is related to the joint position of the segment. Must be a child of the root bone.
    String jointBone_;
    /// Name of the bone that is related to the target position of the segment. Must be a child of the joint bone.
    String targetBone_;
};

/// Character skeleton.
class CharacterSkeleton : public Resource
{
    URHO3D_OBJECT(CharacterSkeleton, Resource);

public:
    /// Container of 2-segment.
    using Segment2Map = HashMap<StringHash, CharacterSkeletonSegment2>;

public:
    /// Construct.
    CharacterSkeleton(Context* context);
    /// Destruct.
    virtual ~CharacterSkeleton();
    /// Register object factory.
    static void RegisterObject(Context* context);

    /// Load resource from stream. May be called from a worker thread. Return true if successful.
    virtual bool BeginLoad(Deserializer& source);

    /// Load from an XML element. Return true if successful.
    bool BeginLoad(const XMLElement& source);

    /// Get 2-segments.
    const Segment2Map& GetSegments2() const { return segments2_; }

private:
    /// Model name.
    String modelName_;
    /// Container of 2-segments.
    Segment2Map segments2_;
};

/// Blend animations and return result.
SharedPtr<Animation> BlendAnimations(Model& model, CharacterSkeleton* skeleton,
    const PODVector<Animation*>& animations,
    const PODVector<float>& weights, const PODVector<float>& offsets, const PODVector<float>& timestamps);

/// Key frame of 2-segment animation track.
// #TODO Rename members!
struct CharacterAnimationSegment2KeyFrame
{
    /// Key frame time.
    float time_ = 0.0f;
    /// Heel position.
    Vector3 heelPosition_;
    /// Direction of knee.
    Vector3 kneeDirection_;
    /// Fix for thigh rotation.
    Quaternion thighRotationFix_;
    /// Fix for calf rotation.
    Quaternion calfRotationFix_;
    /// Heel rotation in local space.
    Quaternion heelRotationLocal_;
    /// Heel rotation in world space.
    Quaternion heelRotationWorld_;
};

/// Animation track of single 2-segment.
// #TODO Rename members!
struct CharacterAnimationSegment2Track
{
    /// Name of the segment.
    String name_;
    /// Base direction is used for resolving joint angle.
    Vector3 initialDirection_;
    /// Key frames.
    Vector<CharacterAnimationSegment2KeyFrame> keyFrames_;

    /// Get length of track.
    float GetLength() const { return keyFrames_.Empty() ? 0 : keyFrames_.Back().time_ - keyFrames_.Front().time_; }
    /// Get index of key frame.
    void GetKeyFrameIndex(float time, unsigned& index) const
    {
        if (time < 0.0f)
            time = 0.0f;

        if (index >= keyFrames_.Size())
            index = keyFrames_.Size() - 1;

        // Check for being too far ahead
        while (index && time < keyFrames_[index].time_)
            --index;

        // Check for being too far behind
        while (index < keyFrames_.Size() - 1 && time >= keyFrames_[index + 1].time_)
            ++index;
    }
    /// Sample frame at specified time.
    CharacterAnimationSegment2KeyFrame SampleFrame(float time, unsigned& frame) const
    {
        GetKeyFrameIndex(time, frame);
        const unsigned nextFrame = frame + 1 < keyFrames_.Size() ? frame + 1 : 0;
        const CharacterAnimationSegment2KeyFrame* keyFrame = &keyFrames_[frame];
        const CharacterAnimationSegment2KeyFrame* nextKeyFrame = &keyFrames_[nextFrame];
        float timeInterval = nextKeyFrame->time_ - keyFrame->time_;
        if (timeInterval < 0.0f)
            timeInterval += GetLength();

        float t = timeInterval > 0.0f ? (time - keyFrame->time_) / timeInterval : 1.0f;

        CharacterAnimationSegment2KeyFrame result;
        result.time_ = time;
        result.heelPosition_ = Lerp(keyFrame->heelPosition_, nextKeyFrame->heelPosition_, t);
        result.kneeDirection_ = Lerp(keyFrame->kneeDirection_, nextKeyFrame->kneeDirection_, t);
        result.thighRotationFix_ = keyFrame->thighRotationFix_.Slerp(nextKeyFrame->thighRotationFix_, t);
        result.calfRotationFix_ = keyFrame->calfRotationFix_.Slerp(nextKeyFrame->calfRotationFix_, t);
        result.heelRotationLocal_ = keyFrame->heelRotationLocal_.Slerp(nextKeyFrame->heelRotationLocal_, t);
        result.heelRotationWorld_ = keyFrame->heelRotationWorld_.Slerp(nextKeyFrame->heelRotationWorld_, t);
        return result;
    }
};

/// Character Animation.
class CharacterAnimation : public Resource
{
    URHO3D_OBJECT(CharacterAnimation, Resource);

public:
    /// Map of tracks for 2-segments.
    using Segment2TrackMap = HashMap<StringHash, CharacterAnimationSegment2Track>;
public:
    /// Construct.
    CharacterAnimation(Context* context);
    /// Destruct.
    virtual ~CharacterAnimation();
    /// Register object factory.
    static void RegisterObject(Context* context);

    /// Load resource from stream. May be called from a worker thread. Return true if successful.
    virtual bool BeginLoad(Deserializer& source);
    /// Save resource. Return true if successful.
    virtual bool Save(Serializer& dest) const;

    /// Load from an XML element. Return true if successful.
    bool Load(const XMLElement& source);
    /// Save to an XML element. Return true if successful.
    bool Save(XMLElement& dest) const;

    /// Import animation using model and skeleton.
    bool ImportAnimation(CharacterSkeleton& characterSkeleton, Model& model, Animation& animation);
    /// Find track by name.
    CharacterAnimationSegment2Track* FindTrack(const String& name) const { return segments2_[name]; }

private:
    /// Tracks for 2-segments.
    Segment2TrackMap segments2_;

};

/// Character Animation Controller.
class CharacterAnimationController : public AnimationController
{
    URHO3D_OBJECT(CharacterAnimationController, AnimationController);

public:
    /// Construct.
    CharacterAnimationController(Context* context);
    /// Destruct.
    virtual ~CharacterAnimationController();
    /// Register object factory.
    static void RegisterObject(Context* context);

    /// Set target transform of segment.
    void SetTargetTransform(StringHash segment, const Matrix3x4& transform);
    /// Set amount of transformation applied to rotation of target bone.
    void SetTargetRotationAmount(StringHash segment, float rotationAmount);
    /// Set target bone local/global rotation balance.
    void SetTargetRotationBalance(StringHash segment, float globalFactor);
    /// Clean up segment configuration.
    void CleanSegment2(StringHash segment);

    /// Update the animations. Is called from HandleScenePostUpdate().
    virtual void Update(float timeStep) override;
    /// Apply animations.
    void ApplyAnimation();

    /// Set skeleton attribute.
    void SetSkeletonAttr(const ResourceRef& value);
    /// Return skeleton attribute.
    ResourceRef GetSkeletonAttr() const;

private:
    /// Get character animation.
    CharacterAnimation* GetCharacterAnimation(const String& animationName);
    /// Apply animation.
    void ApplyAnimation(AnimatedModel* animatedModel);
    /// Update 2-segment.
    void UpdateSegment2(AnimatedModel* animatedModel, const CharacterSkeletonSegment2& segment);
    /// State of 2-segment.
    struct Segment2State
    {
        /// Transform of target point.
        Matrix3x4 targetTransform_;
        /// Adjusts amount of transformation applied to rotation of target bone.
        float targetRotationAmount_ = 0.0f;
        /// Adjusts balance between local and global rotation of target bone.
        float globalRotationFactor_ = 0.0f;
    };

private:
    /// Skeleton.
    SharedPtr<CharacterSkeleton> skeleton_;
    /// Cached animations.
    HashMap<StringHash, SharedPtr<CharacterAnimation>> animationCache_;
    /// States.
    HashMap<StringHash, Segment2State> segment2states_;

};

/// Register script API
void RegisterCharacterAnimatorScriptAPI(asIScriptEngine* engine);

}

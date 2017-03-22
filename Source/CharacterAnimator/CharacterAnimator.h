#pragma once

#include <Urho3D/Container/ArrayPtr.h>
#include <Urho3D/Graphics/AnimationController.h>
#include <Urho3D/Resource/Resource.h>
#include <Urho3D/Resource/XMLFile.h>
#include <Urho3D/Scene/LogicComponent.h>
#include <memory>

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

/// Type of character skeleton segment
enum class CharacterSkeletonSegmentType
{
    /// Track of the root node. Relative position and rotation are stored.
    Root,
    /// Track of the limb. Target position, 1D-rotations and joint direction are stored.
    Limb,
    /// Track of the node chain. Array of rotations and target position are stored.
    Chain
};

/// Get character skeleton segment type by name.
CharacterSkeletonSegmentType GetCharacterSkeletonSegmentType(const String& name);

struct CharacterSkeletonSegment;

/// Character skeleton segment data.
class CharacterSkeletonSegmentData
{
public:
    /// Create by type.
    static CharacterSkeletonSegmentData* Create(CharacterSkeletonSegmentType type);
    /// Destruct.
    virtual ~CharacterSkeletonSegmentData() {}
    /// Reset state.
    virtual void Reset();
    /// Merge state with specified weight to this.
    virtual void Merge(const CharacterSkeletonSegmentData& other, float weight) = 0;
    /// Apply segment animations to segment.
    virtual void Apply(const Matrix3x4& rootTransform, CharacterSkeletonSegment& dest) = 0;

protected:
    /// Weight.
    float accumulatedWeight_ = 0.0f;
};

/// Character skeleton segment.
struct CharacterSkeletonSegment
{
    /// Segment name.
    String name_;
    /// Segment type.
    CharacterSkeletonSegmentType type_;
    /// Names of bones in segment.
    Vector<String> boneNames_;

    /// Bones of the segment.
    PODVector<Bone*> bones_;
    /// Nodes of the segment.
    PODVector<Node*> nodes_;
    /// Global initial positions.
    PODVector<Vector3> globalPositions_;
    /// Global initial rotations.
    PODVector<Quaternion> globalRotations_;

    /// Data.
    std::shared_ptr<CharacterSkeletonSegmentData> data_;
};

/// Character skeleton root segment data.
class CharacterSkeletonRootSegmentData : public CharacterSkeletonSegmentData
{
public:
    /// Position.
    Vector3 position_;
    /// Rotation.
    Quaternion rotation_;

public:
    /// @see CharacterSkeletonSegmentData::Reset
    virtual void Reset() override;
    /// @see CharacterSkeletonSegmentData::Merge
    virtual void Merge(const CharacterSkeletonSegmentData& other, float weight) override;
    /// @see CharacterSkeletonSegmentData::Apply
    virtual void Apply(const Matrix3x4& rootTransform, CharacterSkeletonSegment& dest) override;

};

/// Character skeleton chain segment data.
class CharacterSkeletonChainSegmentData : public CharacterSkeletonSegmentData
{
public:
    /// Target position.
    Vector3 position_;
    /// Segment rotations.
    PODVector<Quaternion> rotations_;

public:
    /// @see CharacterSkeletonSegmentData::Reset
    virtual void Reset() override;
    /// @see CharacterSkeletonSegmentData::Merge
    virtual void Merge(const CharacterSkeletonSegmentData& other, float weight) override;
    /// @see CharacterSkeletonSegmentData::Apply
    virtual void Apply(const Matrix3x4& rootTransform, CharacterSkeletonSegment& dest) override;
};

/// Character skeleton limb segment data.
class CharacterSkeletonLimbSegmentData : public CharacterSkeletonSegmentData
{
public:
    /// Target position.
    Vector3 position_;
    /// Limb direction.
    Vector3 direction_;
    /// First segment rotation.
    float rotation0_ = 0.0f;
    /// Second segment rotation.
    float rotation1_ = 0.0f;
    /// Target segment rotation.
    Quaternion rotation2_;

public:
    /// @see CharacterSkeletonSegmentData::Reset
    virtual void Reset() override;
    /// @see CharacterSkeletonSegmentData::Merge
    virtual void Merge(const CharacterSkeletonSegmentData& other, float weight) override;
    /// @see CharacterSkeletonSegmentData::Apply
    virtual void Apply(const Matrix3x4& rootTransform, CharacterSkeletonSegment& dest) override;
};

/// Character skeleton.
class CharacterSkeleton : public Resource
{
    URHO3D_OBJECT(CharacterSkeleton, Resource);

public:
    /// Construct.
    CharacterSkeleton(Context* context) : Resource(context) {}
    /// Destruct.
    virtual ~CharacterSkeleton() {}
    /// Register object factory.
    static void RegisterObject(Context* context);
    /// Load resource from stream. May be called from a worker thread. Return true if successful.
    virtual bool BeginLoad(Deserializer& source);
    /// Load from an XML element. Return true if successful.
    bool LoadXML(const XMLElement& source);

    /// Get segments.
    const Vector<CharacterSkeletonSegment>& GetSegments() const { return segments_; }
    /// Find segment by name.
    const CharacterSkeletonSegment* FindSegment(const String& name) const;

    /// Allocate segment data.
    bool AllocateSegmentData(Vector<CharacterSkeletonSegment>& segmentsData, Skeleton& skeleton, const Matrix3x4& baseTransform);

private:
    /// Segments.
    Vector<CharacterSkeletonSegment> segments_;
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

/// Bone of character skeleton.
struct CharacterSkeletonBone
{
    /// Bone.
    Bone* bone_;
    /// Node.
    Node* node_;
    /// Global initial position.
    Vector3 globalInitialPosition_;
    /// Global initial rotation.
    Quaternion globalInitialRotation_;
};

/// Character animation track base class.
class CharacterAnimationTrack : public RefCounted
{
public:
    /// Create track by type.
    static SharedPtr<CharacterAnimationTrack> Create(CharacterSkeletonSegmentType type, const String& name);
    /// Construct.
    CharacterAnimationTrack(CharacterSkeletonSegmentType type, const String& name) : type_(type), name_(name) {}
    /// Destruct.
    virtual ~CharacterAnimationTrack() {}
    /// Import frame.
    virtual void ImportFrame(const CharacterSkeletonSegment& segment) = 0;
    /// Merge frame to destination with specified weight.
    virtual void MergeFrame(unsigned firstFrame, unsigned secondFrame, float factor,
        float weight, CharacterSkeletonSegmentData& dest) = 0;
    /// Save XML.
    virtual bool SaveXML(XMLElement& dest) const = 0;
    /// Load XML.
    virtual bool LoadXML(const XMLElement& source) = 0;
    /// Get type string.
    virtual String GetTypeString() const = 0;
    /// Check whether the number of bones is valid.
    virtual bool CheckNumberOfBones(unsigned numBones) const = 0;
    /// Get name.
    String GetName() const { return name_; }

private:
    /// Track type.
    CharacterSkeletonSegmentType type_;
    /// Track name.
    const String name_;
};

/// Root animation track.
class RootAnimationTrack : public CharacterAnimationTrack
{
public:
    /// Construct.
    RootAnimationTrack(const String& name) : CharacterAnimationTrack(CharacterSkeletonSegmentType::Root, name) {}
    /// @see CharacterAnimationTrack::ImportFrame
    virtual void ImportFrame(const CharacterSkeletonSegment& segment) override;
    /// @see CharacterAnimationTrack::MergeFrame
    virtual void MergeFrame(unsigned firstFrame, unsigned secondFrame, float factor,
        float weight, CharacterSkeletonSegmentData& dest) override;
    /// @see CharacterAnimationTrack::SaveXML
    virtual bool SaveXML(XMLElement& dest) const override;
    /// @see CharacterAnimationTrack::LoadXML
    virtual bool LoadXML(const XMLElement& source) override;
    /// @see CharacterAnimationTrack::GetTypeString
    virtual String GetTypeString() const override { return "root"; }
    /// @see CharacterAnimationTrack::CheckNumberOfBones
    virtual bool CheckNumberOfBones(unsigned numBones) const override { return numBones == 1; }

public:
    /// Track.
    Vector<CharacterSkeletonRootSegmentData> track_;
};

/// Chain animation track.
class ChainAnimationTrack : public CharacterAnimationTrack
{
public:
    /// Construct.
    ChainAnimationTrack(const String& name) : CharacterAnimationTrack(CharacterSkeletonSegmentType::Chain, name) {}
    /// @see CharacterAnimationTrack::ImportFrame
    virtual void ImportFrame(const CharacterSkeletonSegment& segment) override;
    /// @see CharacterAnimationTrack::MergeFrame
    virtual void MergeFrame(unsigned firstFrame, unsigned secondFrame, float factor,
        float weight, CharacterSkeletonSegmentData& dest) override;
    /// @see CharacterAnimationTrack::SaveXML
    virtual bool SaveXML(XMLElement& dest) const override;
    /// @see CharacterAnimationTrack::LoadXML
    virtual bool LoadXML(const XMLElement& source) override;
    /// @see CharacterAnimationTrack::GetTypeString
    virtual String GetTypeString() const override { return "chain"; }
    /// @see CharacterAnimationTrack::CheckNumberOfBones
    virtual bool CheckNumberOfBones(unsigned numBones) const override { return numBones > 2; }

public:
    /// Track.
    Vector<CharacterSkeletonChainSegmentData> track_;
};

/// Limb animation track.
class LimbAnimationTrack : public CharacterAnimationTrack
{
public:
    /// Construct.
    LimbAnimationTrack(const String& name) : CharacterAnimationTrack(CharacterSkeletonSegmentType::Limb, name) {}
    /// @see CharacterAnimationTrack::ImportFrame
    virtual void ImportFrame(const CharacterSkeletonSegment& segment) override;
    /// @see CharacterAnimationTrack::MergeFrame
    virtual void MergeFrame(unsigned firstFrame, unsigned secondFrame, float factor,
        float weight, CharacterSkeletonSegmentData& dest) override;
    /// @see CharacterAnimationTrack::SaveXML
    virtual bool SaveXML(XMLElement& dest) const override;
    /// @see CharacterAnimationTrack::LoadXML
    virtual bool LoadXML(const XMLElement& source) override;
    /// @see CharacterAnimationTrack::GetTypeString
    virtual String GetTypeString() const override { return "limb"; }
    /// @see CharacterAnimationTrack::CheckNumberOfBones
    virtual bool CheckNumberOfBones(unsigned numBones) const override { return numBones > 2; }

public:
    /// Track.
    Vector<CharacterSkeletonLimbSegmentData> track_;
};

/// Character Animation.
class CharacterAnimation : public Resource
{
    URHO3D_OBJECT(CharacterAnimation, Resource);

public:
    /// Construct.
    CharacterAnimation(Context* context) : Resource(context) {}
    /// Register object factory.
    static void RegisterObject(Context* context);
    /// Import joint animation.
    bool Import(Animation& animation, Model& model, CharacterSkeleton& rig, const Matrix3x4& transform = Matrix3x4::IDENTITY);

    /// Load resource from stream. May be called from a worker thread. Return true if successful.
    virtual bool BeginLoad(Deserializer& source);
    /// Save resource. Return true if successful.
    virtual bool Save(Serializer& dest) const;

    /// Load from an XML element. Return true if successful.
    bool LoadXML(const XMLElement& source);
    /// Save to an XML element. Return true if successful.
    bool SaveXML(XMLElement& dest) const;

    /// Get tracks.
    const Vector<SharedPtr<CharacterAnimationTrack>>& GetTracks() const { return tracks_; }
    /// Find track by name.
    CharacterAnimationTrack* FindTrack(const String& name) const;
    /// Get index of key frame.
    void GetKeyFrameIndex(float time, unsigned& index) const;
    /// Get index of key frame.
    unsigned GetKeyFrameIndex(float time) const;
    /// Get interpolated key frame indexes.
    void GetKeyFrame(float time, unsigned& firstFrame, unsigned& secondFrame, float& factor) const;
    /// Get animation length.
    float GetLength() const;

private:
    /// Time stamps.
    PODVector<float> timeStamps_;
    /// Tracks.
    Vector<SharedPtr<CharacterAnimationTrack>> tracks_;
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

    /// Set animation transform.
    void SetAnimationTransform(const Matrix3x4& transform);
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

    /// Set skeleton.
    void SetSkeleton(CharacterSkeleton* skeleton);
    /// Set skeleton attribute.
    void SetSkeletonAttr(const ResourceRef& value);
    /// Return skeleton attribute.
    ResourceRef GetSkeletonAttr() const;

private:
    /// Lazy initialize model and hierarchy if needed.
    void UpdateHierarchy();
    /// Get character animation.
    CharacterAnimation* GetCharacterAnimation(const String& animationName);

    /// Update 2-segment.
    void UpdateSegment2(const CharacterSkeletonSegment2& segment);
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
    /// Animated model.
    WeakPtr<AnimatedModel> animatedModel_;
    /// Animated model skeleton.
    Skeleton* animatedModelSkeleton_ = nullptr;
    /// Segment data.
    Vector<CharacterSkeletonSegment> segmentData_;

    /// Cached animations.
    HashMap<StringHash, SharedPtr<CharacterAnimation>> animationCache_;
    /// States.
    HashMap<StringHash, Segment2State> segment2states_;
    /// Animation transform.
    Matrix3x4 animationTransform_;

};

/// Register script API
void RegisterCharacterAnimatorScriptAPI(asIScriptEngine* engine);

}

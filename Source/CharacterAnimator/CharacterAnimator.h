#pragma once

#include <Urho3D/Container/ArrayPtr.h>
#include <Urho3D/Graphics/AnimationController.h>
#include <Urho3D/Resource/Resource.h>
#include <Urho3D/Resource/XMLFile.h>
#include <Urho3D/Scene/LogicComponent.h>
#include <Urho3D/Scene/Node.h>
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
    /// Initial transforms.
    PODVector<Matrix3x4> initialPose_;
    // #TODO Fix duplicate
    /// Global initial positions.
    PODVector<Vector3> globalPositions_;
    /// Global initial rotations.
    PODVector<Quaternion> globalRotations_;
};

/// Character segment data template.
template <CharacterSkeletonSegmentType T>
struct CharacterSegmentDataT
{
    /// Segment type.
    static const CharacterSkeletonSegmentType Type = T;
};

/// Character root segment data.
struct CharacterRootSegmentData : public CharacterSegmentDataT<CharacterSkeletonSegmentType::Root>
{
    /// Position.
    Vector3 position_;
    /// Rotation.
    Quaternion rotation_;

public:
    /// Reset state.
    void Reset();
    /// Blend state with another.
    void Blend(const CharacterRootSegmentData& other, float weight);
    /// Render state from animation space to world space.
    void Render(const Matrix3x4& rootTransform, const Matrix3x4& animTransform,
        const Vector3& basePosition, const Quaternion& baseRotation);
    /// Import state from nodes.
    void Import(const CharacterSkeletonSegment& src);
    /// Export state to nodes.
    void Export(const Matrix3x4& rootTransform, const Matrix3x4& animTransform, CharacterSkeletonSegment& dest) const;
    /// Save to XML.
    bool Save(XMLElement& dest) const;
    /// Load from XML.
    bool Load(const XMLElement& src);
};

/// Character chain segment data.
struct CharacterChainSegmentData : public CharacterSegmentDataT<CharacterSkeletonSegmentType::Chain>
{
    /// Target position.
    Vector3 position_;
    /// Segment rotations.
    PODVector<Quaternion> rotations_;

public:
    /// Reset state.
    void Reset();
    /// Blend state with another.
    void Blend(const CharacterChainSegmentData& other, float weight);
    /// Render state from animation space to world space.
    void Render(const Matrix3x4& rootTransform, const Matrix3x4& animTransform, const CharacterSkeletonSegment& segment);
    /// Import state from nodes.
    void Import(const CharacterSkeletonSegment& src);
    /// Export state to nodes.
    void Export(const Matrix3x4& rootTransform, const Matrix3x4& animTransform, CharacterSkeletonSegment& dest) const;
    /// Save to XML.
    bool Save(XMLElement& dest) const;
    /// Load from XML.
    bool Load(const XMLElement& src);
};

/// Character limb segment data.
struct CharacterLimbSegmentData : public CharacterSegmentDataT<CharacterSkeletonSegmentType::Limb>
{
    /// Total length of limb.
    float length_ = 0.0f;
    /// Target position.
    Vector3 position_;
    /// Limb direction.
    Vector3 direction_;
    /// First segment rotation.
    float rotationA_ = 0.0f;
    /// Second segment rotation.
    float rotationB_ = 0.0f;
    /// Target segment rotation.
    Quaternion rotationC_;

public:
    /// Reset state.
    void Reset();
    /// Blend state with another.
    void Blend(const CharacterLimbSegmentData& other, float weight);
    /// Render state from animation space to world space.
    void Render(const Matrix3x4& rootTransform, const Matrix3x4& animTransform, const CharacterSkeletonSegment& segment);
    /// Import state from nodes.
    void Import(const CharacterSkeletonSegment& src);
    /// Export state to nodes.
    void Export(const Matrix3x4& rootTransform, const Matrix3x4& animTransform, CharacterSkeletonSegment& dest) const;
    /// Save to XML.
    bool Save(XMLElement& dest) const;
    /// Load from XML.
    bool Load(const XMLElement& src);
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
class CharacterAnimationTrack : public Object
{
    URHO3D_OBJECT(CharacterAnimationTrack, Object);

protected:
    /// Construct.
    CharacterAnimationTrack(Context* context, CharacterSkeletonSegmentType type) : Object(context), type_(type) {}

public:
    /// Create track by type.
    static SharedPtr<CharacterAnimationTrack> Create(Context* context, CharacterSkeletonSegmentType type);
    /// Destruct.
    virtual ~CharacterAnimationTrack() {}

    /// Set name.
    void SetName(const String& name) { name_ = name; }
    /// Get name.
    String GetName() const { return name_; }

    /// Get index of key frame.
    void GetKeyFrameIndex(float time, unsigned& index) const;
    /// Get index of key frame.
    unsigned GetKeyFrameIndex(float time) const;
    /// Get interpolated key frame indexes.
    void GetKeyFrame(float time, unsigned& firstFrame, unsigned& secondFrame, float& factor) const;
    /// Get animation length.
    float GetLength() const;

    /// Import frame.
    virtual void ImportFrame(float time, const CharacterSkeletonSegment& segment) = 0;
    /// Save XML.
    virtual bool SaveXML(XMLElement& dest) const = 0;
    /// Load XML.
    virtual bool LoadXML(const XMLElement& source) = 0;

    /// Get type string.
    virtual String GetTypeString() const = 0;
    /// Check whether the number of bones is valid.
    virtual bool CheckNumberOfBones(unsigned numBones) const = 0;

private:
    /// Track type.
    CharacterSkeletonSegmentType type_;
    /// Track name.
    String name_;

protected:
    /// Time stamps.
    PODVector<float> timeStamps_;
};

/// Character animation track template.
template <class TAnimationFrame>
class CharacterAnimationTrackT : public CharacterAnimationTrack
{
public:
    /// Get animation frame by index.
    const TAnimationFrame& GetFrame(unsigned index) const { return track_[index]; }

    /// @see CharacterAnimationTrack::ImportFrame
    virtual void ImportFrame(float time, const CharacterSkeletonSegment& segment) override final
    {
        TAnimationFrame frame;
        frame.Import(segment);
        timeStamps_.Push(time);
        track_.Push(frame);
    }
    /// @see CharacterAnimationTrack::SaveXML
    virtual bool SaveXML(XMLElement& dest) const override final
    {
        assert(track_.Size() == timeStamps_.Size());
        const unsigned numKeys = track_.Size();
        for (unsigned i = 0; i < numKeys; ++i)
        {
            XMLElement child = dest.CreateChild("frame");
            child.SetFloat("time", timeStamps_[i]);
            if (!track_[i].Save(child))
                return false;
        }
        return true;
    }
    /// @see CharacterAnimationTrack::LoadXML
    virtual bool LoadXML(const XMLElement& source) override final
    {
        for (XMLElement child = source.GetChild("frame"); child; child = child.GetNext())
        {
            TAnimationFrame frame;
            const float time = child.GetFloat("time");
            if (!frame.Load(child))
                return false;
            timeStamps_.Push(time);
            track_.Push(frame);
        }
        return true;
    }

protected:
    /// Construct.
    CharacterAnimationTrackT(Context* context) : CharacterAnimationTrack(context, TAnimationFrame::Type) {}
    /// Track frames.
    Vector<TAnimationFrame> track_;
};

/// Root animation track.
class RootAnimationTrack : public CharacterAnimationTrackT<CharacterRootSegmentData>
{
public:
    /// Construct.
    RootAnimationTrack(Context* context) : CharacterAnimationTrackT<CharacterRootSegmentData>(context) {}
    /// @see CharacterAnimationTrack::GetTypeString
    virtual String GetTypeString() const override { return "root"; }
    /// @see CharacterAnimationTrack::CheckNumberOfBones
    virtual bool CheckNumberOfBones(unsigned numBones) const override { return numBones == 1; }
};

/// Chain animation track.
class ChainAnimationTrack : public CharacterAnimationTrackT<CharacterChainSegmentData>
{
public:
    /// Construct.
    ChainAnimationTrack(Context* context) : CharacterAnimationTrackT<CharacterChainSegmentData>(context) {}
    /// @see CharacterAnimationTrack::GetTypeString
    virtual String GetTypeString() const override { return "chain"; }
    /// @see CharacterAnimationTrack::CheckNumberOfBones
    virtual bool CheckNumberOfBones(unsigned numBones) const override { return numBones > 2; }
};

/// Limb animation track.
class LimbAnimationTrack : public CharacterAnimationTrackT<CharacterLimbSegmentData>
{
public:
    /// Construct.
    LimbAnimationTrack(Context* context) : CharacterAnimationTrackT<CharacterLimbSegmentData>(context) {}
    /// @see CharacterAnimationTrack::GetTypeString
    virtual String GetTypeString() const override { return "limb"; }
    /// @see CharacterAnimationTrack::CheckNumberOfBones
    virtual bool CheckNumberOfBones(unsigned numBones) const override { return numBones > 2; }
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

private:
    /// Tracks.
    Vector<SharedPtr<CharacterAnimationTrack>> tracks_;
};

/// Character Segment Effector.
class CharacterEffector : public Component
{
    URHO3D_OBJECT(CharacterEffector, Component);

public:
    /// Construct.
    CharacterEffector(Context* context) : Component(context) {}
    /// Destruct.
    virtual ~CharacterEffector() {}
    /// Register object factory.
    static void RegisterObject(Context* context);

    // /// Return whether the animation is played.
    // virtual bool IsAnimated() const { return true; }
    /// Reset transforms.
    void ResetTransforms(CharacterSkeletonSegment& segment);
    /// Reset animation state.
    virtual void ResetAnimationState() = 0;
    /// Apply animation track.
    virtual void ApplyAnimationTrack(float weight, float time, CharacterAnimationTrack& animationTrack) = 0;
    /// Resolve animations.
    virtual void ResolveAnimations(Node& rootNode, const Matrix3x4& animTransform, CharacterSkeletonSegment& dest) = 0;

    /// Get segment name.
    const String& GetSegmentName() const { return segmentName_; }

private:
    /// @see Component::OnNodeSet
    virtual void OnNodeSet(Node* node) override;

private:
    /// Segment name.
    String segmentName_;
};

/// Character Effector template.
template <class TAnimationTrack, class TAnimationFrame>
class CharacterEffectorT : public CharacterEffector
{
public:
    /// @see CharacterEffector::ResetAnimationState
    virtual void ResetAnimationState() override final
    {
        accumulatedWeight_ = 0.0f;
        animationState_.Reset();
    }
    /// @see CharacterEffector::ApplyAnimationTrack
    virtual void ApplyAnimationTrack(float weight, float time, CharacterAnimationTrack& animationTrack) override final
    {
        if (TAnimationTrack* track = dynamic_cast<TAnimationTrack*>(&animationTrack))
        {
            unsigned firstFrame;
            unsigned secondFrame;
            float factor;
            track->GetKeyFrame(time, firstFrame, secondFrame, factor);

            const float firstWeight = weight * (1 - factor);
            const float secondWeight = weight * factor;

            animationState_.Blend(track->GetFrame(firstFrame), firstWeight / (firstWeight + accumulatedWeight_));
            accumulatedWeight_ += firstWeight;
            animationState_.Blend(track->GetFrame(secondFrame), secondWeight / (secondWeight + accumulatedWeight_));
            accumulatedWeight_ += secondWeight;
        }
    }
    /// @see CharacterEffector::ResolveAnimations
    virtual void ResolveAnimations(Node& rootNode, const Matrix3x4& animTransform, CharacterSkeletonSegment& dest) override
    {
        ResolveAnimationState(effectorState_, animationState_, rootNode, animTransform, dest);
        effectorState_.Export(rootNode.GetWorldTransform(), animTransform, dest);
    }

protected:
    /// Construct.
    CharacterEffectorT(Context* context) : CharacterEffector(context) {}
    /// Resolve animation state.
    virtual void ResolveAnimationState(TAnimationFrame& effectorState, TAnimationFrame& animationState,
        Node& rootNode, const Matrix3x4& animTransform, const CharacterSkeletonSegment& segment) = 0;

private:
    /// Accumulated weight.
    float accumulatedWeight_ = 0.0f;
    /// Current animation state.
    TAnimationFrame animationState_;
    /// Effector state.
    TAnimationFrame effectorState_;
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

    /// Add segment effector.
    void AddEffector(CharacterEffector* effector);
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
    /// Get segment by name.
    CharacterSkeletonSegment* GetSegment(const String& segmentName);
    /// Get effector by name.
    CharacterEffector* GetEffector(const String& segmentName);

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
    /// Cached animation data.
    struct CachedAnimationData
    {
        /// Character animation.
        CharacterAnimation* characterAnimation_ = nullptr;
        /// Animation time.
        float time_ = 0.0f;
        /// Animation weight.
        float weight_ = 0.0f;
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
    /// Current animation data.
    Vector<CachedAnimationData> currentAnimationData_;
    /// Current segments data.
    Vector<Pair<CharacterEffector*, CharacterSkeletonSegment*>> currentSegmentsData_;
    /// States.
    HashMap<StringHash, Segment2State> segment2states_;
    /// Animation transform.
    Matrix3x4 animationTransform_;
    /// Segment controllers.
    Vector<WeakPtr<CharacterEffector>> effectors_;
};

/// Character Root Effector.
class CharacterRootEffector : public CharacterEffectorT<RootAnimationTrack, CharacterRootSegmentData>
{
    URHO3D_OBJECT(CharacterRootEffector, CharacterEffector);

public:
    /// Construct.
    CharacterRootEffector(Context* context) : CharacterEffectorT<RootAnimationTrack, CharacterRootSegmentData>(context) {}
    /// Destruct.
    virtual ~CharacterRootEffector() {}
    /// Register object factory.
    static void RegisterObject(Context* context);

private:
    /// @see CharacterEffectorT::ResolveAnimationState
    virtual void ResolveAnimationState(CharacterRootSegmentData& effectorState, CharacterRootSegmentData& animationState,
        Node& rootNode, const Matrix3x4& animTransform, const CharacterSkeletonSegment& segment) override;

private:
    /// Multiplier of root offset during animation.
    float offsetFactor_ = 1.0f;
    /// Multiplier of root rotation during animation.
    float rotationFactor_ = 1.0f;
};

/// Character Limb Effector.
class CharacterLimbEffector : public CharacterEffectorT<LimbAnimationTrack, CharacterLimbSegmentData>
{
    URHO3D_OBJECT(CharacterLimbEffector, CharacterEffector);

public:
    /// Construct.
    CharacterLimbEffector(Context* context) : CharacterEffectorT<LimbAnimationTrack, CharacterLimbSegmentData>(context) {}
    /// Destruct.
    virtual ~CharacterLimbEffector() {}
    /// Register object factory.
    static void RegisterObject(Context* context);

private:
    /// @see CharacterEffectorT::ResolveAnimationState
    virtual void ResolveAnimationState(CharacterLimbSegmentData& effectorState, CharacterLimbSegmentData& animationState,
        Node& rootNode, const Matrix3x4& animTransform, const CharacterSkeletonSegment& segment) override;

private:
    /// Whether to animate limb target position.
    bool animateTargetPosition_ = false;
    /// Whether to animate limb target rotation.
    bool animateTargetRotation_ = false;
};

/// Character Chain Effector.
class CharacterChainEffector : public CharacterEffectorT<ChainAnimationTrack, CharacterChainSegmentData>
{
    URHO3D_OBJECT(CharacterChainEffector, CharacterEffector);

public:
    /// Construct.
    CharacterChainEffector(Context* context) : CharacterEffectorT<ChainAnimationTrack, CharacterChainSegmentData>(context) {}
    /// Destruct.
    virtual ~CharacterChainEffector() {}
    /// Register object factory.
    static void RegisterObject(Context* context);

private:
    /// @see CharacterEffectorT::ResolveAnimationState
    virtual void ResolveAnimationState(CharacterChainSegmentData& effectorState, CharacterChainSegmentData& animationState,
        Node& rootNode, const Matrix3x4& animTransform, const CharacterSkeletonSegment& segment)
    {
        animationState.Render(rootNode.GetWorldTransform(), animTransform, segment);
        effectorState = animationState;
    }
};

/// Register classes.
void RegisterCharacterAnimator(Context* context);
/// Register script API.
void RegisterCharacterAnimatorScriptAPI(asIScriptEngine* engine);

}

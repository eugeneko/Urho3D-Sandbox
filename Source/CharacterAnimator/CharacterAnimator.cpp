#include "CharacterAnimator.h"

#include <Urho3D/AngelScript/APITemplates.h>
#include <Urho3D/Container/Ptr.h>
#include <Urho3D/Core/Context.h>
#include <Urho3D/Graphics/AnimatedModel.h>
#include <Urho3D/Graphics/Animation.h>
#include <Urho3D/Graphics/AnimationController.h>
#include <Urho3D/Graphics/AnimationState.h>
#include <Urho3D/Graphics/DebugRenderer.h>
#include <Urho3D/IO/FileSystem.h>
#include <Urho3D/IO/Log.h>
#include <Urho3D/Math/MathDefs.h>
#include <Urho3D/Math/Sphere.h>
#include <Urho3D/Resource/ResourceCache.h>
#include <Urho3D/Scene/Serializable.h>

#include <algorithm>

namespace Urho3D
{

static const char* characterAnimatorCategory = "CharacterAnimator";

SharedPtr<Animation> BlendAnimations(Model& model, CharacterSkeleton* skeleton,
    const PODVector<Animation*>& animations,
    const PODVector<float>& weights, const PODVector<float>& offsets, const PODVector<float>& timestamps)
{
    // Create and setup node
    Node node(model.GetContext());
    AnimatedModel* animatedModel = node.CreateComponent<AnimatedModel>();
    animatedModel->SetModel(&model);

    // Create controller for better blending if skeleton is passed
    CharacterAnimationController* animationController = node.CreateComponent<CharacterAnimationController>();
    animationController->SetSkeletonAttr(ResourceRef(XMLFile::GetTypeStatic(), skeleton->GetName()));
    for (unsigned i = 0; i < animations.Size(); ++i)
    {
        const String& animationName = animations[i]->GetName();
        animationController->Play(animationName, 0, false);
        animationController->SetWeight(animationName, i < weights.Size() ? weights[i] : 1.0f);
    }

    // Get all nodes
    SharedPtr<Animation> result = MakeShared<Animation>(model.GetContext());
    Node* rootNode = animatedModel->GetSkeleton().GetRootBone()->node_;
    if (!rootNode)
        return result;

    PODVector<Node*> nodes;
    rootNode->GetChildren(nodes, true);
    nodes.Push(rootNode);

    // Create tracks
    for (Node* node : nodes)
        if (!node->GetName().Empty())
        {
            AnimationTrack* track = result->CreateTrack(node->GetName());
            track->channelMask_ = CHANNEL_POSITION | CHANNEL_ROTATION | CHANNEL_SCALE;
        }

    // Play animation
    const Vector<AnimationControl>& animationControls = animationController->GetAnimations();
    for (unsigned i = 0; i < timestamps.Size(); ++i)
    {
        const float time = timestamps[i];

        // Reset nodes
        for (Node* node : nodes)
            if (Bone* bone = model.GetSkeleton().GetBone(node->GetName()))
                node->SetTransform(bone->initialPosition_, bone->initialRotation_, bone->initialScale_);

        // Play animation
        for (unsigned j = 0; j < animationControls.Size(); ++j)
        {
            const float timeOffset = j < offsets.Size() ? offsets[j] : 0.0f;
            animationController->SetTime(animationControls[i].name_, time + timeOffset);
        }
        animationController->Update(0);
        node.MarkDirty();

        // Write tracks
        for (Node* node : nodes)
        {
            if (AnimationTrack* track = result->GetTrack(node->GetName()))
            {
                AnimationKeyFrame keyFrame;
                keyFrame.time_ = time;
                keyFrame.position_ = node->GetPosition();
                keyFrame.rotation_ = node->GetRotation();
                keyFrame.scale_ = node->GetScale();
                track->AddKeyFrame(keyFrame);
            }
        }
    }

    return result;
}

/// Merge times of animation tracks.
PODVector<float> MergeAnimationTrackTimes(const PODVector<AnimationTrack*>& tracks)
{
    PODVector<float> result;
    for (AnimationTrack* track : tracks)
    {
        if (track)
        {
            for (unsigned i = 0; i < track->GetNumKeyFrames(); ++i)
                result.Push(track->GetKeyFrame(i)->time_);
        }
    }
    Sort(result.Begin(), result.End());
    result.Erase(RandomAccessIterator<float>(std::unique(result.Begin().ptr_, result.End().ptr_)), result.End());
    return result;
}

/// Non-recursively get children bones by parent name.
PODVector<Bone*> GetChildren(Skeleton& skeleton, const String& parentName)
{
    PODVector<Bone*> result;
    const Bone* parent = skeleton.GetBone(parentName);
    const unsigned parentIndex = parent - skeleton.GetBone((unsigned)0);
    for (unsigned i = 0; i < skeleton.GetNumBones(); ++i)
        if (Bone* bone = skeleton.GetBone(i))
            if (bone->parentIndex_ == parentIndex)
                result.Push(bone);
    return result;
}

/// Names of foot bones.
struct FootBoneNames
{
    String thigh_;
    String calf_;
    String heel_;
};

/// Get names of thigh, calf and heel bones.
FootBoneNames GetFootBones(Skeleton& skeleton, const String& thighName)
{
    FootBoneNames result;
    result.thigh_ = thighName;
    PODVector<Bone*> thighChildren = GetChildren(skeleton, result.thigh_);
    if (!thighChildren.Empty())
    {
        result.calf_ = thighChildren[0]->name_;
        PODVector<Bone*> calfChildren = GetChildren(skeleton, result.calf_);
        if (!calfChildren.Empty())
            result.heel_ = calfChildren[0]->name_;
    }
    return result;
}

/// Rotate node to make child match its position. Returns true if successfully matched.
bool MatchChildPosition(Node& parent, Node& child, const Vector3& newChildPosition)
{
    const Vector3 parentPosition = parent.GetWorldPosition();
    const Vector3 childPosition = child.GetWorldPosition();
    const Quaternion thighRotation(childPosition - parentPosition, newChildPosition - parentPosition);
    parent.SetWorldRotation(thighRotation * parent.GetWorldRotation());

    return Equals(0.0f, (newChildPosition - child.GetWorldPosition()).LengthSquared());
}

/// Animation state of single foot animation.
struct FootAnimationState
{
    /// Thigh position.
    Vector3 thighPosition_;
    /// Heel position.
    Vector3 heelPosition_;
    /// Knee rotation.
    float kneeRotation_;
};

/// Resolved animation state of single foot animation.
struct FootAnimationStateResolved
{
    /// Thigh position.
    Vector3 thighPosition_;
    /// Calf position.
    Vector3 calfPosition_;
    /// Heel position.
    Vector3 heelPosition_;
};

/// Foot animation key frame.
struct FootAnimationKeyFrame
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

/// Foot animation track.
struct FootAnimationTrack
{
    Vector3 initialDirection_;
    /// Key frames.
    Vector<FootAnimationKeyFrame> keyFrames_;
    /// Static ranges.
    Vector<Pair<float, float>> staticRanges_;

    bool IsStatic(float time) const
    {
        for (const Pair<float, float>& range : staticRanges_)
            if (range.first_ <= time && time <= range.second_)
                return true;
        return false;
    }
    float GetMomementRange(float time) const
    {
        float result = 0.0f;
        for (const Pair<float, float>& range : staticRanges_)
            if (time < range.first_)
                return range.first_ - time;
        return staticRanges_.Empty() ? 1.0f : staticRanges_.Front().first_ - time + GetLength();
    }
    /// Get length.
    float GetLength() const
    {
        return keyFrames_.Back().time_ - keyFrames_.Front().time_;
    }
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
    FootAnimationKeyFrame SampleFrame(float time, unsigned& frame) const
    {
        GetKeyFrameIndex(time, frame);
        const unsigned nextFrame = frame + 1 < keyFrames_.Size() ? frame + 1 : 0;
        const FootAnimationKeyFrame* keyFrame = &keyFrames_[frame];
        const FootAnimationKeyFrame* nextKeyFrame = &keyFrames_[nextFrame];
        float timeInterval = nextKeyFrame->time_ - keyFrame->time_;
        if (timeInterval < 0.0f)
            timeInterval += GetLength();

        float t = timeInterval > 0.0f ? (time - keyFrame->time_) / timeInterval : 1.0f;

        FootAnimationKeyFrame result;
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

/// Append track times to result array.
void AppendTrackTimes(PODVector<float>& result, AnimationTrack& track)
{
    for (unsigned i = 0; i < track.GetNumKeyFrames(); ++i)
        result.Push(track.GetKeyFrame(i)->time_);
}

/// Get signed angle between two vectors.
float AngleSigned(const Vector3& lhs, const Vector3& rhs, const Vector3& base)
{
    return lhs.Angle(rhs) * (lhs.CrossProduct(rhs).DotProduct(base) < 0 ? 1 : -1);
}

/// Mix in quaternion to another with specified weight. Returns new weight of result quaternion.
Quaternion MixQuaternion(Quaternion& lhs, const Quaternion& rhs, float weight, float totalWeight)
{
    // Nothing to mix at right
    if (weight < M_EPSILON)
        return lhs;

    // Nothing to mix at left
    if (totalWeight < M_EPSILON)
        return rhs;

    return lhs.Slerp(rhs, weight / (weight + totalWeight));
}

void CreateFootAnimationTrack(FootAnimationTrack& track, Model* model, Animation* animation, const String& thighName,
    const Vector3& velocity, float threshold)
{
    if (!model)
        return;

    // Setup temporary node for playback
    Node node(model->GetContext());
    AnimatedModel* animatedModel = node.CreateComponent<AnimatedModel>();
    animatedModel->SetModel(model);
    AnimationState* animationState = animatedModel->AddAnimationState(animation);
    animationState->SetWeight(1.0f);

    // Parse skeleton
    Skeleton& skeleton = animatedModel->GetSkeleton();
    const FootBoneNames boneNames = GetFootBones(skeleton, thighName);

    // Get related nodes and bones
    Node* thighNode = node.GetChild(boneNames.thigh_, true);
    Node* calfNode = thighNode ? thighNode->GetChild(boneNames.calf_) : nullptr;
    Node* heelNode = calfNode ? calfNode->GetChild(boneNames.heel_) : nullptr;
    Bone* thighBone = skeleton.GetBone(boneNames.thigh_);
    Bone* calfBone = skeleton.GetBone(boneNames.calf_);
    Bone* heelBone = skeleton.GetBone(boneNames.heel_);
    if (!thighNode || !calfNode || !heelNode || !thighBone || !calfBone || !heelBone)
        return;

    // Get initial pose
    track.initialDirection_ = heelNode->GetWorldPosition() - thighNode->GetWorldPosition();

    // Get tracks
    PODVector<float> times;
    AnimationTrack* thighTrack = animation->GetTrack(boneNames.thigh_);
    AnimationTrack* calfTrack = animation->GetTrack(boneNames.calf_);
    AnimationTrack* heelTrack = animation->GetTrack(boneNames.heel_);
    AppendTrackTimes(times, *thighTrack);
    AppendTrackTimes(times, *calfTrack);
    AppendTrackTimes(times, *heelTrack);

    Sort(times.Begin(), times.End());
    times.Erase(RandomAccessIterator<float>(std::unique(times.Begin().ptr_, times.End().ptr_)), times.End());

    // Play animation and convert pose
    const unsigned numKeyFrames = times.Size();
    track.keyFrames_.Resize(numKeyFrames);
    Vector<Vector3> globalPositions;
    float minHeight = M_INFINITY;
    for (unsigned i = 0; i < numKeyFrames; ++i)
    {
        // Play animation
        animationState->SetTime(times[i]);
        animationState->Apply();
        node.MarkDirty();

        // Receive animation data
        const Vector3 thighPosition = thighNode->GetWorldPosition();
        const Vector3 calfPosition = calfNode->GetWorldPosition();
        const Vector3 heelPosition = heelNode->GetWorldPosition();
        const Quaternion thighRotation = thighNode->GetRotation();
        const Quaternion calfRotation = calfNode->GetRotation();

        // Get knee direction
        const Vector3 direction = (heelPosition - thighPosition).Normalized();
        const Vector3 jointProjection = direction * (calfPosition - thighPosition).ProjectOntoAxis(direction) + thighPosition;
        const Vector3 jointDirection = Quaternion(direction, track.initialDirection_) * (calfPosition - jointProjection);

        // Revert to bind pose
        thighNode->SetTransform(thighBone->initialPosition_, thighBone->initialRotation_, thighBone->initialScale_);
        calfNode->SetTransform(calfBone->initialPosition_, calfBone->initialRotation_, calfBone->initialScale_);

        // Try to resolve foot shape and generate fix angles
        MatchChildPosition(*thighNode, *calfNode, calfPosition);
        const Quaternion thighRotationFix = thighNode->GetRotation().Inverse() * thighRotation;
        thighNode->SetRotation(thighNode->GetRotation() * thighRotationFix);

        MatchChildPosition(*calfNode, *heelNode, heelPosition);
        const Quaternion calfRotationFix = calfNode->GetRotation().Inverse() * calfRotation;
        calfNode->SetRotation(calfNode->GetRotation() * calfRotationFix);

        // Make new frame
        track.keyFrames_[i].time_ = times[i];
        track.keyFrames_[i].heelPosition_ = heelPosition;
        track.keyFrames_[i].kneeDirection_ = jointDirection.LengthSquared() > M_EPSILON ? jointDirection.Normalized() : Vector3::FORWARD;
        track.keyFrames_[i].thighRotationFix_ = thighRotationFix;
        track.keyFrames_[i].calfRotationFix_ = calfRotationFix;
        track.keyFrames_[i].heelRotationLocal_ = heelNode->GetRotation();
        track.keyFrames_[i].heelRotationWorld_ = heelNode->GetWorldRotation();
        globalPositions.Push(heelPosition + velocity * times[i]);
        minHeight = Min(minHeight, globalPositions.Back().y_);
    }

    float rangeBegin = -1;
    bool wasStatic = false;
    for (unsigned i = 0; i < numKeyFrames; ++i)
    {
        const bool isStatic = globalPositions[i].y_ < minHeight + threshold;
        if (isStatic && !wasStatic)
            rangeBegin = times[i];
        else if (wasStatic && !isStatic)
            track.staticRanges_.Push(MakePair(rangeBegin, times[i]));
        wasStatic = isStatic;
    }
    if (wasStatic)
        track.staticRanges_.Push(MakePair(rangeBegin, times.Back()));
}

//////////////////////////////////////////////////////////////////////////

/// Intersect sphere and sphere. If there is no intersection point, second sphere is moved toward first one.
void IntersectSphereSphereGuaranteed(const Sphere& first, const Sphere& second, float& distance, float& radius)
{
    // http://mathworld.wolfram.com/Sphere-SphereIntersection.html
    const float R = first.radius_;
    const float r = second.radius_;
    const float d = Min(R + r, (second.center_ - first.center_).Length());
    radius = Sqrt(Max(0.0f, (-d + r - R) * (-d - r + R) * (-d + r + R) * (d + r + R))) / (2 * d);
    distance = Sqrt(R * R - radius * radius);
}

/// Compute joint orientation.
Vector3 ComputeJointOrientation(const Vector3& baseOrientation, const Vector3& baseDirection, const Vector3& currentDirection,
    const Matrix3x4& rootTransform)
{
    const Quaternion rootRotation = rootTransform.Rotation();
    const Quaternion rotation(baseDirection, rootRotation.Inverse() * currentDirection);
    return rootRotation * rotation * baseOrientation;
}

/// Return vector clamped by sphere.
Vector3 ClampVector(const Vector3& vec, const Vector3& origin, float radius)
{
    const Vector3 delta = vec - origin;
    const float length = delta.Length();
    return length <= radius ? vec : origin + delta / length * radius;
}

/// Get bone index.
unsigned GetBoneIndex(Skeleton& skeleton, const String& boneName)
{
    Bone* bone = skeleton.GetBone(boneName);
    return bone ? bone - &skeleton.GetBones().Front() : skeleton.GetNumBones();
}

/// Compute global transforms of the skeleton - implementation.
void ComputeGlobalTransforms(Skeleton& skeleton, unsigned boneIndex, const Matrix3x4& baseTransform,
    PODVector<Matrix3x4>& transforms, PODVector<bool>& dirty)
{
    // Don't compute twice
    if (!dirty[boneIndex])
        return;

    // Compute parent if dirty
    const Bone* bone = skeleton.GetBone(boneIndex);
    if (bone->parentIndex_ != boneIndex)
        ComputeGlobalTransforms(skeleton, bone->parentIndex_, baseTransform, transforms, dirty);

    // Compute global transform
    const Matrix3x4 localTransform(bone->initialPosition_, bone->initialRotation_, bone->initialScale_);
    const Matrix3x4 parentTransform = bone->parentIndex_ == boneIndex ? baseTransform : transforms[bone->parentIndex_];
    transforms[boneIndex] = parentTransform * localTransform;
    dirty[boneIndex] = false;
}

/// Compute global transforms of the skeleton.
PODVector<Matrix3x4> ComputeGlobalTransforms(Skeleton& skeleton, const Matrix3x4& baseTransform)
{
    const unsigned numBones = skeleton.GetNumBones();
    PODVector<Matrix3x4> transforms(numBones);
    PODVector<bool> dirty(numBones);
    for (unsigned i = 0; i < numBones; ++i)
        ComputeGlobalTransforms(skeleton, i, baseTransform, transforms, dirty);
    return transforms;
}

/// Project point onto segment.
Vector3 ProjectPointOntoSegment(const Vector3& point, const Vector3& from, const Vector3& to)
{
    const Vector3 direction = (to - from).Normalized();
    return direction * (point - from).ProjectOntoAxis(direction) + from;
}

/// Make vector orthogonal to axis.
Vector3 OrthogonalizeVector(const Vector3& vec, const Vector3& axis)
{
    return axis.CrossProduct(vec).CrossProduct(axis).Normalized();
}

/// Resolve knee position.
Vector3 ResolveKneePosition(const Vector3& thighPosition, const Vector3& targetHeelPosition, const Vector3& jointDirection,
    float thighLength, float calfLength)
{
    float distance;
    float radius;
    IntersectSphereSphereGuaranteed(Sphere(thighPosition, thighLength), Sphere(targetHeelPosition, calfLength), distance, radius);
    const Vector3 direction = (targetHeelPosition - thighPosition).Normalized();
    const Vector3 orthoJointDirection = OrthogonalizeVector(jointDirection, direction);
    return thighPosition + direction * distance + orthoJointDirection * radius;
}

/// Get rotation axis and angle.
void GetAxisAngle(const Quaternion& rotation, Vector3& axis, float& angle)
{
    angle = 2 * Acos(rotation.w_);

    const float ww = rotation.w_ * rotation.w_;
    axis.x_ = rotation.x_ / sqrt(1 - ww);
    axis.y_ = rotation.y_ / sqrt(1 - ww);
    axis.z_ = rotation.z_ / sqrt(1 - ww);
}

/// Get rotation axis.
Vector3 GetAxis(const Quaternion& rotation)
{
    Vector3 axis;
    float angle;
    GetAxisAngle(rotation, axis, angle);
    return axis;
}

/// Get rotation angle.
float GetAngle(const Quaternion& rotation)
{
    Vector3 axis;
    float angle;
    GetAxisAngle(rotation, axis, angle);
    return Fract(angle / 360) * 360;
}

/// Linear interpolation of angles.
float LerpAngle(float lhs, float rhs, float t)
{
    lhs = Fract(lhs / 360) * 360;
    rhs = Fract(rhs / 360) * 360;

    const float diff = Abs(lhs - rhs);
    return (diff < 180)
        ? Lerp(lhs, rhs, t)
        : (lhs > rhs)
            ? Lerp(lhs - 360, rhs, t)
            : Lerp(lhs, rhs - 360, t);
}

/**
Decompose the rotation on to 2 parts.
1. Twist - rotation around the "direction" vector
2. Swing - rotation around axis that is perpendicular to "direction" vector
The rotation can be composed back by
rotation = swing * twist

has singularity in case of swing_rotation close to 180 degrees rotation.
if the input quaternion is of non-unit length, the outputs are non-unit as well
otherwise, outputs are both unit

http://stackoverflow.com/questions/3684269/component-of-a-quaternion-rotation-around-an-axis
*/
void GetSwingTwist(const Quaternion& rotation, const Vector3& direction, Quaternion& swing, Quaternion& twist)
{
    Vector3 ra( rotation.x_, rotation.y_, rotation.z_ ); // rotation axis
    Vector3 p = ra.ProjectOntoAxis(direction) * ra; // return projection v1 on to v2  (parallel component)
    twist = Quaternion(rotation.w_, p.x_, p.y_, p.z_).Normalized();
    swing = rotation * twist.Conjugate();
}

/// Compute timestamps of animation.
PODVector<float> GatherTimeStamps(Animation& animation)
{
    PODVector<float> timeStamps;
    for (const auto& item : animation.GetTracks())
    {
        const AnimationTrack& track = item.second_;
        if (!timeStamps.Empty() && track.keyFrames_.Size() != timeStamps.Size())
        {
            URHO3D_LOGERROR("CharacterAnimation has inconsistent number of key frames in tracks");
            return PODVector<float>();
        }

        if (!track.keyFrames_.Empty())
        {
            timeStamps.Resize(track.keyFrames_.Size());
            for (unsigned i = 0; i < track.GetNumKeyFrames(); ++i)
                timeStamps[i] = track.keyFrames_[i].time_;
        }
    }
    if (timeStamps.Empty())
        URHO3D_LOGWARNING("CharacterAnimation is empty");
    return timeStamps;
}

//////////////////////////////////////////////////////////////////////////
CharacterSkeletonSegmentType GetCharacterSkeletonSegmentType(const String& name)
{
    if (name == "root")
        return CharacterSkeletonSegmentType::Root;
    else if (name == "chain")
        return CharacterSkeletonSegmentType::Chain;
    else if (name == "limb")
        return CharacterSkeletonSegmentType::Limb;
    return CharacterSkeletonSegmentType::Root;
}

void CharacterRootSegmentData::Prepare(const CharacterSkeletonSegment& /*segment*/)
{

}

//////////////////////////////////////////////////////////////////////////
void CharacterRootSegmentData::Reset()
{
    position_ = Vector3::ZERO;
    rotation_ = Quaternion::IDENTITY;
}

void CharacterRootSegmentData::Blend(const CharacterRootSegmentData& other, float weight)
{
    position_ = position_.Lerp(other.position_, weight);
    rotation_ = rotation_.Slerp(other.rotation_, weight);
}

void CharacterRootSegmentData::Import(const CharacterSkeletonSegment& src)
{
    position_ = src.nodes_[0]->GetWorldPosition() - src.initialPose_[0].Translation();
    rotation_ = src.nodes_[0]->GetWorldRotation() * src.initialPose_[0].Rotation().Inverse();
}

void CharacterRootSegmentData::Export(
    const Matrix3x4& rootTransform, const Quaternion& animRotation, CharacterSkeletonSegment& dest) const
{
    const Quaternion transformedRotation = animRotation.Inverse() * rotation_ * animRotation;
    const Vector3 initialPosition = dest.initialPose_[0].Translation();
    const Quaternion initialRotation = dest.initialPose_[0].Rotation();

    Node& node = *dest.nodes_[0];
    node.SetWorldPosition(rootTransform * (initialPosition + animRotation * position_));
    node.SetWorldRotation(rootTransform.Rotation() * transformedRotation * initialRotation);
}

bool CharacterRootSegmentData::Save(XMLElement& dest) const
{
    dest.SetVector3("position", position_);
    dest.SetQuaternion("rotation", rotation_);
    return true;
}

bool CharacterRootSegmentData::Load(const XMLElement& src)
{
    position_ = src.GetVector3("position");
    rotation_ = src.GetQuaternion("rotation");
    return true;
}

//////////////////////////////////////////////////////////////////////////
void CharacterChainSegmentData::Prepare(const CharacterSkeletonSegment& segment)
{
    rotations_.Resize(segment.bones_.Size());
}

void CharacterChainSegmentData::Reset()
{
    position_ = Vector3::ZERO;
    for (Quaternion& rotation : rotations_)
        rotation = Quaternion::IDENTITY;
}

void CharacterChainSegmentData::Blend(const CharacterChainSegmentData& other, float weight)
{
    position_ = position_.Lerp(other.position_, weight);
    rotations_.Resize(other.rotations_.Size());
    for (unsigned i = 0; i < rotations_.Size(); ++i)
        rotations_[i] = rotations_[i].Slerp(other.rotations_[i], weight);
}

void CharacterChainSegmentData::Import(const CharacterSkeletonSegment& src)
{
    position_ = src.nodes_.Back()->GetWorldPosition();
    rotations_.Resize(src.nodes_.Size());
    for (unsigned i = 0; i < src.nodes_.Size(); ++i)
        rotations_[i] = src.nodes_[i]->GetWorldRotation() * src.initialPose_[i].Rotation().Inverse();
}

void CharacterChainSegmentData::Export(
    const Matrix3x4& rootTransform, const Quaternion& animRotation, CharacterSkeletonSegment& dest) const
{
    for (unsigned i = 0; i < dest.nodes_.Size(); ++i)
    {
        const Quaternion transformedRotation = animRotation.Inverse() * rotations_[i] * animRotation;
        dest.nodes_[i]->SetWorldRotation(rootTransform.Rotation() * transformedRotation * dest.initialPose_[i].Rotation());
    }
}

bool CharacterChainSegmentData::Save(XMLElement& dest) const
{
    dest.SetVector3("position", position_);
    for (unsigned j = 0; j < rotations_.Size(); ++j)
        dest.CreateChild("bone").SetQuaternion("rotation", rotations_[j]);
    return true;
}

bool CharacterChainSegmentData::Load(const XMLElement& src)
{
    position_ = src.GetVector3("position");
    for (XMLElement bone = src.GetChild("bone"); bone; bone = bone.GetNext())
        rotations_.Push(bone.GetQuaternion("rotation"));
    return true;
}

//////////////////////////////////////////////////////////////////////////
void CharacterLimbSegmentData::Prepare(const CharacterSkeletonSegment& /*segment*/)
{

}

void CharacterLimbSegmentData::Reset()
{
    position_ = Vector3::ZERO;
    rotation_ = 0.0f;
    rotationA_ = 0.0f;
    rotationB_ = 0.0f;
    rotationC_ = Quaternion::IDENTITY;
}

void CharacterLimbSegmentData::Blend(const CharacterLimbSegmentData& other, float weight)
{
    position_ = position_.Lerp(other.position_, weight);
    rotation_ = LerpAngle(rotation_, other.rotation_, weight);
    rotationA_ = LerpAngle(rotationA_, other.rotationA_, weight);
    rotationB_ = LerpAngle(rotationB_, other.rotationB_, weight);
    rotationC_ = rotationC_.Slerp(other.rotationC_, weight);
}

void CharacterLimbSegmentData::Import(const CharacterSkeletonSegment& src)
{
    Node* nodeA = src.nodes_[0];
    Node* nodeB = src.nodes_[1];
    Node* nodeC = src.nodes_[2];

    const Bone& boneA = *src.bones_[0];
    const Bone& boneB = *src.bones_[1];
    const Bone& boneC = *src.bones_[2];

    const Matrix3x4 transformA = nodeA->GetWorldTransform();
    const Matrix3x4 transformB = nodeB->GetWorldTransform();
    const Matrix3x4 transformC = nodeC->GetWorldTransform();

    const Matrix3x4& initialA = src.initialPose_[0];
    const Matrix3x4& initialB = src.initialPose_[1];
    const Matrix3x4& initialC = src.initialPose_[2];

    // Get position
    const float firstLength = (transformB.Translation() - transformA.Translation()).Length();
    const float secondLength = (transformC.Translation() - transformB.Translation()).Length();
    position_ = (transformC.Translation() - initialA.Translation()) / (firstLength + secondLength);

    // Get reference bend direction
    const Vector3 initialDirection = initialC.Translation() - initialA.Translation();
    const Vector3 currentDirection = transformC.Translation() - transformA.Translation();
    const Vector3 referenceBendDirection = OrthogonalizeVector(
        Quaternion(initialDirection, currentDirection) * src.jointDirection_, currentDirection);

    // Get actual bend direction
    const Vector3 positionBproj = ProjectPointOntoSegment(
        transformB.Translation(), transformA.Translation(), transformC.Translation());
    const Vector3 actualBendDirection = OrthogonalizeVector(
        transformB.Translation() - positionBproj, currentDirection);

    // Compute and revert limb rotation
    const Quaternion bendDirectionRotation(referenceBendDirection, actualBendDirection);
    rotation_ = GetAngle(bendDirectionRotation);
    rotation_ *= Sign(currentDirection.DotProduct(GetAxis(bendDirectionRotation)));
    nodeA->SetWorldRotation(bendDirectionRotation.Inverse() * nodeA->GetWorldRotation());

    // Revert transforms to initial pose
    const Matrix3x4 adjustedInitialA = nodeA->GetParent()->GetWorldTransform()
        * Matrix3x4(boneA.initialPosition_, boneA.initialRotation_, boneA.initialScale_);
    const Matrix3x4 adjustedInitialB = adjustedInitialA
        * Matrix3x4(boneB.initialPosition_, boneB.initialRotation_, boneB.initialScale_);
    const Matrix3x4 adjustedInitialC = adjustedInitialB
        * Matrix3x4(boneC.initialPosition_, boneC.initialRotation_, boneC.initialScale_);

    if (!MatchChildPosition(*nodeA, *nodeB, adjustedInitialB.Translation()))
        URHO3D_LOGWARNING("Failed to resolve thigh-calf segment of foot animation");
    if (!MatchChildPosition(*nodeB, *nodeC, adjustedInitialC.Translation()))
        URHO3D_LOGWARNING("Failed to resolve calf-heel segment of foot animation");

    // Gather local rotations
    const Quaternion beforeA = nodeA->GetWorldRotation();
    const Quaternion beforeB = nodeB->GetWorldRotation();

    const Quaternion pureRotationA = nodeA->GetRotation() * boneA.initialRotation_.Inverse();
    const Quaternion pureRotationB = nodeB->GetRotation() * boneB.initialRotation_.Inverse();

    rotationA_ = GetAngle(pureRotationA);
    rotationB_ = GetAngle(pureRotationB);
    rotationC_ = transformC.Rotation() * initialC.Rotation().Inverse();

    // Flip angles if needed
    const Vector3 directionAB = nodeB->GetPosition().Normalized();
    const Vector3 directionBC = nodeC->GetPosition().Normalized();
    const Vector3 axisAB = GetAxis(pureRotationA);
    const Vector3 axisBC = GetAxis(pureRotationB);

//     Quaternion swingA, twistA, swingB, twistB;
//     GetSwingTwist(pureRotationA, directionAB, swingA, twistA);
//     GetSwingTwist(pureRotationB, directionBC, swingB, twistB);
//     rotationA_ = GetAngle(twistA);
//     rotationB_ = GetAngle(twistB);

    rotationA_ *= Sign(axisAB.DotProduct(directionAB));
    rotationB_ *= Sign(axisBC.DotProduct(directionBC));

    // #TODO
    {
        CharacterSkeletonSegment& dest = const_cast<CharacterSkeletonSegment&>(src);

        // Revert nodes to initial pose
        nodeA->SetRotation(Quaternion(rotationA_, boneB.initialPosition_) * boneA.initialRotation_);
        nodeB->SetRotation(Quaternion(rotationB_, boneC.initialPosition_) * boneB.initialRotation_);
        nodeC->SetRotation(boneC.initialRotation_);
//         nodeA->SetWorldRotation(initialA.Rotation());
//         nodeB->SetWorldRotation(initialB.Rotation());
//         nodeC->SetWorldRotation(initialC.Rotation());

        // Apply local segment rotations
//         const Vector3 directionAB = (nodeB->GetWorldPosition() - nodeA->GetWorldPosition()).Normalized();
//         const Vector3 directionBC = (nodeC->GetWorldPosition() - nodeB->GetWorldPosition()).Normalized();
//         const Quaternion worldRotationA = Quaternion(rotationA_, directionAB) * nodeA->GetWorldRotation();
//         const Quaternion worldRotationB = Quaternion(rotationB_, directionBC) * nodeB->GetWorldRotation();
//         nodeA->SetWorldRotation(worldRotationA);
//         nodeB->SetWorldRotation(worldRotationB);

        const Quaternion afterA = nodeA->GetWorldRotation();
        const Quaternion afterB = nodeB->GetWorldRotation();

        int t = 0;
    }
}

void CharacterLimbSegmentData::Export(
    const Matrix3x4& rootTransform, const Quaternion& animRotation, CharacterSkeletonSegment& dest) const
{
    Bone& boneA = *dest.bones_[0];
    Bone& boneB = *dest.bones_[1];
    Bone& boneC = *dest.bones_[2];

    Node& nodeA = *dest.nodes_[0];
    Node& nodeB = *dest.nodes_[1];
    Node& nodeC = *dest.nodes_[2];

    const Matrix3x4& initialA = dest.initialPose_[0];
    const Matrix3x4& initialB = dest.initialPose_[1];
    const Matrix3x4& initialC = dest.initialPose_[2];

    const float thighLength = (nodeA.GetWorldPosition() - nodeB.GetWorldPosition()).Length();
    const float calfLength = (nodeB.GetWorldPosition() - nodeC.GetWorldPosition()).Length();
    const Vector3 worldPosition = rootTransform * (animRotation * (position_ * (thighLength + calfLength)) + initialA.Translation());

    // Revert nodes to initial pose
//     nodeA.SetRotation(boneA.initialRotation_);
//     nodeB.SetRotation(boneB.initialRotation_);
    nodeA.SetRotation(Quaternion(rotationA_, boneB.initialPosition_) * boneA.initialRotation_);
    nodeB.SetRotation(Quaternion(rotationB_, boneC.initialPosition_) * boneB.initialRotation_);
    nodeC.SetRotation(boneC.initialRotation_);

    // Apply local segment rotations
//     const Vector3 directionAB = (nodeB.GetWorldPosition() - nodeA.GetWorldPosition()).Normalized();
//     const Vector3 directionBC = (nodeC.GetWorldPosition() - nodeB.GetWorldPosition()).Normalized();
//     const Quaternion worldRotationA = Quaternion(rotationA_, directionAB) * nodeA.GetWorldRotation();
//     const Quaternion worldRotationB = Quaternion(rotationB_, directionBC) * nodeB.GetWorldRotation();
//     nodeA.SetWorldRotation(worldRotationA);
//     nodeB.SetWorldRotation(worldRotationB);

    // Compute bend direction
    const Vector3 initialDirection = rootTransform.RotationMatrix() * (initialC.Translation() - initialA.Translation());
    const Vector3 currentDirection = worldPosition - nodeA.GetWorldPosition();
    const Quaternion limbRotation(initialDirection, currentDirection);
    const Vector3 bendDirection = OrthogonalizeVector(
        limbRotation * rootTransform.Rotation() * dest.jointDirection_, currentDirection);

    // Resolve limb
    const Vector3 worldPos0 = nodeA.GetWorldPosition();
    const Vector3 worldPos2 = ClampVector(worldPosition, worldPos0, thighLength + calfLength);
    const Vector3 worldPos1 = ResolveKneePosition(worldPos0, worldPos2, bendDirection, thighLength, calfLength);

    // Apply limb bending
    if (!MatchChildPosition(nodeA, nodeB, worldPos1))
        URHO3D_LOGWARNING("Failed to resolve thigh-calf segment of foot animation");

    if (!MatchChildPosition(nodeB, nodeC, worldPos2))
        URHO3D_LOGWARNING("Failed to resolve calf-heel segment of foot animation");

    // Apply limb rotation
    nodeA.SetWorldRotation(Quaternion(rotation_, worldPos2 - worldPos0) * nodeA.GetWorldRotation());

//     if (dest.name_ == "RightFoot")
//         URHO3D_LOGINFOF("%f %f | %f %f | %f %f",
//             worldPosition.x_, worldPosition.y_, worldPos2.x_, worldPos2.y_,
//             nodeC.GetWorldPosition().x_, nodeC.GetWorldPosition().y_);
//     if (dest.name_ == "RightFoot")
//         URHO3D_LOGINFOF("%f %f | %f %f",
//             nodeC.GetWorldPosition().x_, nodeC.GetWorldPosition().y_,
//             rootTransform.Translation().x_, position_.z_);

    // Apply target rotation
    nodeC.SetWorldRotation(rootTransform.Rotation() * animRotation.Inverse() * rotationC_ * animRotation * initialC.Rotation());
}

bool CharacterLimbSegmentData::Save(XMLElement& dest) const
{
    dest.SetVector3("position", position_);
    dest.SetFloat("rotation", rotation_);
    dest.SetFloat("rotation0", rotationA_);
    dest.SetFloat("rotation1", rotationB_);
    dest.SetQuaternion("rotation2", rotationC_);
    return true;
}

bool CharacterLimbSegmentData::Load(const XMLElement& src)
{
    position_ = src.GetVector3("position");
    rotation_ = src.GetFloat("rotation");
    rotationA_ = src.GetFloat("rotation0");
    rotationB_ = src.GetFloat("rotation1");
    rotationC_ = src.GetQuaternion("rotation2");
    return true;
}

//////////////////////////////////////////////////////////////////////////
void CharacterSkeleton::RegisterObject(Context* context)
{
    context->RegisterFactory<CharacterSkeleton>();
}

bool CharacterSkeleton::BeginLoad(Deserializer& source)
{
    SharedPtr<XMLFile> xmlFile = MakeShared<XMLFile>(context_);
    if (xmlFile->Load(source))
        return LoadXML(xmlFile->GetRoot());
    return false;
}

bool CharacterSkeleton::LoadXML(const XMLElement& source)
{
    for (XMLElement child = source.GetChild(); child; child = child.GetNext())
    {
        CharacterSkeletonSegment segment;
        segment.name_ = child.GetAttribute("name");
        segment.type_ = GetCharacterSkeletonSegmentType(child.GetName());
        segment.boneNames_ = child.GetAttribute("bones").Split(' ');
        segment.jointDirection_ = child.GetVector3("joint");
        segments_.Push(segment);
    }
    return true;
}

HashMap<String, Vector3> CharacterSkeleton::ComputeDirections(Animation& animation, Model& model) const
{
    const PODVector<float> timeStamps = GatherTimeStamps(animation);

    // Setup temporary node for playback
    Node node(context_);
    AnimatedModel* animatedModel = node.CreateComponent<AnimatedModel>();
    animatedModel->SetModel(&model);
    AnimationState* animationState = animatedModel->AddAnimationState(&animation);
    animationState->SetWeight(1.0f);
    Skeleton& skeleton = animatedModel->GetSkeleton();

    Vector<CharacterSkeletonSegment> segmentsData;
    AllocateSegmentData(segmentsData, skeleton, Matrix3x4::IDENTITY);

    // Collect directions
    HashMap<String, Vector3> directions;
    for (unsigned i = 0; i < timeStamps.Size(); ++i)
    {
        // Play animation
        animationState->SetTime(timeStamps[i]);
        animationState->Apply();
        node.MarkDirty();

        // Read tracks
        for (const CharacterSkeletonSegment& segment : segmentsData)
        {
            if (segment.type_ != CharacterSkeletonSegmentType::Limb)
                continue;
            const Node& nodeA = *segment.nodes_[0];
            const Node& nodeB = *segment.nodes_[1];
            const Node& nodeC = *segment.nodes_[2];
            const Quaternion initialRotationA = segment.initialPose_[0].Rotation();
            const Vector3 projB = ProjectPointOntoSegment(nodeB.GetWorldPosition(), nodeA.GetWorldPosition(), nodeC.GetWorldPosition());
            const Vector3 direction = nodeB.GetWorldPosition() - projB;
            directions[segment.name_] += initialRotationA * nodeA.GetWorldRotation().Inverse() * direction.Normalized();
        }
    }
    for (auto& segment : directions)
        segment.second_.Normalize();
    return directions;
}

const CharacterSkeletonSegment* CharacterSkeleton::FindSegment(const String& name) const
{
    for (const CharacterSkeletonSegment& segment : segments_)
        if (segment.name_ == name)
            return &segment;
    return nullptr;
}

bool CharacterSkeleton::AllocateSegmentData(Vector<CharacterSkeletonSegment>& segmentsData,
    Skeleton& skeleton, const Matrix3x4& baseTransform) const
{
    // #TODO Prevent crashes if AllocateSegmentData returns false
    const PODVector<Matrix3x4> globalTransforms = ComputeGlobalTransforms(skeleton, baseTransform);

    const unsigned numSegments = segments_.Size();
    segmentsData.Resize(numSegments);
    for (unsigned i = 0; i < numSegments; ++i)
    {
        // Allocate segment
        const CharacterSkeletonSegment& segment = segments_[i];
        segmentsData[i] = segment;

        // Fill bone data
        const unsigned numBones = segment.boneNames_.Size();
        segmentsData[i].bones_.Resize(numBones);
        segmentsData[i].jointDirection_ = baseTransform.RotationMatrix() * segment.jointDirection_;
        segmentsData[i].nodes_.Resize(numBones);
        segmentsData[i].initialPose_.Resize(numBones);
        for (unsigned j = 0; j < numBones; ++j)
        {
            const unsigned boneIndex = GetBoneIndex(skeleton, segment.boneNames_[j]);
            Bone* bone = skeleton.GetBone(boneIndex);
            if (!bone)
            {
                URHO3D_LOGERRORF("Cannot find bone %s of segment %s", segment.boneNames_[j].CString(), segment.name_.CString());
                return false;
            }

            segmentsData[i].bones_[j] = bone;
            segmentsData[i].nodes_[j] = bone->node_;
            segmentsData[i].initialPose_[j] = globalTransforms[boneIndex];
        }
    }
    return true;
}

//////////////////////////////////////////////////////////////////////////
SharedPtr<CharacterAnimationTrack> CharacterAnimationTrack::Create(Context* context, CharacterSkeletonSegmentType type)
{
    switch (type)
    {
    case Urho3D::CharacterSkeletonSegmentType::Root:
        return MakeShared<RootAnimationTrack>(context);
    case Urho3D::CharacterSkeletonSegmentType::Chain:
        return MakeShared<ChainAnimationTrack>(context);
    case Urho3D::CharacterSkeletonSegmentType::Limb:
        return MakeShared<LimbAnimationTrack>(context);
    default:
        return nullptr;
    }
}

void CharacterAnimationTrack::GetKeyFrameIndex(float time, unsigned& index) const
{
    if (time < 0.0f)
        time = 0.0f;

    if (index >= timeStamps_.Size())
        index = timeStamps_.Size() - 1;

    // Check for being too far ahead
    while (index && time < timeStamps_[index])
        --index;

    // Check for being too far behind
    while (index < timeStamps_.Size() - 1 && time >= timeStamps_[index + 1])
        ++index;
}

unsigned CharacterAnimationTrack::GetKeyFrameIndex(float time) const
{
    unsigned index = 0;
    GetKeyFrameIndex(time, index);
    return index;
}

void CharacterAnimationTrack::GetKeyFrame(float time, unsigned& firstFrame, unsigned& secondFrame, float& factor) const
{
    firstFrame = GetKeyFrameIndex(time);
    secondFrame = firstFrame + 1 < timeStamps_.Size() ? firstFrame + 1 : 0;
    const float timeInterval = Max(0.0f, timeStamps_[secondFrame] - timeStamps_[firstFrame]);
    factor = timeInterval > 0.0f ? (time - timeStamps_[firstFrame]) / timeInterval : 1.0f;
}

float CharacterAnimationTrack::GetLength() const
{
    return timeStamps_.Size() > 0 ? timeStamps_.Back() - timeStamps_.Front() : 0.0f;
}

//////////////////////////////////////////////////////////////////////////
void CharacterAnimation::RegisterObject(Context* context)
{
    context->RegisterFactory<CharacterAnimation>();
}

bool CharacterAnimation::BeginLoad(Deserializer& source)
{
    SharedPtr<XMLFile> xmlFile = MakeShared<XMLFile>(context_);
    if (xmlFile->Load(source))
        return LoadXML(xmlFile->GetRoot());
    return false;
}

bool CharacterAnimation::LoadXML(const XMLElement& source)
{
    for (XMLElement child = source.GetChild(); child; child = child.GetNext())
    {
        const CharacterSkeletonSegmentType type = GetCharacterSkeletonSegmentType(child.GetName());
        if (SharedPtr<CharacterAnimationTrack> track = CharacterAnimationTrack::Create(context_, type))
        {
            track->SetName(child.GetAttribute("name"));
            if (!track->LoadXML(child))
                return false;
            tracks_.Push(track);
        }
    }
    return true;
}

bool CharacterAnimation::SaveXML(XMLElement& dest) const
{
    for (CharacterAnimationTrack* track : tracks_)
    {
        XMLElement node = dest.CreateChild(track->GetTypeString());
        node.SetAttribute("name", track->GetName());
        if (!track->SaveXML(node))
            return false;
    }

    return true;
}

CharacterAnimationTrack* CharacterAnimation::FindTrack(const String& name) const
{
    for (CharacterAnimationTrack* track : tracks_)
        if (track->GetName() == name)
            return track;
    return nullptr;
}

bool CharacterAnimation::Save(Serializer& dest) const
{
    SharedPtr<XMLFile> xml(new XMLFile(context_));
    XMLElement rootNode = xml->CreateRoot("animation");
    SaveXML(rootNode);
    return xml->Save(dest);
}

#if 0
bool CharacterAnimation::ImportAnimation(CharacterSkeleton& characterSkeleton, Model& model, Animation& animation)
{
    // Setup temporary node for playback
    Node node(context_);
    AnimatedModel* animatedModel = node.CreateComponent<AnimatedModel>();
    animatedModel->SetModel(&model);
    AnimationState* animationState = animatedModel->AddAnimationState(&animation);
    animationState->SetWeight(1.0f);
    Skeleton& skeleton = animatedModel->GetSkeleton();

    // Create tracks of 2-segments
    Segment2TrackMap segments2;
    for (const CharacterSkeleton::Segment2Map::KeyValue& elem : characterSkeleton.GetSegments2())
    {
        const CharacterSkeletonSegment2& joint = elem.second_;

        // Get nodes and bones of segment
        Node* segmentRootNode = node.GetChild(joint.rootBone_, true);
        Node* segmentJointNode = segmentRootNode ? segmentRootNode->GetChild(joint.jointBone_) : nullptr;
        Node* segmentTargetNode = segmentJointNode ? segmentJointNode->GetChild(joint.targetBone_) : nullptr;
        Bone* thighBone = skeleton.GetBone(segmentRootNode->GetName());
        Bone* calfBone = skeleton.GetBone(segmentJointNode->GetName());
        Bone* heelBone = skeleton.GetBone(segmentTargetNode->GetName());
        if (!segmentRootNode || !segmentJointNode || !segmentTargetNode || !thighBone || !calfBone || !heelBone)
        {
            URHO3D_LOGERRORF("Failed to load 2-segment '%s' of character skeleton: root='%s', joint='%s', target='%s'",
                joint.name_.CString(), joint.rootBone_.CString(), joint.jointBone_.CString(), joint.targetBone_.CString());
            return false;
        }

        // Get initial pose
        CharacterAnimationSegment2Track track;
        track.name_ = joint.name_;
        track.initialDirection_ = segmentTargetNode->GetWorldPosition() - segmentRootNode->GetWorldPosition();

        // Get sample times
        const PODVector<float> sampleTimes = MergeAnimationTrackTimes(
        {
            animation.GetTrack(joint.rootBone_),
            animation.GetTrack(joint.jointBone_),
            animation.GetTrack(joint.targetBone_)
        });

        // Play animation and convert pose
        const unsigned numKeyFrames = sampleTimes.Size();
        track.keyFrames_.Resize(numKeyFrames);
        for (unsigned i = 0; i < numKeyFrames; ++i)
        {
            // Play animation
            animationState->SetTime(sampleTimes[i]);
            animationState->Apply();
            node.MarkDirty();

            // Receive animation data
            const Vector3 thighPosition = segmentRootNode->GetWorldPosition();
            const Vector3 calfPosition = segmentJointNode->GetWorldPosition();
            const Vector3 heelPosition = segmentTargetNode->GetWorldPosition();
            const Quaternion thighRotation = segmentRootNode->GetRotation();
            const Quaternion calfRotation = segmentJointNode->GetRotation();

            // Get knee direction
            const Vector3 direction = (heelPosition - thighPosition).Normalized();
            const Vector3 jointProjection = direction * (calfPosition - thighPosition).ProjectOntoAxis(direction) + thighPosition;
            const Vector3 jointDirection = Quaternion(direction, track.initialDirection_) * (calfPosition - jointProjection);

            // Revert to bind pose
            segmentRootNode->SetTransform(thighBone->initialPosition_, thighBone->initialRotation_, thighBone->initialScale_);
            segmentJointNode->SetTransform(calfBone->initialPosition_, calfBone->initialRotation_, calfBone->initialScale_);

            // Try to resolve foot shape and generate fix angles
            MatchChildPosition(*segmentRootNode, *segmentJointNode, calfPosition);
            const Quaternion thighRotationFix = segmentRootNode->GetRotation().Inverse() * thighRotation;
            segmentRootNode->SetRotation(segmentRootNode->GetRotation() * thighRotationFix);

            MatchChildPosition(*segmentJointNode, *segmentTargetNode, heelPosition);
            const Quaternion calfRotationFix = segmentJointNode->GetRotation().Inverse() * calfRotation;
            segmentJointNode->SetRotation(segmentJointNode->GetRotation() * calfRotationFix);

            // Make new frame
            track.keyFrames_[i].time_ = sampleTimes[i];
            track.keyFrames_[i].heelPosition_ = heelPosition;
            track.keyFrames_[i].kneeDirection_ = jointDirection.LengthSquared() > M_EPSILON ? jointDirection.Normalized() : Vector3::FORWARD;
            track.keyFrames_[i].thighRotationFix_ = thighRotationFix;
            track.keyFrames_[i].calfRotationFix_ = calfRotationFix;
            track.keyFrames_[i].heelRotationLocal_ = segmentTargetNode->GetRotation();
            track.keyFrames_[i].heelRotationWorld_ = segmentTargetNode->GetWorldRotation();
        }
        segments2.Populate(elem.first_, track);
    }
    segments2_.Insert(segments2);
    return true;
}
#endif

bool CharacterAnimation::Import(Animation& animation, Model& model, CharacterSkeleton& rig, const Matrix3x4& transform)
{
    // Read timestamps
    const PODVector<float> timeStamps = GatherTimeStamps(animation);
    if (timeStamps.Empty())
        return false;

    // Setup temporary node for playback
    Node node(context_);
    node.SetTransform(transform);
    AnimatedModel* animatedModel = node.CreateComponent<AnimatedModel>();
    animatedModel->SetModel(&model);
    AnimationState* animationState = animatedModel->AddAnimationState(&animation);
    animationState->SetWeight(1.0f);
    Skeleton& skeleton = animatedModel->GetSkeleton();

    // Create tracks
    for (const CharacterSkeletonSegment& segment : rig.GetSegments())
    {
        tracks_.Push(CharacterAnimationTrack::Create(context_, segment.type_));
        tracks_.Back()->SetName(segment.name_);
    }

    // Gather segments data
    Vector<CharacterSkeletonSegment> segmentData;
    rig.AllocateSegmentData(segmentData, skeleton, transform);

    // Collect tracks
    for (unsigned i = 0; i < timeStamps.Size(); ++i)
    {
        // Play animation
        animationState->SetTime(timeStamps[i]);
        animationState->Apply();
        node.MarkDirty();

        // Read tracks
        for (unsigned j = 0; j < segmentData.Size(); ++j)
            if (tracks_[j])
                tracks_[j]->ImportFrame(timeStamps[i], segmentData[j]);
        // #TODO Backup
    }

    return true;
}

//////////////////////////////////////////////////////////////////////////
void CharacterEffector::RegisterObject(Context* context)
{
    URHO3D_COPY_BASE_ATTRIBUTES(Component);
    URHO3D_ATTRIBUTE("Segment", String, segmentName_, String::EMPTY, AM_DEFAULT);
    URHO3D_ATTRIBUTE("Is Animated", bool, animated_, false, AM_DEFAULT);
}

void CharacterEffector::PrepareEffector(CharacterSkeletonSegment& segment)
{
    for (Bone* bone : segment.bones_)
        bone->node_->SetTransform(bone->initialPosition_, bone->initialRotation_, bone->initialScale_);
}

//////////////////////////////////////////////////////////////////////////
CharacterAnimationController::CharacterAnimationController(Context* context)
    : AnimationController(context)
{
}

CharacterAnimationController::~CharacterAnimationController()
{
}

void CharacterAnimationController::RegisterObject(Context* context)
{
    context->RegisterFactory<CharacterAnimationController>(characterAnimatorCategory);
    URHO3D_COPY_BASE_ATTRIBUTES(AnimationController);
    URHO3D_MIXED_ACCESSOR_ATTRIBUTE("Skeleton", GetSkeletonAttr, SetSkeletonAttr, ResourceRef, ResourceRef(XMLFile::GetTypeStatic()), AM_DEFAULT);
    URHO3D_ATTRIBUTE("Enable Animation", bool, animationEnabled_, true, AM_DEFAULT);
    URHO3D_ATTRIBUTE("Animation Rotation", Quaternion, animationRotation_, Quaternion::IDENTITY, AM_DEFAULT);
    URHO3D_ATTRIBUTE("Revert Transform", bool, revertAnimationTransform_, true, AM_DEFAULT);
}

void CharacterAnimationController::SetAnimationRotation(const Quaternion& rotation)
{
    if (animationRotation_ != rotation)
    {
        animationRotation_ = rotation;
        MarkDirty();
    }
}

void CharacterAnimationController::SetTargetTransform(StringHash segment, const Matrix3x4& transform)
{
    Segment2State& state = segment2states_[segment];
    state.targetTransform_ = transform;
}

void CharacterAnimationController::SetTargetRotationAmount(StringHash segment, float rotationAmount)
{
    Segment2State& state = segment2states_[segment];
    state.targetRotationAmount_ = rotationAmount;
}

void CharacterAnimationController::SetTargetRotationBalance(StringHash segment, float globalFactor)
{
    Segment2State& state = segment2states_[segment];
    state.globalRotationFactor_ = globalFactor;
}

void CharacterAnimationController::CleanSegment2(StringHash segment)
{
    segment2states_.Erase(segment);
}

void CharacterAnimationController::Update(float timeStep)
{
    CheckIntegrity();
    UpdateHierarchy();
    AnimationController::Update(timeStep);
    if (animatedModel_ && animationEnabled_)
    {
        animatedModel_->ApplyAnimation();
        ApplyAnimation();
    }
}

void CharacterAnimationController::ApplyAnimation()
{
    // Gather animation data
    currentAnimationData_.Clear();
    for (const AnimationControl& animationControl : GetAnimations())
        if (CharacterAnimation* characterAnimation = GetCharacterAnimation(animationControl.name_))
            if (AnimationState* animationState = GetAnimationState(animationControl.name_))
            {
                CachedAnimationData data;
                data.characterAnimation_ = characterAnimation;
                data.time_ = animationState->GetTime();
                data.weight_ = animationState->GetWeight();
                currentAnimationData_.Push(data);
            }

    // Animate segments
    for (auto& elem : segments_)
    {
        CharacterSkeletonSegment& segment = elem.first_;
        CharacterEffector& effector = *elem.second_;
        if (!effector.IsEnabledEffective())
            continue;

        effector.PrepareEffector(segment);
        if (effector.IsAnimated() && !currentAnimationData_.Empty())
        {
            effector.ResetEffector();
            for (auto& animationData : currentAnimationData_)
                if (CharacterAnimationTrack* track = animationData.characterAnimation_->FindTrack(segment.name_))
                    effector.ApplyAnimation(animationData.weight_, animationData.time_, *track);
        }
        effector.ResolveEffector(*this, segment);
    }

    // Revert animation transform
    if (revertAnimationTransform_)
    {
        rootNode_->SetRotation(animationRotation_.Inverse() * rootNode_->GetRotation());
        rootNode_->SetPosition(animationRotation_.Inverse() * rootNode_->GetPosition());
    }
}

void CharacterAnimationController::SetSkeleton(CharacterSkeleton* skeleton)
{
    skeleton_ = skeleton;
    animatedModel_.Reset();
}

void CharacterAnimationController::SetSkeletonAttr(const ResourceRef& value)
{
    ResourceCache* cache = GetSubsystem<ResourceCache>();
    skeleton_ = cache->GetResource<CharacterSkeleton>(value.name_);
}

Urho3D::ResourceRef CharacterAnimationController::GetSkeletonAttr() const
{
    return ResourceRef(XMLFile::GetTypeStatic(), GetResourceName(skeleton_));
}

CharacterEffector* CharacterAnimationController::CreateEffector(Node& node, CharacterSkeletonSegmentType type)
{
    switch (type)
    {
    case CharacterSkeletonSegmentType::Root:
        return node.CreateComponent<CharacterRootEffector>(LOCAL);
    case CharacterSkeletonSegmentType::Limb:
        return node.CreateComponent<CharacterLimbEffector>(LOCAL);
    case CharacterSkeletonSegmentType::Chain:
        return node.CreateComponent<CharacterChainEffector>(LOCAL);
    default:
        return nullptr;
    }
}

void CharacterAnimationController::CheckIntegrity()
{
    static const char* containerName = "[segments]";

    // Skip check if attributes are not initialized
    if (!node_ || !skeleton_)
        return;

    // If model is expired, mark dirty
    if (!animatedModel_ || !rootNode_)
        MarkDirty();

    // Try to get model
    if (!animatedModel_)
    {
        animatedModel_ = node_->GetComponent<AnimatedModel>(true);
        if (!animatedModel_)
            return;
    }

    // If segments are not initialized, mark dirty
    const unsigned numSegments = skeleton_->GetSegments().Size();
    if (segments_.Size() != numSegments)
        MarkDirty();

    // If any segment is expired, mark dirty
    if (!dirty_)
    {
        for (const auto& segment : segments_)
            if (!segment.second_)
            {
                MarkDirty();
                break;
            }
    }

    if (dirty_)
    {
        rootNode_ = animatedModel_->GetSkeleton().GetRootBone()->node_;

        // Remove old nodes
        if (Node* container = node_->GetChild(containerName))
            node_->RemoveChild(container);

        // Reset node transform
        const Matrix3x4 nodeTransform = node_->GetWorldTransform();
        node_->SetWorldTransform(Vector3::ZERO, animationRotation_, Vector3::ONE);

        // Initialize segments
        // #TODO Rename segmentData_
        skeleton_->AllocateSegmentData(segmentData_, animatedModel_->GetSkeleton(), Matrix3x4(Vector3::ZERO, animationRotation_, 1.0f));

        // Create nodes
        Node* container = node_->CreateChild(containerName, LOCAL);
        segments_.Resize(numSegments);
        for (unsigned i = 0; i < numSegments; ++i)
        {
            const CharacterSkeletonSegment& segment = segmentData_[i];
            Node* child = container->CreateChild(segment.name_, LOCAL);
            CharacterEffector* effector = CreateEffector(*child, segment.type_);
            effector->SetAnimated(true);
            effector->InitializeEffector(segment);

            segments_[i].first_ = segment;
            segments_[i].second_ = effector;
        }

        // Restore node transform
        node_->SetWorldTransform(nodeTransform.Translation(), nodeTransform.Rotation(), nodeTransform.Scale());

        // Reset dirty
        dirty_ = false;
    }
}

void CharacterAnimationController::UpdateHierarchy()
{
    // Get component
//     if (!node_ || !skeleton_)
//         return;
//     if (animatedModel_ && &animatedModel_->GetSkeleton() == animatedModelSkeleton_)
//         return;
//     animatedModel_ = node_->GetComponent<AnimatedModel>();
//     animatedModelSkeleton_ = &animatedModel_->GetSkeleton();
//     if (!animatedModel_)
//         return;
//
//     skeleton_->AllocateSegmentData(segmentData_, *animatedModelSkeleton_, Matrix3x4::IDENTITY);
}

CharacterAnimation* CharacterAnimationController::GetCharacterAnimation(const String& animationName)
{
    if (animationCache_.Contains(animationName))
        return animationCache_[animationName];

    ResourceCache* cache = GetSubsystem<ResourceCache>();
    CharacterAnimation* characterAnimation = cache->GetResource<CharacterAnimation>(ReplaceExtension(animationName, ".char.xml"));
    animationCache_[animationName] = characterAnimation;
    return characterAnimation;
}

CharacterSkeletonSegment* CharacterAnimationController::GetSegment(const String& segmentName)
{
    for (CharacterSkeletonSegment& segment : segmentData_)
        if (segment.name_ == segmentName)
            return &segment;
    return nullptr;
}

void CharacterAnimationController::UpdateSegment2(const CharacterSkeletonSegment2& segment)
{
    // Get nodes and bones
    Skeleton& skeleton = animatedModel_->GetSkeleton();
    Bone* rootBone = skeleton.GetBone(segment.rootBone_);
    Bone* jointBone = skeleton.GetBone(segment.jointBone_);
    Bone* targetBone = skeleton.GetBone(segment.targetBone_);
    if (!rootBone || !jointBone || !targetBone)
    {
        URHO3D_LOGERRORF("Skeleton segment '%s' is not found in the model", segment.name_.CString());
        return;
    }

    Node* thighNode = rootBone->node_;
    Node* calfNode = jointBone->node_;
    Node* heelNode = targetBone->node_;

    thighNode->SetRotationSilent(rootBone->initialRotation_);
    calfNode->SetRotationSilent(jointBone->initialRotation_);
    heelNode->SetRotationSilent(targetBone->initialRotation_);
    thighNode->MarkDirty();

    // Apply animations
    float accumulatedWeight = 0.0f;
    Vector3 baseDirection;
    CharacterAnimationSegment2KeyFrame keyFrame;
    keyFrame.time_ = -1;
    for (const AnimationControl& animationControl : GetAnimations())
    {
        CharacterAnimation* characterAnimation = GetCharacterAnimation(animationControl.name_);
        if (!characterAnimation)
            continue;
        AnimationState* animationState = GetAnimationState(animationControl.name_);
        if (!animationState)
            continue;
        CharacterAnimationSegment2Track* track = 0;// characterAnimation->FindTrack(segment.name_);
        if (!track)
            continue;

        // Apply animation to key frame
        unsigned frameIndex = 0;
        const CharacterAnimationSegment2KeyFrame animationFrame = track->SampleFrame(animationState->GetTime(), frameIndex);
        float factor = animationState->GetWeight();
        keyFrame.heelPosition_ += animationFrame.heelPosition_ * factor;
        keyFrame.kneeDirection_ += animationFrame.kneeDirection_ * factor;
        keyFrame.thighRotationFix_ = MixQuaternion(keyFrame.thighRotationFix_, animationFrame.thighRotationFix_, factor, accumulatedWeight);
        keyFrame.calfRotationFix_ = MixQuaternion(keyFrame.calfRotationFix_, animationFrame.calfRotationFix_, factor, accumulatedWeight);
        keyFrame.heelRotationLocal_ = MixQuaternion(keyFrame.heelRotationLocal_, animationFrame.heelRotationLocal_, factor, accumulatedWeight);
        keyFrame.heelRotationWorld_ = MixQuaternion(keyFrame.heelRotationWorld_, animationFrame.heelRotationWorld_, factor, accumulatedWeight);
        baseDirection += track->initialDirection_ * factor;
        accumulatedWeight += factor;
    }

    const float thighLength = (thighNode->GetWorldPosition() - calfNode->GetWorldPosition()).Length();
    const float calfLength = (calfNode->GetWorldPosition() - heelNode->GetWorldPosition()).Length();

    // Get state
    const Segment2State& state = segment2states_[segment.name_];

    // Resolve limb
    const Matrix3x4 animationToWorld = node_->GetWorldTransform() /** animationTransform_*/;
    const Matrix3x4 worldToAnimation = animationToWorld.Inverse();

    const Vector3 worldPos0 = thighNode->GetWorldPosition();
    const Vector3 worldPos2 = ClampVector(animationToWorld * state.targetTransform_ * keyFrame.heelPosition_, worldPos0, thighLength + calfLength);
    const Vector3 worldJointDirection = worldPos2 - worldPos0;
    const Vector3 animJointDirection = worldToAnimation.RotationMatrix() * worldJointDirection;
    const Vector3 animJointOrientation = Quaternion(baseDirection, animJointDirection) * keyFrame.kneeDirection_;
    const Vector3 worldJointOrientation = animationToWorld.RotationMatrix() * animJointOrientation;
    const Vector3 worldPos1 = ResolveKneePosition(worldPos0, worldPos2, worldJointOrientation, thighLength, calfLength);

    // Apply foot shape
    if (!MatchChildPosition(*thighNode, *calfNode, worldPos1))
        URHO3D_LOGWARNING("Failed to resolve thigh-calf segment of foot animation");
    Vector3 q1 = keyFrame.thighRotationFix_.EulerAngles();
    thighNode->SetRotation(thighNode->GetRotation() * keyFrame.thighRotationFix_);

    if (!MatchChildPosition(*calfNode, *heelNode, worldPos2))
        URHO3D_LOGWARNING("Failed to resolve calf-heel segment of foot animation");
    Vector3 q2 = keyFrame.calfRotationFix_.EulerAngles();
    calfNode->SetRotation(calfNode->GetRotation() * keyFrame.calfRotationFix_);

    // Resolve heel rotation
    const Quaternion origHeelRotation = calfNode->GetWorldRotation() * keyFrame.heelRotationLocal_;
    const Quaternion fixedHeelRotation = node_->GetWorldRotation() * keyFrame.heelRotationWorld_;
    const Quaternion adjustToGoundRotation = Quaternion::IDENTITY.Slerp(state.targetTransform_.Rotation(), state.targetRotationAmount_);
    heelNode->SetWorldRotation(adjustToGoundRotation * origHeelRotation.Slerp(fixedHeelRotation, state.globalRotationFactor_));

    thighNode->MarkDirty();
}

//////////////////////////////////////////////////////////////////////////
void CharacterRootEffector::RegisterObject(Context* context)
{
    context->RegisterFactory<CharacterRootEffector>(characterAnimatorCategory);
    URHO3D_COPY_BASE_ATTRIBUTES(CharacterEffector);
    URHO3D_ATTRIBUTE("Position", Vector3, effectorState_.position_, Vector3::ZERO, AM_DEFAULT);
    URHO3D_ATTRIBUTE("Rotation", Quaternion, effectorState_.rotation_, Quaternion::IDENTITY, AM_DEFAULT);
}

//////////////////////////////////////////////////////////////////////////
void CharacterLimbEffector::RegisterObject(Context* context)
{
    context->RegisterFactory<CharacterLimbEffector>(characterAnimatorCategory);
    URHO3D_COPY_BASE_ATTRIBUTES(CharacterEffector);
    URHO3D_ATTRIBUTE("Position", Vector3, effectorState_.position_, Vector3::ZERO, AM_DEFAULT);
    URHO3D_ATTRIBUTE("Base Rotation", float, effectorState_.rotation_, 0.0f, AM_DEFAULT);
    URHO3D_ATTRIBUTE("Twirl 1", float, effectorState_.rotationA_, 0.0f, AM_DEFAULT);
    URHO3D_ATTRIBUTE("Twirl 2", float, effectorState_.rotationB_, 0.0f, AM_DEFAULT);
}

//////////////////////////////////////////////////////////////////////////
void CharacterChainEffector::RegisterObject(Context* context)
{
    context->RegisterFactory<CharacterChainEffector>(characterAnimatorCategory);
    URHO3D_COPY_BASE_ATTRIBUTES(CharacterEffector);
}

//////////////////////////////////////////////////////////////////////////
void RegisterCharacterAnimator(Context* context)
{
    CharacterSkeleton::RegisterObject(context);
    CharacterAnimation::RegisterObject(context);
    CharacterAnimationController::RegisterObject(context);

    CharacterEffector::RegisterObject(context);
    CharacterRootEffector::RegisterObject(context);
    CharacterLimbEffector::RegisterObject(context);
    CharacterChainEffector::RegisterObject(context);
}

void DumpCharacterSkeletonDirections(const String& skeletonName, const String& animationName, const String& modelName)
{
    ResourceCache* cache = GetScriptContext()->GetSubsystem<ResourceCache>();
    if (!cache)
        return;
    CharacterSkeleton* skeleton = cache->GetResource<CharacterSkeleton>(skeletonName);
    Animation* animation = cache->GetResource<Animation>(animationName);
    Model* model = cache->GetResource<Model>(modelName);
    if (!skeleton || !animation || !model)
        return;

    const HashMap<String, Vector3> directions = skeleton->ComputeDirections(*animation, *model);
    for (const auto& item : directions)
        URHO3D_LOGINFOF("%s: (%f, %f, %f)", item.first_.CString(), item.second_.x_, item.second_.y_, item.second_.z_);
}

bool ImportCharacterAnimation(const String& outputName,
    const String& animationName, const String& skeletonName, const String& modelName, const Matrix3x4& animationTransform)
{
    ResourceCache* cache = GetScriptContext()->GetSubsystem<ResourceCache>();
    if (!cache)
        return false;

    SharedPtr<Model> model(cache->GetResource<Model>(modelName));
    SharedPtr<Animation> animation(cache->GetResource<Animation>(animationName));
    SharedPtr<CharacterSkeleton> characterSkeleton(cache->GetResource<CharacterSkeleton>(skeletonName));
    if (!model)
    {
        URHO3D_LOGWARNINGF("Failed to load model %s", modelName.CString());
        return false;
    }
    if (!animation)
    {
        URHO3D_LOGWARNINGF("Failed to load animation %s", animationName.CString());
        return false;
    }
    if (!characterSkeleton)
    {
        URHO3D_LOGWARNINGF("Failed to load skeleton %s", skeletonName.CString());
        return false;
    }

    CharacterAnimation characterAnimation(GetScriptContext());
    characterAnimation.SetName(outputName);
    if (!characterAnimation.Import(*animation, *model, *characterSkeleton, animationTransform))
    {
        URHO3D_LOGWARNING("Failed to import animation");
        return false;
    }

    if (!characterAnimation.SaveFile(outputName))
    {
        URHO3D_LOGWARNINGF("Failed to save destination file %s", outputName.CString());
        return false;
    }
    return true;
}

void OverrideModelScale(const String& outputName, const String& modelName, const Vector3& scale)
{
    ResourceCache* cache = GetScriptContext()->GetSubsystem<ResourceCache>();
    if (!cache)
        return;

    SharedPtr<Model> model(cache->GetResource<Model>(modelName));
    if (!model)
        return;

    Bone* rootBone = model->GetSkeleton().GetRootBone();
    if (!rootBone)
        return;

    rootBone->initialPosition_ *= scale / rootBone->initialScale_;
    rootBone->initialScale_ = scale;
    model->SaveFile(outputName);
}

void OverrideAnimationScale(const String& outputName, const String& animationName, const String& modelName, const Vector3& scale)
{
    ResourceCache* cache = GetScriptContext()->GetSubsystem<ResourceCache>();
    if (!cache)
        return;

    SharedPtr<Model> model(cache->GetResource<Model>(modelName));
    SharedPtr<Animation> animation(cache->GetResource<Animation>(animationName));
    if (!model || !animation)
        return;

    Bone* rootBone = model->GetSkeleton().GetRootBone();
    if (!rootBone)
        return;

    AnimationTrack* track = animation->GetTrack(rootBone->name_);
    if (!track)
        return;

    const bool hasPosition = !!(track->channelMask_ & CHANNEL_POSITION);
    const bool hasScale = !!(track->channelMask_ & CHANNEL_SCALE);

    // Update position
    if (hasPosition)
    {
        for (AnimationKeyFrame& keyFrame : track->keyFrames_)
            keyFrame.position_ *= scale / (hasScale ? keyFrame.scale_ : Vector3::ONE);
    }

    // Write scale
    track->channelMask_ |= CHANNEL_SCALE;
    for (AnimationKeyFrame& keyFrame : track->keyFrames_)
        keyFrame.scale_ = scale;

    animation->SaveFile(outputName);
}

void ResetRootAnimationTrackPosition(const String& outputName, const String& animationName, const String& modelName, float from, float to)
{
    ResourceCache* cache = GetScriptContext()->GetSubsystem<ResourceCache>();
    if (!cache)
        return;

    SharedPtr<Model> model(cache->GetResource<Model>(modelName));
    SharedPtr<Animation> animation(cache->GetResource<Animation>(animationName));
    if (!model || !animation)
        return;

    Bone* rootBone = model->GetSkeleton().GetRootBone();
    if (!rootBone)
        return;

    AnimationTrack* track = animation->GetTrack(rootBone->name_);
    if (!track)
        return;

    if (track->channelMask_ & CHANNEL_POSITION)
    {
        for (AnimationKeyFrame& keyFrame : track->keyFrames_)
            if (keyFrame.time_ >= from && keyFrame.time_ <= to)
                keyFrame.position_ = rootBone->initialPosition_;
    }

    animation->SaveFile(outputName);
}

void TrimAnimation(const String& outputName, const String& animationName, float from, float to)
{
    if (from > to)
        return;

    ResourceCache* cache = GetScriptContext()->GetSubsystem<ResourceCache>();
    if (!cache)
        return;

    SharedPtr<Animation> animation(cache->GetResource<Animation>(animationName));
    if (!animation)
        return;

    for (unsigned i = 0; i < animation->GetNumTracks(); ++i)
        if (AnimationTrack* track = animation->GetTrack(i))
        {
            unsigned fromIndex = 0;
            track->GetKeyFrameIndex(from, fromIndex);
            unsigned toIndex = 0;
            track->GetKeyFrameIndex(to, toIndex);
            track->keyFrames_.Erase(toIndex + 1, track->keyFrames_.Size() - toIndex - 1);
            track->keyFrames_.Erase(0, fromIndex);

            const float baseTime = track->keyFrames_.Front().time_;
            for (AnimationKeyFrame& keyFrame : track->keyFrames_)
                keyFrame.time_ -= baseTime;

            // Needed only once
            if (i == 0)
                animation->SetLength(track->keyFrames_.Back().time_ - track->keyFrames_.Front().time_);
        }

    animation->SaveFile(outputName);
}

float ComputeSegmentLength(const String& modelName, const String& skeletonName, const String& segmentName)
{
    ResourceCache* cache = GetScriptContext()->GetSubsystem<ResourceCache>();
    if (!cache)
        return 0;

    // Load resources
    SharedPtr<Model> model(cache->GetResource<Model>(modelName));
    SharedPtr<CharacterSkeleton> characterSkeleton(cache->GetResource<CharacterSkeleton>(skeletonName));
    if (!model)
    {
        URHO3D_LOGWARNINGF("Failed to load model %s", modelName.CString());
        return 0;
    }
    if (!characterSkeleton)
    {
        URHO3D_LOGWARNINGF("Failed to load skeleton %s", skeletonName.CString());
        return 0;
    }

    // Find segment
    const CharacterSkeletonSegment* segment = characterSkeleton->FindSegment(segmentName);
    if (!segment)
    {
        URHO3D_LOGWARNINGF("Can't find segment %s", segmentName.CString());
        return 0;
    }
    if (segment->type_ != CharacterSkeletonSegmentType::Limb || segment->boneNames_.Size() < 3)
    {
        URHO3D_LOGWARNINGF("Only limb segments allowed");
        return 0;
    }

    // Find bones
    Skeleton& skeleton = model->GetSkeleton();
    PODVector<Matrix3x4> globalTransforms = ComputeGlobalTransforms(skeleton, Matrix3x4::IDENTITY);
    Bone* boneA = skeleton.GetBone(segment->boneNames_[0]);
    Bone* boneB = skeleton.GetBone(segment->boneNames_[1]);
    Bone* boneC = skeleton.GetBone(segment->boneNames_[2]);
    if (!boneA || !boneB || !boneC)
    {
        URHO3D_LOGWARNINGF("Can't find bone");
        return 0;
    }

    // Compute length
    const Vector3 positionA = globalTransforms[boneA - skeleton.GetBones().Buffer()].Translation();
    const Vector3 positionB = globalTransforms[boneB - skeleton.GetBones().Buffer()].Translation();
    const Vector3 positionC = globalTransforms[boneC - skeleton.GetBones().Buffer()].Translation();
    return (positionC - positionB).Length() + (positionB - positionA).Length();
}

void CharacterAnimationController_SetTargetTransform(const String& segment, const Matrix3x4& transform,
    CharacterAnimationController* characterAnimationController)
{
    characterAnimationController->SetTargetTransform(segment, transform);
}

void CharacterAnimationController_SetTargetRotationAmount(const String& segment, float rotationAmount,
    CharacterAnimationController* characterAnimationController)
{
    characterAnimationController->SetTargetRotationAmount(segment, rotationAmount);
}

void CharacterAnimationController_SetTargetRotationBalance(const String& segment, float globalFactor,
    CharacterAnimationController* characterAnimationController)
{
    characterAnimationController->SetTargetRotationBalance(segment, globalFactor);
}

void RegisterCharacterAnimatorScriptAPI(asIScriptEngine* engine)
{
    RegisterResource<CharacterSkeleton>(engine, "CharacterSkeleton");

    engine->RegisterGlobalFunction("void DumpCharacterSkeletonDirections(const String&in, const String&in, const String&in)", asFUNCTION(DumpCharacterSkeletonDirections), asCALL_CDECL);
    engine->RegisterGlobalFunction("void ImportCharacterAnimation(const String&in, const String&in, const String&in, const String&in, const Matrix3x4&in)", asFUNCTION(ImportCharacterAnimation), asCALL_CDECL);
    engine->RegisterGlobalFunction("void OverrideModelScale(const String&in, const String&in, const Vector3&in)", asFUNCTION(OverrideModelScale), asCALL_CDECL);
    engine->RegisterGlobalFunction("void OverrideAnimationScale(const String&in, const String&in, const String&in, const Vector3&in)", asFUNCTION(OverrideAnimationScale), asCALL_CDECL);
    engine->RegisterGlobalFunction("void ResetRootAnimationTrackPosition(const String&in, const String&in, const String&in, float, float)", asFUNCTION(ResetRootAnimationTrackPosition), asCALL_CDECL);
    engine->RegisterGlobalFunction("void TrimAnimation(const String&in, const String&in, float, float)", asFUNCTION(TrimAnimation), asCALL_CDECL);
    engine->RegisterGlobalFunction("float ComputeSegmentLength(const String&in, const String&in, const String&in)", asFUNCTION(ComputeSegmentLength), asCALL_CDECL);

    RegisterComponent<CharacterAnimationController>(engine, "CharacterAnimationController");
    RegisterSubclass<CharacterAnimationController, AnimationController>(engine, "AnimationController", "CharacterAnimationController");
    engine->RegisterObjectMethod("CharacterAnimationController", "void SetTargetTransform(const String&in, const Matrix3x4&in)", asFUNCTION(CharacterAnimationController_SetTargetTransform), asCALL_CDECL_OBJLAST);
    engine->RegisterObjectMethod("CharacterAnimationController", "void SetTargetRotationAmount(const String&in, float)", asFUNCTION(CharacterAnimationController_SetTargetRotationAmount), asCALL_CDECL_OBJLAST);
    engine->RegisterObjectMethod("CharacterAnimationController", "void SetTargetRotationBalance(const String&in, float)", asFUNCTION(CharacterAnimationController_SetTargetRotationBalance), asCALL_CDECL_OBJLAST);
}

}

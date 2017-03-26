#include "CharacterAnimator.h"

#include <Urho3D/AngelScript/APITemplates.h>
#include <Urho3D/Core/Context.h>
#include <Urho3D/Container/Ptr.h>
#include <Urho3D/Graphics/AnimatedModel.h>
#include <Urho3D/Graphics/Animation.h>
#include <Urho3D/Graphics/AnimationState.h>
#include <Urho3D/Graphics/AnimationController.h>
#include <Urho3D/Graphics/DebugRenderer.h>
#include <Urho3D/IO/Log.h>
#include <Urho3D/Math/MathDefs.h>
#include <Urho3D/Math/Sphere.h>
#include <Urho3D/Resource/ResourceCache.h>

#include <algorithm>

namespace Urho3D
{

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

/// Resolve knee position.
Vector3 ResolveKneePosition(const Vector3& thighPosition, const Vector3& targetHeelPosition, const Vector3& jointDirection,
    float thighLength, float calfLength)
{
    float distance;
    float radius;
    IntersectSphereSphereGuaranteed(Sphere(thighPosition, thighLength), Sphere(targetHeelPosition, calfLength), distance, radius);
    const Vector3 direction = (targetHeelPosition - thighPosition).Normalized();
    const Vector3 orthoJointDirection = direction.CrossProduct(jointDirection).CrossProduct(direction).Normalized();
    return thighPosition + direction * distance + orthoJointDirection * radius;
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
    return angle;
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

//////////////////////////////////////////////////////////////////////////
CharacterSkeletonSegmentData* CharacterSkeletonSegmentData::Create(CharacterSkeletonSegmentType type)
{
    switch (type)
    {
    case CharacterSkeletonSegmentType::Root:
        return new CharacterSkeletonRootSegmentData();
    case CharacterSkeletonSegmentType::Chain:
        return new CharacterSkeletonChainSegmentData();
    case CharacterSkeletonSegmentType::Limb:
        return new CharacterSkeletonLimbSegmentData();
    default:
        return nullptr;
    }
}

void CharacterSkeletonSegmentData::Reset()
{
    accumulatedWeight_ = 0.0f;
}

//////////////////////////////////////////////////////////////////////////
void CharacterSkeletonRootSegmentData::Reset()
{
    CharacterSkeletonSegmentData::Reset();
    position_ = Vector3::ZERO;
    rotation_ = Quaternion::IDENTITY;
}

void CharacterSkeletonRootSegmentData::Merge(const CharacterSkeletonSegmentData& other, float weight)
{
    const CharacterSkeletonRootSegmentData& rhs = static_cast<const CharacterSkeletonRootSegmentData&>(other);
    const float balance = weight / (weight + accumulatedWeight_);
    position_ = position_.Lerp(rhs.position_, balance);
    rotation_ = rotation_.Slerp(rhs.rotation_, balance);
    accumulatedWeight_ += weight;
}

void CharacterSkeletonRootSegmentData::Apply(const Matrix3x4& rootTransform, CharacterSkeletonSegment& dest)
{
    const Vector3 localPosition = dest.globalPositions_[0] + position_;
    const Quaternion localRotation = rotation_ * dest.globalRotations_[0];
    const Matrix3x4 localTransform(localPosition, localRotation, dest.nodes_[0]->GetScale());
    const Matrix3x4& parentTransform = dest.nodes_[0]->GetParent()->GetWorldTransform();
    dest.nodes_[0]->SetTransform(parentTransform.Inverse() * rootTransform * localTransform);
}

//////////////////////////////////////////////////////////////////////////
void CharacterSkeletonChainSegmentData::Reset()
{
    CharacterSkeletonSegmentData::Reset();
    position_ = Vector3::ZERO;
    for (Quaternion& rotation : rotations_)
        rotation = Quaternion::IDENTITY;
}

void CharacterSkeletonChainSegmentData::Merge(const CharacterSkeletonSegmentData& other, float weight)
{
    const CharacterSkeletonChainSegmentData& rhs = static_cast<const CharacterSkeletonChainSegmentData&>(other);
    const float balance = weight / (weight + accumulatedWeight_);
    position_ = position_.Lerp(rhs.position_, balance);
    rotations_.Resize(rhs.rotations_.Size());
    for (unsigned i = 0; i < rotations_.Size(); ++i)
        rotations_[i] = rotations_[i].Slerp(rhs.rotations_[i], balance);
    accumulatedWeight_ += weight;
}

void CharacterSkeletonChainSegmentData::Apply(const Matrix3x4& rootTransform, CharacterSkeletonSegment& dest)
{
    for (unsigned i = 0; i < dest.nodes_.Size(); ++i)
        dest.nodes_[i]->SetWorldRotation(rootTransform.Rotation() * rotations_[i] * dest.globalRotations_[i]);
}

//////////////////////////////////////////////////////////////////////////
void CharacterSkeletonLimbSegmentData::Reset()
{
    CharacterSkeletonSegmentData::Reset();
    position_ = Vector3::ZERO;
    direction_ = Vector3::ZERO;
    rotationA_ = 0.0f;
    rotationB_ = 0.0f;
    rotationC_ = Quaternion::IDENTITY;
}

void CharacterSkeletonLimbSegmentData::Merge(const CharacterSkeletonSegmentData& other, float weight)
{
    const CharacterSkeletonLimbSegmentData& rhs = static_cast<const CharacterSkeletonLimbSegmentData&>(other);
    const float balance = weight / (weight + accumulatedWeight_);
    position_ = position_.Lerp(rhs.position_, balance);
    direction_ = direction_.Lerp(rhs.direction_, balance);
    rotationA_ = Lerp(rotationA_, rhs.rotationA_, balance);
    rotationB_ = Lerp(rotationB_, rhs.rotationB_, balance);
    rotationC_ = rotationC_.Slerp(rhs.rotationC_, balance);
    accumulatedWeight_ += weight;
}

void CharacterSkeletonLimbSegmentData::Apply(const Matrix3x4& rootTransform, CharacterSkeletonSegment& dest)
{
    Node* node0 = dest.nodes_[0];
    Node* node1 = dest.nodes_[1];
    Node* node2 = dest.nodes_[2];

    const Matrix3x4& initialA = dest.initialTransforms_[0];
    const Matrix3x4& initialB = dest.initialTransforms_[1];
    const Matrix3x4& initialC = dest.initialTransforms_[2];

    const float thighLength = (node0->GetWorldPosition() - node1->GetWorldPosition()).Length();
    const float calfLength = (node1->GetWorldPosition() - node2->GetWorldPosition()).Length();

    // Apply rotation
    const Vector3 directionAB = (initialB.Translation() - initialA.Translation()).Normalized();
    const Vector3 directionBC = (initialC.Translation() - initialB.Translation()).Normalized();
    node0->SetWorldRotation(rootTransform.Rotation() * Quaternion(rotationA_, directionAB) * initialA.Rotation());
    node1->SetWorldRotation(rootTransform.Rotation() * Quaternion(rotationA_, directionAB) * initialB.Rotation());

    // Resolve limb
    const Vector3 worldPos0 = node0->GetWorldPosition();
    const Vector3 worldPos2 = ClampVector(rootTransform.ToMatrix3() * position_ + worldPos0, worldPos0, thighLength + calfLength);
    const Vector3 baseDirection = initialC.Translation() - initialA.Translation();
    const Vector3 worldJointOrientation = ComputeJointOrientation(direction_, baseDirection, worldPos2 - worldPos0, rootTransform);
    const Vector3 worldPos1 = ResolveKneePosition(worldPos0, worldPos2, worldJointOrientation, thighLength, calfLength);

    // Apply foot shape
    if (!MatchChildPosition(*node0, *node1, worldPos1))
        URHO3D_LOGWARNING("Failed to resolve thigh-calf segment of foot animation");

    if (!MatchChildPosition(*node1, *node2, worldPos2))
        URHO3D_LOGWARNING("Failed to resolve calf-heel segment of foot animation");

    // Resolve heel rotation
//     const Quaternion origHeelRotation = calfNode->GetWorldRotation() * keyFrame.heelRotationLocal_;
//     const Quaternion fixedHeelRotation = node_->GetWorldRotation() * keyFrame.heelRotationWorld_;
//     const Quaternion adjustToGoundRotation = Quaternion::IDENTITY.Slerp(state.targetTransform_.Rotation(), state.targetRotationAmount_);
//     heelNode->SetWorldRotation(adjustToGoundRotation * origHeelRotation.Slerp(fixedHeelRotation, state.globalRotationFactor_));
//
//     thighNode->MarkDirty();

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
        segments_.Push(segment);
    }
    return true;
}

const CharacterSkeletonSegment* CharacterSkeleton::FindSegment(const String& name) const
{
    for (const CharacterSkeletonSegment& segment : segments_)
        if (segment.name_ == name)
            return &segment;
    return nullptr;
}

bool CharacterSkeleton::AllocateSegmentData(Vector<CharacterSkeletonSegment>& segmentsData,
    Skeleton& skeleton, const Matrix3x4& baseTransform)
{
    const PODVector<Matrix3x4> globalTransforms = ComputeGlobalTransforms(skeleton, baseTransform);

    const unsigned numSegments = segments_.Size();
    segmentsData.Resize(numSegments);
    for (unsigned i = 0; i < numSegments; ++i)
    {
        // Allocate segment
        const CharacterSkeletonSegment& segment = segments_[i];
        segmentsData[i] = segment;
        segmentsData[i].data_.reset(CharacterSkeletonSegmentData::Create(segment.type_));
        if (!segmentsData[i].data_)
        {
            URHO3D_LOGERRORF("Cannot initialize data of segment %s", segment.name_.CString());
            return false;
        }

        // Fill bone data
        const unsigned numBones = segment.boneNames_.Size();
        segmentsData[i].bones_.Resize(numBones);
        segmentsData[i].nodes_.Resize(numBones);
        segmentsData[i].initialTransforms_.Resize(numBones);
        segmentsData[i].globalPositions_.Resize(numBones);
        segmentsData[i].globalRotations_.Resize(numBones);
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
            segmentsData[i].initialTransforms_[j] = globalTransforms[boneIndex];
            segmentsData[i].globalPositions_[j] = globalTransforms[boneIndex].Translation();
            segmentsData[i].globalRotations_[j] = globalTransforms[boneIndex].Rotation();
        }
    }
    return true;
}

//////////////////////////////////////////////////////////////////////////
SharedPtr<CharacterAnimationTrack> CharacterAnimationTrack::Create(CharacterSkeletonSegmentType type, const String& name)
{
    switch (type)
    {
    case Urho3D::CharacterSkeletonSegmentType::Root:
        return MakeShared<RootAnimationTrack>(name);
    case Urho3D::CharacterSkeletonSegmentType::Chain:
        return MakeShared<ChainAnimationTrack>(name);
    case Urho3D::CharacterSkeletonSegmentType::Limb:
        return MakeShared<LimbAnimationTrack>(name);
    default:
        return nullptr;
    }
}

//////////////////////////////////////////////////////////////////////////
void RootAnimationTrack::ImportFrame(const CharacterSkeletonSegment& segment)
{
    CharacterSkeletonRootSegmentData frame;
    frame.position_ = segment.nodes_[0]->GetWorldPosition() - segment.globalPositions_[0];
    frame.rotation_ = segment.nodes_[0]->GetWorldRotation() * segment.globalRotations_[0].Inverse();
    track_.Push(frame);
}

void RootAnimationTrack::MergeFrame(unsigned firstFrame, unsigned secondFrame, float factor,
    float weight, CharacterSkeletonSegmentData& dest)
{
    const CharacterSkeletonRootSegmentData& first = track_[firstFrame];
    const CharacterSkeletonRootSegmentData& second = track_[secondFrame];
    dest.Merge(first, weight * (1 - factor));
    dest.Merge(second, weight * factor);
}

bool RootAnimationTrack::SaveXML(XMLElement& dest) const
{
    const unsigned numKeys = track_.Size();
    for (unsigned i = 0; i < numKeys; ++i)
    {
        XMLElement child = dest.CreateChild("key");
        child.SetVector3("position", track_[i].position_);
        child.SetQuaternion("rotation", track_[i].rotation_);
    }
    return true;
}

bool RootAnimationTrack::LoadXML(const XMLElement& source)
{
    for (XMLElement child = source.GetChild("key"); child; child = child.GetNext())
    {
        CharacterSkeletonRootSegmentData frame;
        frame.position_ = child.GetVector3("position");
        frame.rotation_ = child.GetQuaternion("rotation");
        track_.Push(frame);
    }
    return true;
}

//////////////////////////////////////////////////////////////////////////
void ChainAnimationTrack::ImportFrame(const CharacterSkeletonSegment& segment)
{
    CharacterSkeletonChainSegmentData frame;
    frame.position_ = segment.nodes_.Back()->GetWorldPosition();
    frame.rotations_.Resize(segment.nodes_.Size());
    for (unsigned i = 0; i < segment.nodes_.Size(); ++i)
        frame.rotations_[i] = segment.nodes_[i]->GetWorldRotation() * segment.globalRotations_[i].Inverse();
    track_.Push(frame);
}

void ChainAnimationTrack::MergeFrame(unsigned firstFrame, unsigned secondFrame, float factor,
    float weight, CharacterSkeletonSegmentData& dest)
{
    const CharacterSkeletonChainSegmentData& first = track_[firstFrame];
    const CharacterSkeletonChainSegmentData& second = track_[secondFrame];
    dest.Merge(first, weight * (1 - factor));
    dest.Merge(second, weight * factor);
}

bool ChainAnimationTrack::SaveXML(XMLElement& dest) const
{
    const unsigned numKeys = track_.Size();
    for (unsigned i = 0; i < numKeys; ++i)
    {
        XMLElement child = dest.CreateChild("frame");
        child.SetVector3("position", track_[i].position_);
        for (unsigned j = 0; j < track_[i].rotations_.Size(); ++j)
            child.CreateChild("bone").SetQuaternion("rotation", track_[i].rotations_[j]);
    }
    return true;
}

bool ChainAnimationTrack::LoadXML(const XMLElement& source)
{
    for (XMLElement child = source.GetChild("frame"); child; child = child.GetNext())
    {
        CharacterSkeletonChainSegmentData frame;
        frame.position_ = child.GetVector3("position");
        for (XMLElement bone = child.GetChild("bone"); bone; bone = bone.GetNext())
            frame.rotations_.Push(bone.GetQuaternion("rotation"));
        track_.Push(frame);
    }
    return true;
}

//////////////////////////////////////////////////////////////////////////
void LimbAnimationTrack::ImportFrame(const CharacterSkeletonSegment& segment)
{
    CharacterSkeletonLimbSegmentData frame;

    Node* nodeA = segment.nodes_[0];
    Node* nodeB = segment.nodes_[1];
    Node* nodeC = segment.nodes_[2];

    const Matrix3x4 transformA = nodeA->GetWorldTransform();
    const Matrix3x4 transformB = nodeB->GetWorldTransform();
    const Matrix3x4 transformC = nodeC->GetWorldTransform();

    const Matrix3x4& initialA = segment.initialTransforms_[0];
    const Matrix3x4& initialB = segment.initialTransforms_[1];
    const Matrix3x4& initialC = segment.initialTransforms_[2];

    // Get relative position
    frame.position_ = transformC.Translation() - transformA.Translation();

    // Get joint bending direction
    const Vector3 positionAproj = ProjectPointOntoSegment(transformB.Translation(), transformA.Translation(), transformC.Translation());
    const Vector3 initialDirection = initialC.Translation() - initialA.Translation();
    const Vector3 newDirection = transformC.Translation() - transformA.Translation();
    frame.direction_ = Quaternion(newDirection, initialDirection) * (transformB.Translation() - positionAproj).Normalized();

    // Revert transforms to initial
    const Vector3 zeroDelta = nodeA->GetWorldPosition() - initialA.Translation();
    if (!MatchChildPosition(*nodeA, *nodeB, initialB.Translation() + zeroDelta))
        URHO3D_LOGWARNING("Failed to resolve thigh-calf segment of foot animation");
    if (!MatchChildPosition(*nodeB, *nodeC, initialC.Translation() + zeroDelta))
        URHO3D_LOGWARNING("Failed to resolve calf-heel segment of foot animation");

    // Gather rotations
    const Quaternion pureRotationA = nodeA->GetWorldRotation() * initialA.Rotation().Inverse();
    const Quaternion pureRotationB = nodeB->GetWorldRotation() * initialB.Rotation().Inverse();
    const Quaternion pureRotationC = nodeC->GetWorldRotation() * initialC.Rotation().Inverse();

    frame.rotationA_ = GetAngle(pureRotationA);
    frame.rotationB_ = GetAngle(pureRotationB);
    frame.rotationC_ = transformC.Rotation() * initialC.Rotation().Inverse();

    // Flip angles if needed
    const Vector3 directionAB = (nodeB->GetWorldPosition() - nodeA->GetWorldPosition()).Normalized();
    const Vector3 directionBC = (nodeC->GetWorldPosition() - nodeB->GetWorldPosition()).Normalized();
    const Vector3 axisAB = GetAxis(pureRotationA);
    const Vector3 axisBC = GetAxis(pureRotationB);

    frame.rotationA_ *= Sign(directionAB.DotProduct(axisAB));
    frame.rotationB_ *= Sign(directionBC.DotProduct(axisBC));

    track_.Push(frame);
}

void LimbAnimationTrack::MergeFrame(unsigned firstFrame, unsigned secondFrame, float factor,
    float weight, CharacterSkeletonSegmentData& dest)
{
    const CharacterSkeletonLimbSegmentData& first = track_[firstFrame];
    const CharacterSkeletonLimbSegmentData& second = track_[secondFrame];
    dest.Merge(first, weight * (1 - factor));
    dest.Merge(second, weight * factor);
}

bool LimbAnimationTrack::SaveXML(XMLElement& dest) const
{
    const unsigned numKeys = track_.Size();
    for (unsigned i = 0; i < numKeys; ++i)
    {
        XMLElement child = dest.CreateChild("frame");
        child.SetVector3("position", track_[i].position_);
        child.SetVector3("direction", track_[i].direction_);
        child.SetFloat("rotation0", track_[i].rotationA_);
        child.SetFloat("rotation1", track_[i].rotationB_);
        child.SetQuaternion("rotation2", track_[i].rotationC_);
    }
    return true;
}

bool LimbAnimationTrack::LoadXML(const XMLElement& source)
{
    for (XMLElement child = source.GetChild("frame"); child; child = child.GetNext())
    {
        CharacterSkeletonLimbSegmentData frame;
        frame.position_ = child.GetVector3("position");
        frame.direction_ = child.GetVector3("direction");
        frame.rotationA_ = child.GetFloat("rotation0");
        frame.rotationB_ = child.GetFloat("rotation1");
        frame.rotationC_ = child.GetQuaternion("rotation2");
        track_.Push(frame);
    }
    return true;
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
    // Read timestamps
    if (XMLElement timestampsNode = source.GetChild("timestamps"))
        for (XMLElement child = timestampsNode.GetChild("frame"); child; child = child.GetNext())
            timeStamps_.Push(child.GetFloat("time"));

    // Read tracks
    for (XMLElement child = source.GetChild(); child; child = child.GetNext())
        if (child.GetName() != "timestamps")
        {
            const CharacterSkeletonSegmentType type = GetCharacterSkeletonSegmentType(child.GetName());
            if (SharedPtr<CharacterAnimationTrack> track = CharacterAnimationTrack::Create(type, child.GetAttribute("name")))
            {
                track->LoadXML(child);
                tracks_.Push(track);
            }
        }
    return true;
}

bool CharacterAnimation::SaveXML(XMLElement& dest) const
{
    unsigned numFrames = timeStamps_.Size();

    // Save timestamps
    if (XMLElement node = dest.CreateChild("timestamps"))
        for (unsigned i = 0; i < numFrames; ++i)
            node.CreateChild("frame").SetFloat("time", timeStamps_[i]);

    // Save tracks
    for (CharacterAnimationTrack* track : tracks_)
    {
        XMLElement node = dest.CreateChild(track->GetTypeString());
        node.SetAttribute("name", track->GetName());
        track->SaveXML(node);
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

void CharacterAnimation::GetKeyFrameIndex(float time, unsigned& index) const
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

unsigned CharacterAnimation::GetKeyFrameIndex(float time) const
{
    unsigned index = 0;
    GetKeyFrameIndex(time, index);
    return index;
}

void CharacterAnimation::GetKeyFrame(float time, unsigned& firstFrame, unsigned& secondFrame, float& factor) const
{
    firstFrame = GetKeyFrameIndex(time);
    secondFrame = firstFrame + 1 < timeStamps_.Size() ? firstFrame + 1 : 0;
    const float timeInterval = Max(0.0f, timeStamps_[secondFrame] - timeStamps_[firstFrame]);
    factor = timeInterval > 0.0f ? (time - timeStamps_[firstFrame]) / timeInterval : 1.0f;
}

float CharacterAnimation::GetLength() const
{
    return timeStamps_.Size() > 0 ? timeStamps_.Back() - timeStamps_.Front() : 0.0f;
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
    unsigned numFrames = 0;
    for (const auto& item : animation.GetTracks())
    {
        const AnimationTrack& track = item.second_;
        if (numFrames != 0 && track.keyFrames_.Size() != numFrames)
        {
            URHO3D_LOGERROR("CharacterAnimation has inconsistent number of key frames in tracks");
            return false;
        }

        numFrames = track.keyFrames_.Size();
        if (numFrames > 0)
        {
            timeStamps_.Resize(numFrames);
            for (unsigned i = 0; i < track.GetNumKeyFrames(); ++i)
                timeStamps_[i] = track.keyFrames_[i].time_;
        }
    }
    if (numFrames == 0)
    {
        URHO3D_LOGWARNING("CharacterAnimation is empty");
        return true;
    }

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
        tracks_.Push(CharacterAnimationTrack::Create(segment.type_, segment.name_));

    // Gather segments data
    Vector<CharacterSkeletonSegment> segmentData;
    rig.AllocateSegmentData(segmentData, skeleton, transform);

    // Collect tracks
    for (unsigned i = 0; i < numFrames; ++i)
    {
        // Play animation
        animationState->SetTime(timeStamps_[i]);
        animationState->Apply();
        node.MarkDirty();

        // Read tracks
        for (unsigned j = 0; j < segmentData.Size(); ++j)
            if (tracks_[j])
                tracks_[j]->ImportFrame(segmentData[j]);
        // #TODO Backup
    }

    return true;
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
    context->RegisterFactory<CharacterAnimationController>();
    URHO3D_COPY_BASE_ATTRIBUTES(AnimationController);
    URHO3D_MIXED_ACCESSOR_ATTRIBUTE("Skeleton", GetSkeletonAttr, SetSkeletonAttr, ResourceRef, ResourceRef(XMLFile::GetTypeStatic()), AM_DEFAULT);
}

void CharacterAnimationController::SetAnimationTransform(const Matrix3x4& transform)
{
    animationTransform_ = transform;

    // #TODO Add dirty flag
    animatedModelSkeleton_ = nullptr;
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
    UpdateHierarchy();
    AnimationController::Update(timeStep);
    if (animatedModel_)
    {
        animatedModel_->ApplyAnimation();
        ApplyAnimation();
    }
}

void CharacterAnimationController::ApplyAnimation()
{
    // #TODO Remove it
    animatedModel_->GetSkeleton().Reset();

    // Reset segments data
    for (CharacterSkeletonSegment& segment : segmentData_)
        segment.data_->Reset();

    // Animate segments
    for (const AnimationControl& animationControl : GetAnimations())
    {
        CharacterAnimation* characterAnimation = GetCharacterAnimation(animationControl.name_);
        if (!characterAnimation)
            continue;
        AnimationState* animationState = GetAnimationState(animationControl.name_);
        if (!animationState)
            continue;

        unsigned firstFrame;
        unsigned secondFrame;
        float factor;
        characterAnimation->GetKeyFrame(animationState->GetTime(), firstFrame, secondFrame, factor);

        for (CharacterSkeletonSegment& segment : segmentData_)
            if (CharacterAnimationTrack* track = characterAnimation->FindTrack(segment.name_))
                track->MergeFrame(firstFrame, secondFrame, factor, animationState->GetWeight(), *segment.data_);
    }

    // Apply animations
    for (CharacterSkeletonSegment& segment : segmentData_)
        segment.data_->Apply(node_->GetWorldTransform(), segment);
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

void CharacterAnimationController::UpdateHierarchy()
{
    // Get component
    if (!node_ || !skeleton_)
        return;
    if (animatedModel_ && &animatedModel_->GetSkeleton() == animatedModelSkeleton_)
        return;
    animatedModel_ = node_->GetComponent<AnimatedModel>();
    animatedModelSkeleton_ = &animatedModel_->GetSkeleton();
    if (!animatedModel_)
        return;

    skeleton_->AllocateSegmentData(segmentData_, *animatedModelSkeleton_, animationTransform_);
}

CharacterAnimation* CharacterAnimationController::GetCharacterAnimation(const String& animationName)
{
    if (animationCache_.Contains(animationName))
        return animationCache_[animationName];

    ResourceCache* cache = GetSubsystem<ResourceCache>();
    CharacterAnimation* characterAnimation = cache->GetResource<CharacterAnimation>(animationName.Replaced(".ani", ".xml", false));
    animationCache_[animationName] = characterAnimation;
    return characterAnimation;
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
    const Matrix3x4 animationToWorld = node_->GetWorldTransform() * animationTransform_;
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

    RegisterComponent<CharacterAnimationController>(engine, "CharacterAnimationController");
    RegisterSubclass<CharacterAnimationController, AnimationController>(engine, "AnimationController", "CharacterAnimationController");
    engine->RegisterObjectMethod("CharacterAnimationController", "void SetAnimationTransform(const Matrix3x4&in)", asMETHOD(CharacterAnimationController, SetAnimationTransform), asCALL_THISCALL);
    engine->RegisterObjectMethod("CharacterAnimationController", "void SetTargetTransform(const String&in, const Matrix3x4&in)", asFUNCTION(CharacterAnimationController_SetTargetTransform), asCALL_CDECL_OBJLAST);
    engine->RegisterObjectMethod("CharacterAnimationController", "void SetTargetRotationAmount(const String&in, float)", asFUNCTION(CharacterAnimationController_SetTargetRotationAmount), asCALL_CDECL_OBJLAST);
    engine->RegisterObjectMethod("CharacterAnimationController", "void SetTargetRotationBalance(const String&in, float)", asFUNCTION(CharacterAnimationController_SetTargetRotationBalance), asCALL_CDECL_OBJLAST);
}

}

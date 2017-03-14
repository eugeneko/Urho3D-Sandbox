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

//////////////////////////////////////////////////////////////////////////
CharacterSkeleton::CharacterSkeleton(Context* context)
    : Resource(context)
{
}

CharacterSkeleton::~CharacterSkeleton()
{
}

void CharacterSkeleton::RegisterObject(Context* context)
{
    context->RegisterFactory<CharacterSkeleton>();
}

bool CharacterSkeleton::BeginLoad(Deserializer& source)
{
    SharedPtr<XMLFile> xmlFile = MakeShared<XMLFile>(context_);
    if (xmlFile->Load(source))
        return BeginLoad(xmlFile->GetRoot());
    return false;
}

bool CharacterSkeleton::BeginLoad(const XMLElement& source)
{
    modelName_ = source.GetChild("model").GetAttribute("name");
    if (modelName_.Empty())
    {
        URHO3D_LOGERROR("CharacterSkeleton model name mustn't be empty");
        return false;
    }

    for (XMLElement segmentNode = source.GetChild("segment2"); !segmentNode.IsNull(); segmentNode = segmentNode.GetNext("segment2"))
    {
        CharacterSkeletonSegment2 segment;
        segment.name_ = segmentNode.GetAttribute("name");
        segment.rootBone_ = segmentNode.GetAttribute("root");
        segment.jointBone_ = segmentNode.GetAttribute("joint");
        segment.targetBone_ = segmentNode.GetAttribute("target");

        if (segment.name_.Empty())
        {
            URHO3D_LOGERROR("CharacterSkeleton 2-segment name mustn't be empty");
            return false;
        }

        if (segment.rootBone_.Empty() || segment.jointBone_.Empty() || segment.targetBone_.Empty())
        {
            URHO3D_LOGERROR("CharacterSkeleton 2-segment bones names mustn't be empty");
            return false;
        }

        segments2_.Populate(segment.name_, segment);
    }
    return true;
}

//////////////////////////////////////////////////////////////////////////
CharacterAnimation::CharacterAnimation(Context* context)
    : Resource(context)
{
}

CharacterAnimation::~CharacterAnimation()
{
}

void CharacterAnimation::RegisterObject(Context* context)
{
    context->RegisterFactory<CharacterAnimation>();
}

bool CharacterAnimation::BeginLoad(Deserializer& source)
{
    SharedPtr<XMLFile> xmlFile = MakeShared<XMLFile>(context_);
    if (xmlFile->Load(source))
        return Load(xmlFile->GetRoot());
    return false;
}

bool CharacterAnimation::Load(const XMLElement& source)
{
    for (XMLElement segmentNode = source.GetChild("segment2"); !segmentNode.IsNull(); segmentNode = segmentNode.GetNext("segment2"))
    {
        CharacterAnimationSegment2Track track;
        track.name_ = segmentNode.GetAttribute("name");
        track.initialDirection_ = segmentNode.GetVector3("baseDirection");

        for (XMLElement keyFrameNode = segmentNode.GetChild("keyFrame");
            !keyFrameNode.IsNull();
            keyFrameNode = keyFrameNode.GetNext("keyFrame"))
        {
            CharacterAnimationSegment2KeyFrame keyFrame;
            keyFrame.time_ = keyFrameNode.GetFloat("time");
            keyFrame.heelPosition_ = keyFrameNode.GetVector3("targetPosition");
            keyFrame.kneeDirection_ = keyFrameNode.GetVector3("jointOrientation");
            keyFrame.thighRotationFix_ = keyFrameNode.GetQuaternion("rootRotation");
            keyFrame.calfRotationFix_ = keyFrameNode.GetQuaternion("jointRotation");
            keyFrame.heelRotationLocal_ = keyFrameNode.GetQuaternion("targetRotation");
            keyFrame.heelRotationWorld_ = keyFrameNode.GetQuaternion("targetRotationWorld");
            track.keyFrames_.Push(keyFrame);
        }
        segments2_.Populate(track.name_, track);
    }
    return true;
}

bool CharacterAnimation::Save(XMLElement& dest) const
{
    for (const Segment2TrackMap::KeyValue& elem : segments2_)
    {
        XMLElement segment2Node = dest.CreateChild("segment2");
        const CharacterAnimationSegment2Track& track = elem.second_;
        segment2Node.SetAttribute("name", track.name_);
        segment2Node.SetVector3("baseDirection", track.initialDirection_);
        for (const CharacterAnimationSegment2KeyFrame& keyFrame : track.keyFrames_)
        {
            XMLElement keyFrameNode = segment2Node.CreateChild("keyFrame");
            keyFrameNode.SetFloat("time", keyFrame.time_);
            keyFrameNode.SetVector3("targetPosition", keyFrame.heelPosition_);
            keyFrameNode.SetVector3("jointOrientation", keyFrame.kneeDirection_);
            keyFrameNode.SetQuaternion("rootRotation", keyFrame.thighRotationFix_);
            keyFrameNode.SetQuaternion("jointRotation", keyFrame.calfRotationFix_);
            keyFrameNode.SetQuaternion("targetRotation", keyFrame.heelRotationLocal_);
            keyFrameNode.SetQuaternion("targetRotationWorld", keyFrame.heelRotationWorld_);
        }
    }
    return true;
}

bool CharacterAnimation::Save(Serializer& dest) const
{
    SharedPtr<XMLFile> xml(new XMLFile(context_));
    XMLElement animationNode = xml->CreateRoot("animation");

    Save(animationNode);
    return xml->Save(dest);
}

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
    AnimationController::Update(timeStep);
    AnimatedModel* animatedModel = node_->GetComponent<AnimatedModel>();
    animatedModel->ApplyAnimation();
    ApplyAnimation(animatedModel);
}

void CharacterAnimationController::ApplyAnimation()
{
    AnimatedModel* animatedModel = node_->GetComponent<AnimatedModel>();
    ApplyAnimation(animatedModel);
}

void CharacterAnimationController::ApplyAnimation(AnimatedModel* animatedModel)
{
    if (skeleton_ && !skeleton_->GetSegments2().Empty())
    {
        for (const CharacterSkeleton::Segment2Map::KeyValue& elem : skeleton_->GetSegments2())
            UpdateSegment2(animatedModel, elem.second_);
    }
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

CharacterAnimation* CharacterAnimationController::GetCharacterAnimation(const String& animationName)
{
    if (animationCache_.Contains(animationName))
        return animationCache_[animationName];

    ResourceCache* cache = GetSubsystem<ResourceCache>();
    CharacterAnimation* characterAnimation = cache->GetResource<CharacterAnimation>(animationName.Replaced(".ani", ".xml", false));
    animationCache_[animationName] = characterAnimation;
    return characterAnimation;
}

void CharacterAnimationController::UpdateSegment2(AnimatedModel* animatedModel, const CharacterSkeletonSegment2& segment)
{
    // Get nodes and bones
    Skeleton& skeleton = animatedModel->GetSkeleton();
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
        CharacterAnimationSegment2Track* track = characterAnimation->FindTrack(segment.name_);
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

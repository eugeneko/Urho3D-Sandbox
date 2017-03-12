project (Player)

set (TARGET_NAME Player)
set (SOURCE_FILES Main.cpp Urho3DPlayer.h Urho3DPlayer.cpp)
setup_main_executable ()
target_link_libraries (Player CharacterAnimator)

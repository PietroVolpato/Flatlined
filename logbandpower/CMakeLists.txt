cmake_minimum_required(VERSION 3.0.2)
project(logbandpower)

find_package(catkin REQUIRED COMPONENTS
  message_generation
  rospy
  std_msgs
)

catkin_package(
 CATKIN_DEPENDS message_runtime rospy std_msgs
)

include_directories(
  ${catkin_INCLUDE_DIRS}
)

catkin_install_python(PROGRAMS
    scripts/bandpower.py
    scripts/thresholding.py 
    scripts/recorder.py 
  DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)
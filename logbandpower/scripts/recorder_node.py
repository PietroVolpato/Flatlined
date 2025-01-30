import rospy
import rosbag
from rosneuro_msgs.msg import NeuroFrame, NeuroEvent

##########################################################

class RecorderNode:
    def __init__(self):
        rospy.init_node('recorder_node', anonymous=True)

            # Obtain parameters
        self.filepath = rospy.get_param('~filepath', '$(find logbandpower)/results/')
        
            # Create the bag file
        self.bag = rosbag.Bag(self.filepath + 'recorded_data.bag', 'w')
        rospy.loginfo("Recorder: starting the bag file")

            # Subscribers
        self.filtered_sub = rospy.Subscriber('/eeg/filtered', NeuroFrame, self.filtered_callback)
        self.bandpower_sub = rospy.Subscriber('/eeg/bandpower', NeuroFrame, self.bandpower_callback)
        self.events_sub = rospy.Subscriber('/events/bus', NeuroEvent, self.events_callback)

        rospy.on_shutdown(self.shutdown_hook)

    #----------------------------------------------------#

    def filtered_callback(self, msg):
        self.bag.write('/eeg/filtered', msg)

    def bandpower_callback(self, msg):
        self.bag.write('/eeg/bandpower', msg)

    def events_callback(self, msg):
         self.bag.write('/events/bus', msg)

    def shutdown_hook(self):
      self.bag.close()
      rospy.loginfo("Recorder: closing the bag file")
      
    def run(self):
        rospy.spin()

##########################################################

if __name__ == '__main__':
    try:
        recorder = RecorderNode()
        recorder.run()
    except rospy.ROSInterruptException:
        pass
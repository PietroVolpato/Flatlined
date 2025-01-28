import rospy
from rosneuro_msgs.msg import NeuroFrame
from rosneuro_msgs.msg import NeuroEvent

##########################################################

class ThresholdingNode:
    def __init__(self):
        rospy.init_node('thresholding_node')

            # Obtain parameters
        self.selected_channel = rospy.get_param('~selected_channel', 0)
        self.threshold = rospy.get_param('~threshold', 10.0)
        rospy.loginfo("Thresholding node called on channel #%d with threshold %f", self.selected_channel, self.threshold)

        self.threshold_crossed = False

            # Subscriber and Publisher
        self.bandpower_sub = rospy.Subscriber("/eeg/bandpower", NeuroFrame, self.bandpower_callback)
        self.event_pub = rospy.Publisher("/events/bus", NeuroEvent, queue_size=10)

    #----------------------------------------------------#

    def bandpower_callback(self, msg):
        
        if len(msg.eeg.data) <= self.selected_channel:
            rospy.logerr("\t-> Selected channel index out of bounds!")
            return
        
        current_value = msg.eeg.data[self.selected_channel]

        if current_value > self.threshold and not self.threshold_crossed:
            self.publish_event(self.selected_channel)
            self.threshold_crossed = True
        elif current_value < self.threshold:
            self.threshold_crossed = False

    #----------------------------------------------------#
    
    def publish_event(self, channel):
        event_msg = NeuroEvent()
        event_msg.family = channel
        event_msg.description = 'Crossed the threshold'
        self.event_pub.publish(event_msg)
        rospy.loginfo("Threshold crossed on channel #%d ---> Event published!", channel)

##########################################################


def main ():
    thresholding_node = ThresholdingNode()
    rospy.spin()

if __name__ == '__main__':
    main()
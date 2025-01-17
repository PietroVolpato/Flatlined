import rospy
from std_msgs.msgs import NeuroFrame
from rosnuero_msgs import NeuroEvent

def main ():
    rospy.init_node('thresholding')

    rospy.Subscriber('/eeg/bandpower', NeuroFrame, callback)

    pub = rospy.Publisher('/event/bus', NeuroEvent, queue_size=1)
import rospy
import numpy as np
from rosneuro_msgs.msg import NeuroFrame
from std_msgs.msg import Float32MultiArray

class BandpowerNode:
    def __init__(self):
        # Initialize ROS node
        rospy.init_node('bandpower', anonymous=True)

        # Parameters
        self.rate = int(rospy.get_param('rate', 16))  # Publishing rate

        # Publisher and Subscriber
        self.pub = rospy.Publisher('/eeg/bandpower', NeuroFrame, queue_size=1)
        self.sub = rospy.Subscriber('/eeg/filtered', NeuroFrame, self.callback)

        # Buffer and state variables
        self.buffer = None
        self.seq = 0
        self.new_data = False
        self.current_frame = None

        # ROS rate
        self.loop_rate = rospy.Rate(self.rate)

    def callback(self, data: NeuroFrame):
        """Callback function to receive NeuroFrame messages."""
        self.current_frame = data
        self.new_data = True

    def initialize_buffer(self, nchannels, buffer_size):
        """Initialize the buffer for the first time."""
        self.buffer = np.zeros((nchannels, buffer_size))

    def update_buffer(self, data):
        """Update the buffer with new EEG data."""
        eeg_data = np.array(data.eeg.data).reshape((data.eeg.info.nchannels, data.eeg.info.nsamples))

        # Remove oldest samples and append new data
        self.buffer = np.hstack((self.buffer[:, data.eeg.info.nsamples:], eeg_data))

    def compute_bandpower(self):
        """Compute the bandpower for each channel in the buffer."""
        # Calculate power and apply log10
        power = np.multiply(self.buffer, self.buffer)
        log_power = np.log10(power)

        return np.mean(log_power, axis=1)

    @staticmethod
    def generate_neuroframe_message(bandpower_values, old_message):
        msg = NeuroFrame()
        msg.eeg.info = old_message.eeg.info
        msg.eeg.data = bandpower_values.flatten().tolist()
    
        return msg

    def run(self):
        """Main loop of the ROS node."""
        while not rospy.is_shutdown():
            if self.new_data:
                # Process new data
                data = self.current_frame

                # Initialize buffer if needed
                if self.seq == 0:
                    buffer_size = data.sr  # Assuming buffer size matches sampling rate
                    self.initialize_buffer(data.eeg.info.nchannels, buffer_size)

                # Update buffer
                self.update_buffer(data)
                self.seq += 1

                # Compute bandpower if buffer is sufficiently filled
                if self.seq * data.eeg.info.nsamples >= data.sr:
                   bandpower_values = self.compute_bandpower()
                   msg = self.generate_neuroframe_message(bandpower_values, self.current_frame)
                   self.pub.publish(msg)

                self.new_data = False

            self.loop_rate.sleep()

if __name__ == '__main__':
    try:
        node = BandpowerNode()
        node.run()
    except rospy.ROSInterruptException:
        pass

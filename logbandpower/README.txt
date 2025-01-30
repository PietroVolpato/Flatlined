Member contribution:
Gnocato Margherita e Lo Faro Alessio: remapping of the topic and fixing of the acquisition and filterchain node
Kovachev Zlatko: creation of the thresholding node and recording node
Volpato Pietro: creation of the bandpower node

HOW TO RUN: 
Open a terminal, build the package logbandpower and run roslaunch logbandpower logbandpower.launch
Note that since the data must not be included, it is necessary to change the path of the data in the acquisition node on the launch file in the correspondence of the "devarg" parameter.

WHAT WE HAVE DONE:
we decide to use the acquisition, filterchain nodes already developed in the ROsNeuro package while we opted for a recasting of the bandpower and thresholding nodes.
Thw working of this last two nodes is the same as in the RosNeuro package but they have been modified for making them suitable for the assignment.
We decide not to develop a node for the recording of the filtered data but to use the rosbag package that can be launched directly from the launch file.
All the values needed by the nodes can be adjusted from the launch file by changing the parameters passed by the arguments.
The recorded results are saved in the results folder
Member contribution:
    -Gnocato Margherita e Lo Faro Alessio: remapping of the topic and fixing of the acquisition and filterchain node
    
    -Kovachev Zlatko: creation of the thresholding node and recorder node
    
    -Volpato Pietro: creation of the bandpower node


HOW TO RUN: 
    - Open a terminal, build the package logbandpower and run "roslaunch logbandpower logbandpower.launch"
    - Note that since the data must not be included, it is necessary to change the path of the data in the acquisition node on the launch file in the correspondence of the "devarg" parameter.

WHAT WE HAVE DONE:

We opted to utilize the acquisition and filterchain nodes available in the RosNeuro package while we recreated the bandpower and thresholding nodes. 
Although these two nodes function similarly to those in the RosNeuro package, they were modified to suit our assignment. 
Additionally, we developed a node for data recording using the rosbag package. Node parameters can be modified through the launch file by adjusting the argument-passed parameters. 
The recorded data is saved in the results folder.
# Reversible Cellular Automata Simulator

### Synopsis

This software project is designed to simulate reversible cellular automata \(RCA\), calculate particular data for each of these simulations, and optionally save the resuting data and the corresponding image of the system. The software itself is currently in the form of an OS X application. The systems simulated allow for most parameters to be adjusted by the user, and the width of the simulation space is periodic, meaning the indices wrap around.

This project was adapted into an OS X application using the [Qt application development framework](https://wiki.qt.io/About_Qt). 

### Installation

Most simply, this software can be installed to a local machine \(running Mac OS X 10.7 Lion or later\) by simply selecting the **Download ZIP** option at the top of this repository page and then downloading the files into any directory. From there, the simulation can be started with a simple double-click on the **RCA_Simulator.app** OS X application file.

To recompile the code, the provided makefile can be used. Therefore, after the project code is edited \(likely in RCASimulator.cpp\), the application can be recreated by simply typing `make` from the Mac Terminal after navigating to the project directory. Because this application utilizes the Qt framework, recompilation of the source code will require that [Qt is downloaded](https://www.qt.io/download/) to the local machine. Additionally, [GMP](https://gmplib.org/), the GNU multiple precision arithmetic library, will need to be installed.

### Usage

After opening the application \(RCA\_Simulator.app\) by either double-clicking its icon or right-clicking it and selecting _Open_ from the pop-up menu, a new input window should appear. This input window contains various fields of information used to set up the particular RCA system parameters that the user wants to run. For example, the user specifies a space width and the number of time steps to simulate. 
After clicking the accept button at the bottom of the window, the simulation will run in the background. Two new windows will pop up. The background window is a "log window" of information about the simulation as it runs, displaying calculations and status messages. The foreground window is a simple message box that asks if the user would like to save the calculated data and the image of the RCA that was just simulated. 

If the user chooses to save the results of the simulation \(by selecting _yes_ in the previous window\), another window will pop up after a few seconds. This new window will show a directory hierarchy for the local machine. The user can traverse this directory hierarchy to specify a particular existing directory in which to save the resulting simulation data/image.

### Credit

This software project utilized development tools/libraries from the following sources:

1. [Qt Application Development Framework](https://www.qt.io/)
2. [GMP Multiple Precision Arithmetic Library](https://gmplib.org/)
3. [Andrew Noske's Qt Custom Input Dialog](http://www.andrewnoske.com/wiki/Code_-_qt_custom_input_dialog#mw-navigation)
4. [CImg Template Image Processing Toolkit](http://cimg.eu/)

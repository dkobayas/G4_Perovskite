# G4_Perovskite
Geant4 simulation package for development of perovskite semiconductor radiation detector

## == Quick Start Guide ==

#### First of all, please clone this repository via git clone command;

  `$ git clone https://github.com/dkobayas/G4_Perovskite.git`

  note: you have to be registered as the github user.

#### Then, you can find the base directory, G4_Perovskite.

  Please move into the directory.

  `$ cd G4Pvsk`

  You can run the cmake command, and find the Makefile generated with your environment.

  `$ cmake .`

#### Now, you can compile this package.
  
  Please do make
  
  `$ make`
  
  You can find and execute the executable file, G4Pvsk.
  
  Finally, please execute with this run macro.
  
  `$ ./G4Pvsk run1.mac`
  
  I hope you can find the processing lines, and the result output as root file, Perovskite.root.
  
  note: You have to setup ROOT for processing with output as root file.

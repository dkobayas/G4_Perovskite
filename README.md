# G4_Perovskite
Geant4 simulation package for development of perovskite semiconductor radiation detector

## == Quick Start Guide ==

#### First of all, please clone this repository via git clone command;

  `$ git clone https://github.com/dkobayas/G4_Perovskite.git`

  note: you have to be registered as the github user.

#### Then, you can find the base directory, G4_Perovskite.

  Let's move into this directory and make a new directory for built process.

  `$ cd G4_Perovskite`

  `$ mkdir built`
  
  `$ cd built`

  Plase move to the built directory, and copy Mkaefile from G4Pvsk to built

  $ cp ../G4Pvsk/CMakeList.txt .

#### Now, you can compile this package.
  
  `$ cmake .`
  
  Then, you can find generated Makefile built for your environment.
  
  Please do make
  
  `$ make`
  
  You can find and execute the executable file, G4Pvsk.
  
  Please copy run1.mac(one of run macros) from G4Pvsk to built.
  
  `$ cp ../G4Pvsk/run1.mac .`
  
  Finally, please execute with this run macro.
  
  `$ ./G4Pvsk run1.mac`
  
  I hope you can find the processing lines, and the result output as root file, Perovskite.root.
  
  note: You have to setup ROOT for processing with output as root file.

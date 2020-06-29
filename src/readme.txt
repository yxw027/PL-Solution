RTKDEMO ver.0.1  2006/01/29 by T.TAKASU (ttaka@gpspp.sakura.ne.jp)

1.Description

This is an implementation of demonstration for the real-time kinematic (RTK)
positioning algorithm with GPS/GNSS. In the algorithm, parameter adjustment
is done by the Extended Kalman-Filter using double-differenced phase
observables. Integer ambiguity is resolved by LAMBDA/MLAMBDA method.
Cycle-slip/outlier detection is not implemented. The reference satellite is
fixed to the highest elevation one at the first epoch.


2.Files

testrtk.m : test driver for RTKDEMO
rtkdemo.m : real-time kinematic (RTK) positioning demo.
mlambda.m : LAMBDA/MLAMBDA integer ambiguity resolution
readrnx.m : read rinex observation data/navigation message files


3.Environment

Matlab(5.x,6.x,7.x) or Octave(2.1.x) + Octave-forge(matlab-compatible libs)


4.More Information

http://gpspp.sakura.ne.jp


5.History

2006/01/29 ver.0.1  new release


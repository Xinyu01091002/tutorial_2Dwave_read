# Basilisk 2D wave by Xinyu
This is a test code in Basilisk C (http://basilisk.fr/Front%20Page) to simulate a two-phase surfave water wave breaking with AMR. 
Currently goals are 
- [x] Initialize the phase and velocity field with surface elevation and surface velocity potential.
- [x] Adjust the vertical resolution by modifying the metric factor fm, cm.
- [x] Compare against the third order stoke wave cases by Popinet (http://basilisk.fr/sandbox/popinet/wave.c).
- [ ] Kick off the simulation for a wave group
- [ ] Compare against the Kinematics with OceanWave3D

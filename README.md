# CAPRICEP
New TSP

All components are uploaded. Type "realTimeTester" to start the application. Supporting materials will be ready soon. (Wed Nov 11 11:31:20 JST 2020)

## Reference

Kawahara, H. and Yatabe, K.: Cascaded all-pass filterswith randomized center frequencies and phase polarity for　acoustic and speech measurement and data augmentation,arXiv:2010.13185(2020).

### Abstract

We introduce a new member of TSP (Time Stretched Pulse) for acoustic and speech measurement infrastructure, 
based on a simple all-pass filter and systematic randomization. 
This new infrastructure fundamentally upgrades our previous measurement procedure, 
which enables simultaneous measurement of multiple attributes, including non-linear ones
without requiring extra filtering nor post-processing. 
Our new proposal establishes a theoretically solid, flexible, and extensible foundation in acoustic measurement.
Moreover, it is general enough to provide versatile research tools for other fields, such as biological signal analysis.
We illustrate using acoustic measurements and data augmentation as representative examples among various prospective applications. 
We open-sourced MATLAB implementation. 
It consists of an interactive and real-time acoustic tool, MATLAB functions, and supporting materials.

## Interactive and real-time acoustic measurement tool
Currently the following links to an FVN-based tool. Will be updated to CAPRICEP-based one today (05.Nov.2020).
[![Watch the video](fvnAcousticTool.jpg)](https://youtu.be/yT7N1MITKzY)


## Sample files

The following files sound pervcptually indistinguishable (at least for me).

Please refer to Fig.5 of the reference.

 ### The worst SNR (-1dB) file filtered by CAPRICEP

capricepWorstSNRminus1dB.wav

### The best SNR (6dB) file filtered by CAPRICEP

capricepBestSNR6dB.wav

### The original file (level adjusted to CAPRICEP samples)

originalLevelAdjusted.wav

## Memo (working)

### All components are uploaded. Type "realTimeTester" to start the application
### Added off-line analysis function
### Uploaded unit-CAPRICEPs for 100, 200, 400, 800, and 1600 ms

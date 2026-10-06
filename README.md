# DSurfTomo-4D

DSurfTomo-4D is an extended version of the direct surface-wave tomography framework DSurfTomo for time-lapse monitoring of subsurface velocity variations.
The original DSurfTomo method directly inverts surface-wave dispersion measurements for three-dimensional shear-wave velocity structures without constructing intermediate phase-velocity maps.
In this work, the original framework is further developed to enable time-lapse imaging through common-path differential surface-wave tomography.
The DSurfTomo-4D framework constructs differential travel-time data from common observations between two epochs and directly inverts for relative velocity variations, aiming to reduce spurious temporal signals associated with inconsistent observation geometries.

## Relationship with original DSurfTomo

This repository is developed based on the original DSurfTomo code released by Fang et al. (2015). 

The present version introduces additional modules for time-lapse inversion and differential velocity-change recovery.
The original DSurfTomo code and Taipei Basin dataset are available at:
https://github.com/HongjianFang/DSurfTomo.git

Users interested in the fundamental tomography methodology should refer to:
Fang, H., Yao, H., Zhang, H., Huang, Y. C., & van der Hilst, R. D. (2015). Direct inversion of surface wave dispersion for three-dimensional shallow crustal structure based on ray tracing: methodology and application. Geophysical Journal International, 201(3), 1251-1263. 

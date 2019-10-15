# modelling-and-optimizing-the-geometry-of-acoustic-diffusers-using-FDTD
This repository contains the MATLAB files for modelling and optimizing the geometry of acoustic diffusers using Finite-Difference Time-Domain method, as part of the outcomes of my Master's dissertation

The .m files should be grouped and used together according to their functions:
"Yuqing_dc_calculation.m" - calls "Yuqing_FDTD_func.m" - calls "Yuqing_qrs.m"
They together perform room acoustics simulation. You can simulate and visualize sound propagation and diffusion, or measure the diffusion coefficient for a known diffuser with these files.

"Yuqing_opti.m" - calls "Yuqing_opti_dc_calculation_func.m" - calls "Yuqing_opti_FDTD_func.m" - calls "Yuqing_qrs.m"
They together perform diffuser geometry optimization. With these files, you can optimize the geometry of a diffuser with selected number of wells for a desired bandwidth.

There is also an animation of sound propagation and diffusion produced by the room acoustics model in .mp4.

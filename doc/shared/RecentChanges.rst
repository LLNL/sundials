**Major Features**

**New Features and Enhancements**

Users may now disable interpolated output in ARKODE steppers by passing
``ARK_INTERP_NONE`` to the ``SetInterpolantType`` functions:

* :c:func:`ARKStepSetInterpolantType`,
* :c:func:`ERKStepSetInterpolantType`,
* :c:func:`MRIStepSetInterpolantType`, and
* :c:func:`SPRKStepSetInterpolantType`.

When interpolation is disabled, rootfinding is not supported, implicit methods
must use the trivial predictor (the default option), and interpolation at stop
times cannot be used (interpolating at stop times is disabled by default). With
interpolation disabled, ``Evolve`` functions called in ``ARK_NORMAL`` mode will
return at or past the requested output time (setting a stop time may still be
used to halt the integrator at a specific time). Disabling interpolation will
reduce the memory footprint of an integrator by two or more state vectors
(depending on the interpolant type and degree) which can be beneficial when
interpolation is not needed e.g., when integrating to a final time without
output in between or using an explicit fast time scale integrator with MRIStep.

**Bug Fixes**

Updated the CMake variable ``HIP_PLATFORM`` default to ``amd`` as the previous
default, ``hcc``, is no longer recognized in ROCm 5.7.0 or newer. The new
default is also valid in older version of ROCm (at least back to version 4.3.1).

Fixed a bug in the HIP execution policies where ``WARP_SIZE`` would not be set
with ROCm 6.0.0 or newer.

**Deprecation Notices**

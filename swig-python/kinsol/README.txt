To use the KINSOL Python interfaces, you must have have python 3 and
numpy >= 1.15 installed. Additionally, you should have already built and
installed the KINSOL C library. **Note, you must build KINSOL with
double precision and 32-bit indices**.

The commands below will build and install a basic configuration of
KINSOL/SUNDIALS for use with the Python interfaces:

$ tar xvzf kinsol-5.0.0.tar.gz
$ cd kinsol-5.0.0
$ mkdir build && cd build
$ cmake -DCMAKE_INSTALL_PREFIX=/path/to/kinsol-installation -DCMAKE_BUILD_TYPE=Release -DSUNDIALS_INDEX_SIZE=32 ../
$ make
$ make install

Finally, to build the python interfaces for a desired package complete
the commands below (in the example below we build the KINSOL interface).

$ cd kinsol
$ export SUNDIALS_ROOT=/path/to/kinsol-installation
$ make kinsol
$ python setup.py install


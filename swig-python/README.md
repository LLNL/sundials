# SUNDIALS Python Interfaces

The SUNDIALS Python interfaces are generated with SWIG and support for Numpy.
So far only the KINSOL package, the generic SUNDIALS modules, and the serial
N_Vector are interfaced.


## Building the Interfaces

To use the SUNDIALS Python interfaces, you must have have python 3 and
numpy >= 1.15 installed. Additionally, you should have already built and
installed the SUNDIALS C library. **Note, you must build SUNDIALS with
double precision and 32-bit indices**.

Finally, to build the python interfaces for a desired package complete
the commands below (in the example below we build the KINSOL interface).

```bash
$ cd kinsol
$ export SUNDIALS_ROOT=/path/to/sundials-installation
$ make kinsol
```

## Using the Interface

Simply import the generated python module in your code (e.g. `import kinsol`).

### KINSOL

See the example `kinsol/kin_example.py`.

## SUNDIALS Developers

To generate the interfaces, it is required to have SWIG as well as the
`numpy.i` SWIG module that is packaged with the numpy source. To
generate the interfaces for a package, run:

```bash
$ cd <package>
$ export SUNDIALS_ROOT=/path/to/sundials-installation
$ make <package>_swig
```

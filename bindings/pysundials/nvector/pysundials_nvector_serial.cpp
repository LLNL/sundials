#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/vector.h>

#include <nvector/nvector_serial.h>
#include <sundials/sundials_core.h>

#include "pysundials_types.hpp"

namespace nb = nanobind;

void bind_nvector_serial(nb::module_& m)
{
#include "pysundials_nvector_serial_generated.hpp"
}

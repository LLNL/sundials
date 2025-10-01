from __future__ import annotations
import litgen
from litgen.internal.adapted_types import AdaptedFunction
from srcmlcpp.cpp_types import CppPublicProtectedPrivate


def ensure_return_policy_reference_for_pointers(func: AdaptedFunction) -> None:
    options = func.options

    pointer_types = options.sundials_pointer_types

    if not func.return_value_policy:
        if str(func.cpp_adapted_function.return_type) in pointer_types:
            func.return_value_policy = "reference"

# TODO(CJB): Since litgen is GPLv3, this script might have to be GPLv3.
# Will need to determine if this is the case or not.
# The outputs of the script, i.e. the generated code, are definitely
# not subject to GPLv3 though, and can use our standard license.

import copy
import re
from typing import Optional
from codemanip import code_utils
from srcmlcpp.cpp_types import CppParameter
from litgen.internal.adapt_function_params._lambda_adapter import LambdaAdapter
from litgen.internal.adapted_types import AdaptedFunction, AdaptedParameter


def adapt_array_pointer_to_std_vector(
    adapted_function: AdaptedFunction,
) -> Optional[LambdaAdapter]:
    """
    Adapts functions that take pointer_type* (pointer to pointer_type) parameters to instead accept
    std::vector<pointer_type> in the Python interface. Handles multiple such parameters in one function.
    """
    options = adapted_function.options
    function_name = adapted_function.cpp_adapted_function.function_name

    pointer_type_patterns = options.fn_params_array_pointers_to_std_vector

    # Map each parameter name to the pointer_type_pattern it matches
    param_to_pattern = {}
    for pointer_type_pattern in pointer_type_patterns:
        for param in adapted_function.adapted_parameters():
            cpp_type = param.cpp_element().full_type()
            param_name = param.cpp_element().decl.decl_name
            if re.match(pointer_type_pattern, str(param.cpp_element().decl)):
                param_to_pattern[param_name] = pointer_type_pattern

    if not param_to_pattern:
        return None

    lambda_adapter = LambdaAdapter()
    lambda_adapter.new_function_infos = copy.deepcopy(adapted_function.cpp_adapted_function)
    old_function_params: list[AdaptedParameter] = adapted_function.adapted_parameters()
    new_function_params: list[CppParameter] = []

    for old_param in old_function_params:
        old_cpp_type = old_param.cpp_element().full_type()
        param_name = old_param.cpp_element().decl.decl_name
        if param_name in param_to_pattern:
            base_type = remove_pointer(old_cpp_type)

            if base_type.endswith("1d"):
                dimensions = "*"
                base_type = base_type.replace("1d", "")
            elif base_type.endswith("2d"):
                dimensions = "**"
                base_type = base_type.replace("2d", "")

            new_param = copy.deepcopy(old_param.cpp_element())
            new_decl = new_param.decl
            new_decl.cpp_type.modifiers = []
            new_decl.cpp_type.specifiers = []
            new_decl.cpp_type.typenames = [f"std::vector<{base_type}>"]
            new_decl.initial_value_code = ""
            new_function_params.append(new_param)

            lambda_adapter.lambda_input_code += f"{base_type}{dimensions} {param_name}_ptr = reinterpret_cast<{base_type}{dimensions}>( {param_name}.empty() ? nullptr : {param_name}.data() );\n"
            lambda_adapter.adapted_cpp_parameter_list.append(f"{param_name}_ptr")
        else:
            new_function_params.append(old_param.cpp_element())
            lambda_adapter.adapted_cpp_parameter_list.append(param_name)

    lambda_adapter.new_function_infos.parameter_list.parameters = new_function_params
    lambda_adapter.lambda_name = (
        adapted_function.cpp_adapted_function.function_name + f"_adapt_arr_ptr_to_std_vector"
    )

    return lambda_adapter


def remove_pointer(cpp_type: str) -> str:
    """
    Removes pointer and const qualifiers from a C++ type string, e.g. 'const Foo *' -> 'Foo'.
    """
    # Remove 'const' and '*' and extra spaces
    t = cpp_type.replace("const", "").replace("*", "").strip()
    # Remove any duplicate spaces
    t = re.sub(r"\s+", " ", t)
    return t

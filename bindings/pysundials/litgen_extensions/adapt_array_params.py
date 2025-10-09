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
from .utils import is_array_param


def adapt_array_pointer_to_std_vector(
    adapted_function: AdaptedFunction,
) -> Optional[LambdaAdapter]:
    """
    Adapts functions that take pointer_type* (pointer to pointer_type) parameters to instead accept
    std::vector<pointer_type> in the Python interface. Handles multiple such parameters in one function.
    """
    options = adapted_function.options
    function_name = adapted_function.cpp_adapted_function.function_name

    has_array_param = False
    for param in adapted_function.adapted_parameters():
        if is_array_param(param.cpp_element().decl.decl_name):
            has_array_param = True

    if not has_array_param:
        return None

    lambda_adapter = LambdaAdapter()
    lambda_adapter.new_function_infos = copy.deepcopy(adapted_function.cpp_adapted_function)
    old_function_params: list[AdaptedParameter] = adapted_function.adapted_parameters()
    new_function_params: list[CppParameter] = []

    for old_param in old_function_params:
        param_name = old_param.cpp_element().decl.decl_name
        if is_array_param(param_name):
            base_type = " ".join(old_param.cpp_element().decl.cpp_type.typenames)
            ptr_dimension = old_param.cpp_element().decl.cpp_type.modifiers.count("*")

            if param_name.endswith("_1d"):
                if ptr_dimension != 1:
                    raise RuntimeWarning(
                        f"while processing {function_name} encountered pointer dimension ({ptr_dimension}) and param suffix (_1d) mismatch for param {old_param.cpp_element().full_type()} {param_name}"
                    )
                dimensions = "*"
            elif param_name.endswith("_2d"):
                if ptr_dimension != 2:
                    raise RuntimeWarning(
                        f"while processing {function_name} encountered pointer dimension ({ptr_dimension}) and param suffix (_2d) mismatch for param {old_param.cpp_element().full_type()} {param_name}"
                    )
                dimensions = "**"
            elif param_name.endswith("_3d"):
                if ptr_dimension != 3:
                    raise RuntimeWarning(
                        f"while processing {function_name} encountered pointer dimension ({ptr_dimension}) and param suffix (_3d) mismatch for param {old_param.cpp_element().full_type()} {param_name}"
                    )
                dimensions = "***"


            if is_float_type(base_type):
                new_param = copy.deepcopy(old_param.cpp_element())
                new_decl = new_param.decl
                new_decl.cpp_type.modifiers = []
                new_decl.cpp_type.specifiers = []
                new_decl.cpp_type.typenames = [f"nb::ndarray<{base_type}, nb::numpy, nb::ndim<1>, nb::c_contig>"]
                new_decl.initial_value_code = ""
                new_function_params.append(new_param)

                lambda_adapter.lambda_input_code += f"{base_type}{dimensions} {param_name}_ptr = reinterpret_cast<{base_type}{dimensions}>({param_name}.data());\n"
                lambda_adapter.adapted_cpp_parameter_list.append(f"{param_name}_ptr")
            else:
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


def is_float_type(cpp_type: str):
    return cpp_type in ["float", "double", "sunrealtype"]

def count_stars(cpp_type: str):
    """
    Counts the number of *
    """
    # Remove 'const' and '*' and extra spaces
    num_pointers = cpp_type.count("*")
    return num_pointers

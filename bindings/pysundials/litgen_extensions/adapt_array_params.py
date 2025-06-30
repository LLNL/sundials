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


def adapt_nvector_pointer_to_vector(adapted_function: AdaptedFunction) -> Optional[LambdaAdapter]:
    """
    Adapts functions that take N_Vector* (pointer to N_Vector) parameters to instead accept
    std::vector<N_Vector> in the Python interface.
    """
    options = adapted_function.options
    function_name = adapted_function.cpp_adapted_function.function_name

    # Find parameters that are N_Vector* (or const N_Vector*)
    nvector_ptr_params = []
    for param in adapted_function.adapted_parameters():
        cpp_type = param.cpp_element().full_type()
        # Match "N_Vector *" or "const N_Vector *"
        if re.match(r"(const\s+)?N_Vector\s*\*", cpp_type):
            nvector_ptr_params.append(param)

    if not nvector_ptr_params:
        return None
    
    lambda_adapter = LambdaAdapter()
    lambda_adapter.new_function_infos = copy.deepcopy(adapted_function.cpp_adapted_function)
    old_function_params: list[AdaptedParameter] = adapted_function.adapted_parameters()
    new_function_params: list[CppParameter] = []

    for old_param in old_function_params:
        cpp_type = old_param.cpp_element().full_type()
        if old_param in nvector_ptr_params:
            # Replace N_Vector* with std::vector<N_Vector>
            new_param = copy.deepcopy(old_param.cpp_element())
            new_decl = new_param.decl
            # Remove pointer and const qualifiers
            new_decl.cpp_type.modifiers = []
            new_decl.cpp_type.specifiers = []
            new_decl.cpp_type.typenames = ["std::vector<N_Vector>"]
            new_decl.initial_value_code = ""
            new_function_params.append(new_param)

            param_name = old_param.cpp_element().decl.decl_name
            param_vec_name = param_name + "_vec"
            lambda_adapter.lambda_input_code += (
                f"N_Vector* {param_name}_ptr = {param_name}.empty() ? nullptr : {param_name}.data();\n"
            )
            lambda_adapter.adapted_cpp_parameter_list.append(f"{param_name}_ptr")
        else:
            new_function_params.append(old_param.cpp_element())
            lambda_adapter.adapted_cpp_parameter_list.append(old_param.cpp_element().decl.decl_name)

    lambda_adapter.new_function_infos.parameter_list.parameters = new_function_params
    lambda_adapter.lambda_name = (
        adapted_function.cpp_adapted_function.function_name + "_adapt_nvector_ptr_to_vector"
    )

    return lambda_adapter



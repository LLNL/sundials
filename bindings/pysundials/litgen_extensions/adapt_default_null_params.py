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
from .utils import match_regex


def adapt_default_arg_pointer_with_default_null(
    adapted_function: AdaptedFunction,
) -> Optional[LambdaAdapter]:

    def is_param_optional(param, param__regex):
        if not code_utils.does_match_regex(param__regex, param.cpp_element().decl.decl_name):
            return False
        return True

    options = adapted_function.options

    function_name = adapted_function.cpp_adapted_function.function_name

    modifiable_fns = options.fn_params_optional_with_default_null.keys()
    modifiable_fns__regex = code_utils.join_string_by_pipe_char(modifiable_fns)

    if not code_utils.does_match_regex(modifiable_fns__regex, function_name):
        return None

    needs_adapt = False

    for old_adapted_param in adapted_function.adapted_parameters():
        if code_utils.does_match_regex(modifiable_fns__regex, function_name):
            needs_adapt = True

    if not needs_adapt:
        return None

    lambda_adapter = LambdaAdapter()

    lambda_adapter.new_function_infos = copy.deepcopy(adapted_function.cpp_adapted_function)

    old_function_params: list[AdaptedParameter] = adapted_function.adapted_parameters()

    new_function_params: list[CppParameter] = []

    for old_adapted_param in old_function_params:
        fn_key = match_regex(modifiable_fns__regex, function_name)[0].group()
        modifiable_fns_params__regex = code_utils.join_string_by_pipe_char(
            options.fn_params_optional_with_default_null[fn_key]
        )
        if is_param_optional(old_adapted_param, modifiable_fns_params__regex):
            param_original_type = old_adapted_param.cpp_element().full_type()
            param_name_value = (
                old_adapted_param.cpp_element().decl.decl_name + "_adapt_default_null"
            )
            param_name = old_adapted_param.cpp_element().decl.decl_name
            _i_ = options._indent_cpp_spaces()

            new_param = copy.deepcopy(old_adapted_param.cpp_element())
            new_decl = new_param.decl
            new_decl.cpp_type.modifiers = []
            new_decl.cpp_type.specifiers = []
            new_decl.cpp_type.typenames = [f"std::optional<{param_original_type}>"]
            new_decl.initial_value_code = "std::nullopt"

            new_function_params.append(new_param)

            lambda_input_code = f"""
                {param_original_type} {param_name_value} = nullptr;
                if ({param_name}.has_value())
                {_i_}{param_name_value} = {param_name}.value();
            """

            lambda_adapter.lambda_input_code += (
                code_utils.unindent_code(lambda_input_code, flag_strip_empty_lines=True) + "\n"
            )

            #
            # Fill adapted_cpp_parameter_list (those that will call the original C style function)
            #
            lambda_adapter.adapted_cpp_parameter_list.append(f"{param_name_value}")

        else:
            new_function_params.append(old_adapted_param.cpp_element())
            lambda_adapter.adapted_cpp_parameter_list.append(
                old_adapted_param.cpp_element().decl.decl_name
            )

    lambda_adapter.new_function_infos.parameter_list.parameters = new_function_params
    lambda_adapter.lambda_name = (
        adapted_function.cpp_adapted_function.function_name
        + "_adapt_optional_arg_with_default_null"
    )

    return lambda_adapter

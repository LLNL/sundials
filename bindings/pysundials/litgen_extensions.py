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


def generate_code(options, source_code):
    generated_code = litgen.generate_code(options, source_code)
    return generated_code


def load_opt_from_yaml(config_object, module, opt):
    opt_list = []
    all_section = config_object.get("all", [])
    if all_section:
        opt_list.extend(all_section.get(opt, []))
    module_section = config_object.get(module, [])
    if module_section:
        opt_list.extend(module_section.get(opt, []))
    return opt_list


def load_enum_exclusions_from_yaml(config_object, module):
    return load_opt_from_yaml(config_object, module, "enum_exclude_by_name__regex")


def load_class_exclusions_from_yaml(config_object, module):
    return load_opt_from_yaml(config_object, module, "class_exclude_by_name__regex")


def load_fn_exclusions_from_yaml(config_object, module):
    return load_opt_from_yaml(config_object, module, "fn_exclude_by_name__regex")


def load_macro_defines_from_yaml(config_object, module):
    return load_opt_from_yaml(config_object, module, "macro_define_include_by_name__regex")


def load_nullable_params_from_yaml(config_object, module):
    opt = "fn_params_optional_with_default_null"
    opt_dict = {}
    all_section = config_object.get("all", {})
    if all_section:
        opt_dict = opt_dict | all_section.get(opt, {})
    module_section = config_object.get(module, {})
    if module_section:
        opt_dict = opt_dict | module_section.get(opt, {})
    return opt_dict


def strip_sundials_export(code):
    return code.replace("SUNDIALS_EXPORT", "")


def change_long_int_to_long(code):
    return code.replace("long int", "long")


def change_sundials_types(code):
    code = re.sub(r"\bsunrealtype\b", "double", code)
    code = re.sub(r"\bsunindextype\b", "long", code)
    code = re.sub(r"\bsunbooleantype\b", "bool", code)
    return code


def preprocess_header(code):
    code = strip_sundials_export(code)
    code = change_long_int_to_long(code)
    code = change_sundials_types(code)
    return code


def match_regex(regex_str: str, word: str) -> bool:
    if regex_str.startswith("|"):
        regex_str = regex_str[1:]
    if len(regex_str) == 0:
        return False
    if word is None:
        return False
    matches = list(re.finditer(regex_str, word, re.MULTILINE))
    return matches


def is_param_optional(param, param__regex):
    if not code_utils.does_match_regex(param__regex, param.cpp_element().decl.decl_name):
        return False
    return True


def adapt_default_arg_pointer_with_default_null(adapted_function: AdaptedFunction) -> Optional[LambdaAdapter]:
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
            param_name_value = old_adapted_param.cpp_element().decl.decl_name + "_adapt_default_null"
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
            lambda_adapter.adapted_cpp_parameter_list.append(old_adapted_param.cpp_element().decl.decl_name)

    lambda_adapter.new_function_infos.parameter_list.parameters = new_function_params
    lambda_adapter.lambda_name = (
        adapted_function.cpp_adapted_function.function_name + "_adapt_optional_arg_with_default_null"
    )

    return lambda_adapter

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

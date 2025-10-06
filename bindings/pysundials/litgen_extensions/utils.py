# TODO(CJB): Since litgen is GPLv3, this script might have to be GPLv3.
# Will need to determine if this is the case or not.
# The outputs of the script, i.e. the generated code, are definitely
# not subject to GPLv3 though, and can use our standard license.

import litgen
from litgen.internal.adapted_types import AdaptedFunction, AdaptedParameter
import re


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


def load_pointer_types_from_yaml(config_object, module):
    opt = "sundials_pointer_types"
    opt_list = []
    all_section = config_object.get("all", [])
    if all_section:
        opt_list.extend(all_section.get(opt, []))
    module_section = config_object.get(module, [])
    if module_section:
        opt_list.extend(module_section.get(opt, []))
    return opt_list


def match_regex(regex_str: str, word: str) -> bool:
    if regex_str.startswith("|"):
        regex_str = regex_str[1:]
    if len(regex_str) == 0:
        return False
    if word is None:
        return False
    matches = list(re.finditer(regex_str, word, re.MULTILINE))
    return matches


def is_array_param(param_name):
    return param_name.endswith("_1d") or param_name.endswith("_2d") or param_name.endswith("_3d")

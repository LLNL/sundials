# TODO(CJB): Since litgen is GPLv3, this script might have to be GPLv3.
# Will need to determine if this is the case or not.
# The outputs of the script, i.e. the generated code, are definitely
# not subject to GPLv3 though, and can use our standard license.

import argparse
import srcmlcpp
import litgen
from codemanip import code_utils
import yaml


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


def load_class_exclusions_from_yaml(config_object, module):
    return load_opt_from_yaml(config_object, module, "class_exclude_by_name__regex")


def load_fn_exclusions_from_yaml(config_object, module):
    return load_opt_from_yaml(config_object, module, "fn_exclude_by_name__regex")


def load_macro_defines_from_yaml(config_object, module):
    return load_opt_from_yaml(config_object, module, "macro_define_include_by_name__regex")


def strip_sundials_export(code):
    return code.replace("SUNDIALS_EXPORT", "")


def preprocess_header(code):
    code = strip_sundials_export(code)
    return code


def main():
    parser = argparse.ArgumentParser(
        description="Generate Python bindings for SUNDIALS using litgen and nanobind."
    )
    parser.add_argument("config_yaml_path", type=str, help="Path to the generate.yaml config file")
    parser.add_argument("--dump-srcml", action="store_true", help="Dump the srcML XML for the parsed headers and exit")
    args = parser.parse_args()

    options = litgen.LitgenOptions()
    options.bind_library = litgen.BindLibraryType.nanobind
    options.python_run_black_formatter = True
    options.python_convert_to_snake_case = False
    options.comments_exclude = True
    options.srcmlcpp_options.code_preprocess_function = preprocess_header
    options.srcmlcpp_options.ignored_warning_parts.append(
        # "ops" functions pointers cause this warning, but we dont care cause we dont need to bind those.
        'A cpp element of type "function_decl" was stored as CppUnprocessed'
    )
    options.srcmlcpp_options.header_filter_preprocess_regions = True
    options.srcmlcpp_options.header_filter_acceptable__regex = "__cplusplus|_h_$|_h$|_H$|_H_$|hpp$|HPP$|hxx$|HXX$|SWIG$"

    config_yaml_path = args.config_yaml_path
    with open(config_yaml_path, "r") as yaml_file:
        config_object = yaml.safe_load(yaml_file).get("modules", [])
    if not config_object:
        raise RuntimeError(f"modules: section not found in {config_yaml_path}")

    # print(config_object)
    for module_name in config_object:
        if module_name == "all":
            continue

        module = config_object.get(module_name)

        options.class_exclude_by_name__regex = code_utils.join_string_by_pipe_char(
            load_class_exclusions_from_yaml(config_object, module_name)
        )

        options.fn_exclude_by_name__regex = code_utils.join_string_by_pipe_char(
            load_fn_exclusions_from_yaml(config_object, module_name)
        )

        options.macro_define_include_by_name__regex = code_utils.join_string_by_pipe_char(
            load_macro_defines_from_yaml(config_object, module_name)
        )

        # TODO(CJB): this does not seem to work
        # options.fn_return_force_policy_reference_for_pointers__regex = code_utils.join_string_by_pipe_char([
        #   "N_VNewEmpty",
        #   "N_VMake"
        # ])

        source_code = ""
        for file_path in module["headers"]:
            with open(file_path, "r") as file:
                source_code = source_code + file.read()
            source_code = source_code + "\n"

        if args.dump_srcml:
            srcmlcpp_options = options.srcmlcpp_options
            cpp_unit = srcmlcpp.code_to_cpp_unit(srcmlcpp_options, source_code)
            with open(f'{module["path"]}.xml', "w") as file:
                file.write(str(cpp_unit))
            continue

        generated_code = litgen.generate_code(options, source_code)

        if "path" in module:
            with open(module["path"], "w") as file:
                file.write(generated_code.glue_code)
                file.write(generated_code.pydef_code)
            # Macros in our headers (like SUNDIALS_EXPORT) throw off litgen
            # and pollute the stub file, so we don't generate it for now.
            # with open('pysundials.pyi', 'w') as file:
            #   file.write(generated_code.stub_code)
        else:
            print(generated_code.glue_code)
            print(generated_code.pydef_code)


if __name__ == "__main__":
    main()

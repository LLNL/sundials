# TODO(CJB): Since litgen is GPLv3, this script might have to be GPLv3.
# Will need to determine if this is the case or not.
# The outputs of the script, i.e. the generated code, are definitely
# not subject to GPLv3 though, and can use our standard license.

import argparse
import srcmlcpp
import litgen
from codemanip import code_utils
import yaml
from litgen_extensions import *


def main():
    parser = argparse.ArgumentParser(
        description="Generate Python bindings for SUNDIALS using litgen and nanobind."
    )
    parser.add_argument("config_yaml_path", type=str, help="Path to the generate.yaml config file")
    parser.add_argument(
        "--dump-srcml",
        action="store_true",
        help="Dump the srcML XML for the parsed headers and exit",
    )
    args = parser.parse_args()

    options = litgen.LitgenOptions()
    options.bind_library = litgen.BindLibraryType.nanobind
    options.python_run_black_formatter = True
    options.python_convert_to_snake_case = False

    # These are types in SUNDIALS which are typedefs to pointers to structs
    # The pointer types that are specific to a package are defined in the respective generate.yaml files.
    options.sundials_pointer_types = [
        "N_Vector",
        "SUNAdaptController"
        "SUNAdjointCheckpointScheme",
        "SUNAdjointStepper",
        "SUNContext",
        "SUNDomEigEstimator",
        "SUNErrHandler",
        "SUNLinearSolver",
        "SUNLogger",
        "SUNMatrix",
        "SUNMemoryHelper",
        "SUNNonlinearSolver",
        "SUNProfiler",
        "SUNStepper"
    ]

    # Don't capture comments from the source for generating Python doc strings
    options.comments_exclude = True

    # Export enum values to the package namespace
    options.enum_export_values = True

    # Enum values should be exported with their prefix. E.g., SUNDATAIOMODE_INMEM is exported for the enum SUNDataIOMode_, instead of just INMEM.
    options.enum_flag_remove_values_prefix = False

    # Allow const char to be nullable
    options.fn_params_const_char_pointer_with_default_null = True

    # Transform inplace modification of values, e.g. int CVodeGetNumSteps(void* cvode_mem, long int* num_steps), to CvodeGetNumSteps(cvode_mem) -> Tuple[int, long int]
    # Litgen original option is fn_params_output_modifiable_immutable_to_return__regex, but we use fn_params_output_modifiable_immutable_to_return__regex_custom
    # since we override the adapt_modifiable_immutable_to_return function adapter
    options.fn_params_output_modifiable_immutable_to_return__regex_custom = r".*"
    options.fn_return_force_policy_reference__callback = ensure_return_policy_reference_for_pointers

    # Force the functions that return pointers to use `nb::rv_policy::reference`
    options.fn_return_force_policy_reference_for_pointers__regex = r".*"

    # Don't create default constructors for any struct
    options.struct_create_default_named_ctor__regex = ""

    # Don't interface any struct or class member directly
    options.member_exclude_by_name__regex = r".*"

    # Our own custom function adapters
    options.fn_custom_adapters = [
        adapt_array_pointer_to_std_vector, # this must go first!
        adapt_modifiable_immutable_to_return,
        adapt_default_arg_pointer_with_default_null,
    ]

    options.srcmlcpp_options.code_preprocess_function = preprocess_header
    options.srcmlcpp_options.ignored_warning_parts.append(
        # "ops" functions pointers cause this warning, but we dont care cause we dont need to bind those.
        'A cpp element of type "function_decl" was stored as CppUnprocessed'
    )
    options.srcmlcpp_options.header_filter_preprocess_regions = True
    options.srcmlcpp_options.header_filter_acceptable__regex = (
        "__cplusplus|_h_$|_h$|_H$|_H_$|hpp$|HPP$|hxx$|HXX$|SWIG$"
    )

    config_yaml_path = args.config_yaml_path
    with open(config_yaml_path, "r") as yaml_file:
        config_object = yaml.safe_load(yaml_file).get("modules", [])
    if not config_object:
        raise RuntimeError(f"modules: section not found in {config_yaml_path}")

    for module_name in config_object:
        if module_name == "all":
            continue

        module = config_object.get(module_name)

        options.sundials_pointer_types.extend(load_pointer_types_from_yaml(
            config_object, module_name
        ))

        options.fn_params_optional_with_default_null = load_nullable_params_from_yaml(
            config_object, module_name
        )

        options.enum_exclude_by_name__regex = code_utils.join_string_by_pipe_char(
            load_enum_exclusions_from_yaml(config_object, module_name)
        )

        options.class_exclude_by_name__regex = code_utils.join_string_by_pipe_char(
            load_class_exclusions_from_yaml(config_object, module_name)
        )

        options.fn_exclude_by_name__regex = code_utils.join_string_by_pipe_char(
            load_fn_exclusions_from_yaml(config_object, module_name)
        )

        options.macro_define_include_by_name__regex = code_utils.join_string_by_pipe_char(
            load_macro_defines_from_yaml(config_object, module_name)
        )

        source_code = ""
        for file_path in module["headers"]:
            with open(file_path, "r") as file:
                source_code = source_code + file.read()
            source_code = source_code + "\n"

        if args.dump_srcml:
            srcmlcpp_options = options.srcmlcpp_options
            cpp_unit = srcmlcpp.code_to_cpp_unit(srcmlcpp_options, source_code)
            with open(f'{module["path"]}.xml', "w") as file:
                file.write(cpp_unit.str_code())
            continue

        generated_code = litgen.generate_code(options, source_code)

        if "path" in module:
            with open(module["path"], "w") as file:
                file.write(generated_code.glue_code)
                file.write(generated_code.pydef_code)
            # TODO(CJB): Not sure how we would combine generated and custom code for stubs
            # with open(f'{module["path"]}.pyi', 'w') as file:
            #   file.write(generated_code.stub_code)
        else:
            print(generated_code.glue_code)
            print(generated_code.pydef_code)


if __name__ == "__main__":
    main()

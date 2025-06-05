import argparse
import litgen
from codemanip import code_utils
import yaml


def generate_code(options, source_code):
    generated_code = litgen.generate_code(options, source_code)
    return generated_code


def load_exclusions_from_yaml(config_object, module):
    exclusions = []
    all_section = config_object.get("all", [])
    if all_section:
        exclusions.extend(all_section.get("fn_exclude_by_name__regex", []))
    module_section = config_object.get(module, [])
    if module_section:
        exclusions.extend(module_section.get("fn_exclude_by_name__regex", []))
    return exclusions


def main():
    parser = argparse.ArgumentParser(
        description="Generate Python bindings for SUNDIALS using litgen and nanobind."
    )
    parser.add_argument("config_yaml_path", type=str, help="Path to the generate.yaml config file")
    args = parser.parse_args()

    options = litgen.LitgenOptions()
    options.bind_library = litgen.BindLibraryType.nanobind
    options.python_run_black_formatter = True
    options.python_convert_to_snake_case = False

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

        options.fn_exclude_by_name__regex = code_utils.join_string_by_pipe_char(
            load_exclusions_from_yaml(config_object, module_name)
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

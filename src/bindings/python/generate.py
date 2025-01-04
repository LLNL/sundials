import argparse
import litgen
from codemanip import code_utils

def generate_code(options, source_code):
  generated_code = litgen.generate_code(options, source_code)
  return generated_code

def main():
  parser = argparse.ArgumentParser(description='Generate Python bindings for SUNDIALS using litgen and nanobind.')
  parser.add_argument('source_code', type=str, help='Path to the source code file')
  parser.add_argument('-o', '--output_path', type=str, help='Path to the output file', default=None)
  args = parser.parse_args()

  with open(args.source_code, 'r') as file:
    source_code = file.read()

  options = litgen.LitgenOptions()
  options.bind_library = litgen.BindLibraryType.nanobind
  options.python_run_black_formatter = True
  options.python_convert_to_snake_case = False

  options.fn_exclude_by_name__regex = code_utils.join_string_by_pipe_char([
    "Free",
    "N_VDestroy",
    "ArrayPointer",
    "BufPack",
    "BufUnpack",
    "VectorArray",
    "Multi",
    "LinearCombination"
  ])

  # TODO(CJB): this does not seem to work
  # options.fn_return_force_policy_reference_for_pointers__regex = code_utils.join_string_by_pipe_char([
  #   "N_VNewEmpty",
  #   "N_VMake"
  # ])

  # options.enum_make_arithmetic__regex = "^N_Vector_ID_Enum$"

  generated_code = litgen.generate_code(options, source_code)

  if args.output_path:
    with open(args.output_path, 'w') as file:
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

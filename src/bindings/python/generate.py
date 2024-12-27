import litgen
from codemanip import code_utils

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
  "Multi"
])

# TODO(CJB): this does not seem to work
# options.fn_return_force_policy_reference_for_pointers__regex = code_utils.join_string_by_pipe_char([
#   "N_VNewEmpty",
#   "N_VMake"
# ])

source_code = ''
with open('/Volumes/Workspace/SUNDIALS/repos/sundials/include/sundials/sundials_nvector.h', 'r') as file:
  source_code = file.read()

generated_code = litgen.generate_code(options, source_code)

# print(generated_code.stub_code)
print(generated_code.pydef_code)
# print(generated_code.glue_code)

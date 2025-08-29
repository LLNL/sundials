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


def strip_sundials_export(code):
    return code.replace("SUNDIALS_EXPORT", "")


def change_long_int_to_long(code):
    return code.replace("long int", "long")


def change_sundials_types(code):
    code = re.sub(r"\bsunrealtype\b", "double", code)
    code = re.sub(r"\bsunindextype\b", "long", code)
    code = re.sub(r"\bsunbooleantype\b", "int", code)
    code = re.sub(r"\bFILE\*\b", "std::shared_ptr<FILE>", code)
    return code


def preprocess_header(code):
    code = strip_sundials_export(code)
    code = change_long_int_to_long(code)
    code = change_sundials_types(code)
    return code

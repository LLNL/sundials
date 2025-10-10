# TODO(CJB): Since litgen is GPLv3, this script might have to be GPLv3.
# Will need to determine if this is the case or not.
# The outputs of the script, i.e. the generated code, are definitely
# not subject to GPLv3 though, and can use our standard license.

import re

def strip_sundials_export(code):
    return code.replace("SUNDIALS_EXPORT", "")


def change_long_int_to_long(code):
    return code.replace("long int", "long")


def preprocess_header(code):
    code = strip_sundials_export(code)
    code = change_long_int_to_long(code)
    return code

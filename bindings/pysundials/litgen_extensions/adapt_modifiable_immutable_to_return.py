# Copied from https://github.com/pthom/litgen/blob/2fafb1e202689df5bb7b0457c5d2b81527fb3829/src/litgen/internal/adapt_function_params/adapt_modifiable_immutable_to_return.py
# and modified to exclude _<1|2|3>d parameters in sundials.


from __future__ import annotations
import copy
from typing import Optional

from munch import Munch  # type: ignore

from codemanip import code_utils

from srcmlcpp.cpp_types import CppParameter, CppType

from litgen.internal.adapt_function_params._lambda_adapter import LambdaAdapter
from litgen.internal.adapted_types import AdaptedFunction, AdaptedParameter

from .utils import is_array_param


def adapt_modifiable_immutable_to_return(
    adapted_function: AdaptedFunction,
) -> Optional[LambdaAdapter]:
    """
    We want to adapt functions params that use modifiable pointer or reference to a type that is immutable in python.
    For example
        `int foo(int* value)`
        `void foo(float& value)`
        `void foo(string& value)`

    But not
        `int foo(int* value_1d)`

    On the C++ side, these params are modifiable by the function.

    Note: immutable data types in python are
        - Int, Float, String (currently handled here)
        - Complex, Bytes (not yet handed here)
        - Tuple (not handled here)

    For the function
        bool SliderInt(const char* label, int* value);
    We will generate an adapter lambda that looks like
        auto SliderInt_adapt_modifiable_immutable_to_return =
        [](const char* label, int value) -> std::tuple<bool, int>
        {
            value_adapt_modifiable = &value;
            bool r = SliderInt(label, &value_adapt_pointer);
            return  std::make_tuple(r, value);
        };
        return SliderInt_adapt_modifiable_immutable_to_return(input);
    """
    options = adapted_function.options

    function_name = adapted_function.cpp_adapted_function.function_name
    if not code_utils.does_match_regex_or_matcher(
        options.fn_params_output_modifiable_immutable_to_return__regex_custom, function_name
    ):
        return None

    needs_adapt = False
    any_param_needs_rv_reference = False

    # is_modifiable_python_immutable_ref_or_pointer only recognizes some built in types
    # but we would like to also transform functions that return SUNDIALS objects as outputs,
    # e.g., SUNErrCode SUNSparseMatrix_ToCSC(SUNMatrix A, SUNMatrix* Bout) becomes
    # std::tuple<SUNErrCode, SUNMatrix> SUNSparseMatrix_ToCSC(SUNMatrix A, SUNMatrix B); . 
    def is_immutable_ref_or_pointer(param):
        a = param.is_modifiable_python_immutable_ref_or_pointer() 
        b = param.cpp_element().decl.cpp_type.name_without_modifier_specifier() in options.sundials_pointer_types
        c = param.cpp_element().decl.cpp_type.modifiers == ["*"]
        if a or (b and c):
            return True, b and c
        return False, False

    for old_adapted_param in adapted_function.adapted_parameters():
        if (
            is_immutable_ref_or_pointer(old_adapted_param)[0]
            or old_adapted_param.is_modifiable_python_immutable_fixed_size_array()
        ) and not is_array_param(old_adapted_param.cpp_element().decl.decl_name):
            needs_adapt = True

    if not needs_adapt:
        return None

    lambda_adapter = LambdaAdapter()

    lambda_adapter.new_function_infos = copy.deepcopy(adapted_function.cpp_adapted_function)

    old_function_params: list[AdaptedParameter] = adapted_function.adapted_parameters()

    new_function_params: list[CppParameter] = []
    new_output_function_params: list[CppParameter] = []
    for old_adapted_param in old_function_params:
        is_imm_ref_or_ptr, needs_rv_reference = is_immutable_ref_or_pointer(old_adapted_param)
        any_param_needs_rv_reference = any_param_needs_rv_reference or needs_rv_reference
        if is_imm_ref_or_ptr:
            is_pointer = old_adapted_param.cpp_element().decl.cpp_type.modifiers == ["*"]

            # For signatures like
            #       void foo(bool * flag = NULL);
            # the python param type will be type Optional[BoxedBool]
            def compute_is_optional_boxed_type(old_adapted_param=old_adapted_param, is_pointer=is_pointer) -> bool:  # type: ignore
                initial_value_cpp = old_adapted_param.cpp_element().decl.initial_value_code
                is_initial_value_null = initial_value_cpp in ["NULL", "nullptr"]
                return is_pointer and is_initial_value_null

            is_optional_type = compute_is_optional_boxed_type()

            #
            # Create new calling param same type, without pointer or reference
            #
            param_name_value = old_adapted_param.cpp_element().decl.decl_name + "_adapt_modifiable"

            new_param = copy.deepcopy(old_adapted_param.cpp_element())
            old_decl = old_adapted_param.cpp_element().decl
            new_decl = new_param.decl
            new_decl.decl_name = param_name_value
            new_decl.cpp_type.typenames = old_decl.cpp_type.typenames
            new_decl.cpp_type.modifiers = []
            new_decl.cpp_type.specifiers = []
            if is_optional_type:
                new_decl.cpp_type.typenames = [f"std::optional<{new_decl.cpp_type.str_code()}>"]
                new_decl.initial_value_code = "std::nullopt"
                new_function_params.append(new_param)

            #
            # Fill lambda_input_code
            #
            """
            lambda_input_code will look like
                value_adapt_modifiable;
            """
            # param_original_type = old_adapted_param.cpp_element().full_type()
            param_original_type = " ".join(old_adapted_param.cpp_element().decl.cpp_type.typenames)
            param_name = old_adapted_param.cpp_element().decl.decl_name
            _i_ = options._indent_cpp_spaces()
            if is_optional_type:
                lambda_input_code = f"""
                    {param_original_type} {param_name_value} = nullptr;
                    if ({param_name}.has_value())
                    {_i_}{param_name_value} = & (*{param_name});
                """
            else:
                if is_pointer:
                    lambda_input_code = f"""
                        {param_original_type} {param_name_value};
                        """
                else:
                    lambda_input_code = f"""
                        {param_original_type} {param_name_value};
                        """

            lambda_adapter.lambda_input_code += (
                code_utils.unindent_code(lambda_input_code, flag_strip_empty_lines=True) + "\n"
            )

            #
            # Fill adapted_cpp_parameter_list (those that will call the original C style function)
            #
            lambda_adapter.adapted_cpp_parameter_list.append(f"&{param_name_value}")

            #
            # Fill new_output_function_params
            #
            new_output_function_params.append(new_param)

        elif old_adapted_param.is_modifiable_python_immutable_fixed_size_array():
            #
            # Create new calling param same type, without pointer or reference
            #
            new_param = copy.deepcopy(old_adapted_param.cpp_element())
            array_type = new_param.decl.cpp_type.name_without_modifier_specifier()
            array_size = new_param.decl.c_array_size_as_int()
            decl_name = new_param.decl.decl_name
            from srcmlcpp import srcmlcpp_main

            new_decl_str = f"std::array<{array_type}, {array_size}> {decl_name}"
            new_decl = srcmlcpp_main.code_first_decl(options.srcmlcpp_options, new_decl_str)
            new_param.decl = new_decl

            new_function_params.append(new_param)

            #
            # Fill lambda_input_code
            #
            param_original_type = old_adapted_param.cpp_element().full_type()
            param_name_value = old_adapted_param.cpp_element().decl.decl_name + "_adapt_modifiable"
            param_name = old_adapted_param.cpp_element().decl.decl_name
            _i_ = options._indent_cpp_spaces()
            lambda_input_code = f"""
                {param_original_type} * {param_name_value} = {param_name}.data();
                """

            lambda_adapter.lambda_input_code += (
                code_utils.unindent_code(lambda_input_code, flag_strip_empty_lines=True) + "\n"
            )

            #
            # Fill adapted_cpp_parameter_list (those that will call the original C style function)
            #
            lambda_adapter.adapted_cpp_parameter_list.append(f"{param_name_value}")

            #
            # Fill new_output_function_params
            #
            new_output_function_params.append(new_param)

        else:
            new_function_params.append(old_adapted_param.cpp_element())
            lambda_adapter.adapted_cpp_parameter_list.append(
                old_adapted_param.cpp_element().decl.decl_name
            )

    lambda_adapter.new_function_infos.parameter_list.parameters = new_function_params
    lambda_adapter.lambda_name = (
        adapted_function.cpp_adapted_function.function_name
        + "_adapt_modifiable_immutable_to_return"
    )

    #
    # Adapt lambda return type
    #

    if any_param_needs_rv_reference:
        adapted_function.return_value_policy = "reference"

    old_function = adapted_function.cpp_adapted_function
    was_void_return_type = old_function.returns_void()

    def new_return_type_and_names() -> tuple[str, str]:
        number_of_return_params = len(new_output_function_params)
        if not was_void_return_type:
            number_of_return_params += 1

        assert number_of_return_params > 0

        if number_of_return_params == 0:
            new_return_type = new_output_function_params[0].decl.cpp_type.str_code()
            new_return_name = new_output_function_params[0].decl.decl_name
            return new_return_type, new_return_name
        else:
            type_list = []
            name_list = []
            if not was_void_return_type:
                type_list.append(old_function.return_type.str_code())
                name_list.append("r")
            for new_output_function_param in new_output_function_params:
                type_list.append(new_output_function_param.full_type())
                name_list.append(new_output_function_param.decl.decl_name)

            if number_of_return_params > 1:
                type_list_str = ", ".join(type_list)
                new_return_type = f"std::tuple<{type_list_str}>"

                name_list_str = ", ".join(name_list)
                new_return_name = f"std::make_tuple({name_list_str})"
            else:
                new_return_type = type_list[0]
                new_return_name = name_list[0]

            return new_return_type, new_return_name

    lambda_adapter.new_function_infos.return_type = CppType(
        lambda_adapter.new_function_infos.return_type
    )
    new_return_types, new_return_names = new_return_type_and_names()
    lambda_adapter.new_function_infos.return_type.typenames = [new_return_types]

    #
    # Fill lambda_template_end
    #
    def fill_lambda_template_end() -> str:
        """
        lambda_template_end will look like:
                bool r = SliderInt(label, &value_adapt_pointer);
                return std::make_tuple(r, value);
        """
        template_code = code_utils.unindent_code(
            """
        {maybe_store_function_output}{lambda_to_call}({params});
        return {output_tuple};
        """,
            flag_strip_empty_lines=True,
        )

        replacements = Munch()

        # fill maybe_store_function_output
        if was_void_return_type:
            replacements.maybe_store_function_output = ""
        else:
            old_return_type = old_function.str_full_return_type()
            replacements.maybe_store_function_output = f"{old_return_type} r = "

        # fill function_name
        if adapted_function.lambda_to_call is not None:
            replacements.lambda_to_call = adapted_function.lambda_to_call
        else:
            if old_function.is_method():
                replacements.lambda_to_call = (
                    "self." + old_function.function_name_with_specialization()
                )
            else:
                replacements.lambda_to_call = (
                    old_function.qualified_function_name_with_specialization()
                )

        # fill params
        replacements.params = ", ".join(lambda_adapter.adapted_cpp_parameter_list)

        # fill output_tuple
        replacements.output_tuple = new_return_names

        code = code_utils.replace_in_string(template_code, replacements)
        return code

    lambda_adapter.lambda_template_end = fill_lambda_template_end()

    return lambda_adapter

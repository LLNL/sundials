#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Programmer(s): Cody Balos and David J. Gardner @ LLNL
# -----------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2024, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------------------
# Module of Python functions that may be useful to SUNDIALS developers writing
# scripts to parse logs produced by SUNLogger.
# -----------------------------------------------------------------------------

import re
import numpy as np

def convert_to_num(s):
    """Try to convert a string to an int or float"""

    try:
        return np.longlong(s)
    except ValueError:
        try:
            return np.double(s)
        except ValueError:
            return s

def parse_logfile_payload(payload, line_number, all_lines, array_indicator="(:)"):
    """
    This function parses the payload of in a SUNDIALS log file line
    into a dictionary. The payload of a SUNDIALS log file line
    is the part after all the [ ] brackets.
    """
    kvpstrs = payload.split(",")
    kvp_dict = {}
    for kvpstr in kvpstrs:
        kvp = kvpstr.split("=")
        if len(kvp) == 1:
            # Check for empty payload
            if not kvp[0].strip():
                continue
            kvp_dict[kvp[0].strip()] = ""
        else:
            key, value = kvp
            values = []
            if array_indicator in key:
                for line in all_lines[line_number + 1 :]:
                    if line.startswith("["):
                        break
                    values.append(np.double(line))
                kvp_dict[key.strip()] = values
            else:
                kvp_dict[key.strip()] = value.strip()
    return kvp_dict


def parse_logfile_line(line, line_number, all_lines):
    """
    This function takes a line from a SUNDIALS log file and parses it into a dictionary.
    A log file line has the form:
      [loglvl][rank][scope][label] key1 = value, key2 = value
    Log file payloads can be multiline if they are an array/vector with one value per line.
    I.e.
      [loglvl][rank][scope][label] y(:)
      y_1
      y_2
      ...
    """
    pattern = re.compile(r'\[(\w+)\]\[(rank \d+)\]\[(.*)\]\[(.*)\](.*)')
    matches = pattern.findall(line)
    line_dict = {}
    if matches:
        line_dict["loglvl"] = matches[0][0]
        line_dict["rank"] = convert_to_num(matches[0][1].split()[1])
        line_dict["scope"] = matches[0][2]
        line_dict["label"] = matches[0][3]
        line_dict["payload"] = parse_logfile_payload(
            matches[0][4], line_number, all_lines
        )
    return line_dict


def log_file_to_list(filename):
    """
    This function takes a SUNDIALS log file and creates a list where each list
    element represents an integrator step attempt.

    E.g.,
      [
          {
              "loglvl": "DEBUG",
              "rank": "rank 0",
              "scope": "<step_scope_txt>",
              "label": "enter-step-attempt-loop",
              "payload": {"step": "0", "h": "1e-06", "q": "1", "t_n": "0"},
          }, ...
      ]
    """
    with open(filename, "r") as logfile:

        # List of step attempts, each entry is a dictionary for one attempt
        step_attempts = []

        # Stack of lists to handle sublists for different log scopes e.g., stage
        # data, algebraic solver data, fast integrator data, etc.
        #
        # stack[-1] is the active list of dictionaries
        # stage[-1][-1] is the active dictionary
        stack = []

        # Make step attempt list the active list
        stack.append(step_attempts)

        # Time level for nested integrators e.g., MRI methods
        level = 0

        # Read the log file
        all_lines = logfile.readlines()

        for line_number, line in enumerate(all_lines):

            line_dict = parse_logfile_line(line.rstrip(), line_number, all_lines)

            if not line_dict:
                continue

            label = line_dict["label"]

            if label == "begin-step-attempt":
                # Set the current time level
                line_dict["payload"]["level"] = level
                # Add new step attempt dictionary to the active list
                stack[-1].append(line_dict["payload"])
                continue
            elif label == "end-step-attempt":
                # Update active step attempt dictionary
                stack[-1][-1].update(line_dict["payload"])
                continue

            if label == "begin-stage":
                # Add stage sublist to the active dictionary
                if "stages" not in stack[-1][-1]:
                    stack[-1][-1]["stages"] = []
                # Make the stage sublist the active list
                stack.append(stack[-1][-1]["stages"])
                # Add new stage dictionary to the active list
                stack[-1].append(line_dict["payload"])
                continue
            elif label == "end-stage":
                # Update the active stage dictionary
                stack[-1][-1].update(line_dict["payload"])
                # Deactivate stage list
                stack.pop()
                continue

            if label == "begin-fast-steps":
                level += 1
                key = f"time-level-{level}"
                # Add fast step sublist to the active dictionary
                if key not in stack[-1][-1]:
                    stack[-1][-1][key] = []
                # Make the fast step sublist the active list
                stack.append(stack[-1][-1][key])
                continue
            elif label == "end-fast-steps":
                level -= 1
                # Deactivate fast step list
                stack.pop()
                continue

            if label == "begin-nonlinear-solve":
                if "nonlinear-solve" not in stack[-1][-1]:
                    stack[-1][-1]["nonlinear-solve"] = {}
                stack[-1][-1]["nonlinear-solve"].update(line_dict["payload"])
                continue
            elif label == "end-nonlinear-solve":
                stack[-1][-1]["nonlinear-solve"].update(line_dict["payload"])
                continue

            if (label == "begin-nonlinear-iterate"):
                if "iterations" not in stack[-1][-1]["nonlinear-solve"]:
                    stack[-1][-1]["nonlinear-solve"]["iterations"] = []
                # Make the stage sublist the active list
                stack.append(stack[-1][-1]["nonlinear-solve"]["iterations"])
                # Add new solver iteration dictionary
                stack[-1].append(line_dict["payload"])
                continue
            elif (label == "end-nonlinear-iterate"):
                # Update the active iteration dictionary
                stack[-1][-1].update(line_dict["payload"])
                stack.pop()
                continue

            if label == "begin-linear-solve":
                if "linear-solve" not in stack[-1][-1]:
                    stack[-1][-1]["linear-solve"] = {}
                stack[-1][-1]["linear-solve"].update(line_dict["payload"])
                continue
            elif label == "end-linear-solve":
                stack[-1][-1]["linear-solve"].update(line_dict["payload"])
                continue

            if (label == "begin-linear-iterate"):
                if "iterations" not in stack[-1][-1]["linear-solve"]:
                    stack[-1][-1]["linear-solve"]["iterations"] = []
                # Make the stage sublist the active list
                stack.append(stack[-1][-1]["linear-solve"]["iterations"])
                # Add new solver iteration dictionary
                stack[-1].append(line_dict["payload"])
                continue
            elif (label == "end-linear-iterate"):
                # Update the active iteration dictionary
                stack[-1][-1].update(line_dict["payload"])
                stack.pop()
                continue

            # Update current step attempt entry with intermediate output
            stack[-1][-1].update(line_dict["payload"])

    return step_attempts


def print_log(log, indent=0):
    """
    This function prints the list of entries from a log file.
    """

    spaces = indent * " "

    for entry in log:
        print(f"{spaces}{{")
        for key in entry:
            if type(entry[key]) is list:
                print(f"{spaces}{key} :")
                print(f"{spaces}[")
                subindent = indent + 2
                print_log(entry[key], indent=subindent)
                print(f"{spaces}]")
            elif type(entry[key]) is dict:
                print(f"{spaces}{key} :")
                print(f"{spaces}{{")
                subindent = indent + 2
                subspaces = subindent * " "
                for subkey in entry[key]:
                    if type(entry[key][subkey]) is list:
                        print(f"{subspaces}{subkey} :")
                        print(f"{subspaces}[")
                        subsubindent = subindent + 2
                        print_log(entry[key][subkey], indent=subsubindent)
                        print(f"{subspaces}]")
                    else:
                        print(f"{subspaces}{subkey} : {entry[key][subkey]}")
                print(f"{spaces}}}")
            else:
                print(f"{spaces}{key} : {entry[key]}")
        print(f"{spaces}}}")


def get_history(log, key, step_status = None, time_range = None,
                step_range = None):
    """
    This function extracts the step/time series of the requested value.
    """

    steps = []
    times = []
    values = []

    for entry in log:

        step = np.longlong(entry['step'])
        time = np.double(entry['t_n'])

        if time_range is not None:
            if time < time_range[0] or time > time_range[1]:
                continue

        if step_range is not None:
            if step < step_range[0] or step > step_range[1]:
                continue

        if step_status is not None:
            if step_status not in entry['status']:
                continue

        if key not in entry:
            continue

        steps.append(step)
        times.append(time)
        values.append(convert_to_num(entry[key]))

        if "stages" in entry:
            for s in entry["stages"]:
                next_level_key = f'time-level-{entry["level"] + 1}'
                if next_level_key in s:
                    sub_steps, sub_times, sub_values = get_history(s[next_level_key], key)
                    steps.extend(sub_steps)
                    times.extend(sub_times)
                    values.extend(sub_values)

    return steps, times, values

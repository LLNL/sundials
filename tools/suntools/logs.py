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
from collections import ChainMap


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
    Parse the payload of a SUNDIALS log file line into a dictionary. The payload
    of a SUNDIALS log file line is the part after all the [ ] brackets.
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
                    if line.startswith("[") or not line.strip():
                        break
                    values.append(np.double(line))
                kvp_dict[key.strip()] = values
            else:
                kvp_dict[key.strip()] = value.strip()
    return kvp_dict


def parse_logfile_line(line, line_number, all_lines):
    """
    Parse a line from a SUNDIALS log file into a dictionary.

    A log file line has the form:
      [loglvl][rank][scope][label] key1 = value, key2 = value
    The log line payload (everything after the brackets) can be multiline with
    one value per line for keys corresponding to an array/vector:
      [loglvl][rank][scope][label] y(:) =
      y_1
      y_2
      ...
    """
    pattern = re.compile(r"\[(\w+)\]\[(rank \d+)\]\[(.*)\]\[(.*)\](.*)")
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


class StepData:
    """
    Helper class for parsing a step attempt from a SUNDIALS log file into a
    hierarchical dictionary where entries may be lists of dictionaries.
    """

    def __init__(self):
        self.container = [ChainMap()]
        self.parent_keys = ["main"]

    def __repr__(self):
        tmp = "Container:"
        for entry in self.container:
            tmp += f"\n  {entry}"
        tmp += "\nParent Keys:"
        for entry in self.parent_keys:
            tmp += f"\n  {entry}"
        return tmp

    def update(self, data):
        """Update the active dictionary"""
        self.container[-1].update(data)

    def open_dict(self, key):
        """Activate a dictionary"""
        self.container[-1][key] = {}
        self.container[-1] = self.container[-1].new_child(self.container[-1][key])

    def close_dict(self):
        """Deactivate the active dictionary"""
        self.container[-1] = self.container[-1].parents

    def open_list(self, key):
        """Activate a list of dictionaries"""
        if key in self.container[-1]:
            self.container.append(ChainMap())
        else:
            self.container[-1][key] = []
            self.container.append(ChainMap())
        self.parent_keys.append(key)

    def close_list(self):
        """Deactivate a the active list"""
        tmp = self.container[-1].maps[0]
        self.container[-2][self.parent_keys[-1]].append(tmp)
        self.parent_keys.pop()
        self.container.pop()

    def get_step(self):
        """Get the step dictionary and reset the container"""
        tmp = self.container.pop().maps[0]
        self.container = [ChainMap()]
        self.parent_keys = ["main"]
        return tmp


def log_file_to_list(filename):
    """
    This function takes a SUNDIALS log file and creates a list where each list
    element represents an integrator step attempt.

    E.g.,
      [
        {
          step   : 1,
          tn     : 0.0,
          stages : [ {stage : 1, tcur : 0.0, ...}, {stage : 2, tcur : 0.5, ...}, ...]
          ...
        }, ...
      ]
    """
    with open(filename, "r") as logfile:

        # List of step attempts, each entry is a dictionary for one attempt
        step_attempts = []

        # Time level for nested integrators e.g., MRI methods
        level = 0

        # Partition for split integrators e.g., operator splitting methods
        partition = 0

        # Read the log file
        all_lines = logfile.readlines()

        # Create instance of helper class for building attempt dictionary
        s = StepData()

        for line_number, line in enumerate(all_lines):

            line_dict = parse_logfile_line(line.rstrip(), line_number, all_lines)

            if not line_dict:
                continue

            label = line_dict["label"]

            if label == "begin-step-attempt":
                line_dict["payload"]["level"] = level
                if level > 0:
                    s.open_list(f"time-level-{level}")
                if partition > 0:
                    s.open_list(f"evolve")
                s.update(line_dict["payload"])
                continue
            elif label == "end-step-attempt":
                s.update(line_dict["payload"])
                if level > 0 or partition > 0:
                    s.close_list()
                else:
                    step_attempts.append(s.get_step())
                continue

            if label == "begin-sequential-method":
                s.open_list("sequential methods")
                s.update(line_dict["payload"])
                continue
            elif label == "end-sequential-method":
                s.update(line_dict["payload"])
                s.close_list()
                continue

            if label == "begin-partition":
                s.open_list("partitions")
                s.update(line_dict["payload"])
                partition += 1
                continue
            elif label == "end-partition":
                s.update(line_dict["payload"])
                s.close_list()
                partition -= 1
                continue

            if label == "begin-nonlinear-solve":
                s.open_dict("nonlinear-solve")
                s.update(line_dict["payload"])
                continue
            elif label == "end-nonlinear-solve":
                s.update(line_dict["payload"])
                s.close_dict()
                continue

            if label == "begin-nonlinear-iterate":
                s.open_list("iterations")
                s.update(line_dict["payload"])
                continue
            elif label == "end-nonlinear-iterate":
                s.update(line_dict["payload"])
                s.close_list()
                continue

            if label == "begin-linear-solve":
                s.open_dict("linear-solve")
                s.update(line_dict["payload"])
                continue
            elif label == "end-linear-solve":
                s.update(line_dict["payload"])
                s.close_dict()
                continue

            if label == "begin-linear-iterate":
                s.open_list("iterations")
                s.update(line_dict["payload"])
                continue
            elif label == "end-linear-iterate":
                s.update(line_dict["payload"])
                s.close_list()
                continue

            if label == "begin-group":
                s.open_list("groups")
                s.update(line_dict["payload"])
                continue
            elif label == "end-group":
                s.update(line_dict["payload"])
                s.close_list()
                continue

            if label == "begin-stage":
                s.open_list("stages")
                s.update(line_dict["payload"])
                continue
            elif label == "end-stage":
                s.update(line_dict["payload"])
                s.close_list()
                continue

            if label == "begin-fast-steps":
                level += 1
                continue
            elif label == "end-fast-steps":
                level -= 1
                continue

            if label == "begin-mass-linear-solve":
                s.open_dict("mass-linear-solve")
                s.update(line_dict["payload"])
                continue
            elif label == "end-mass-linear-solve":
                s.update(line_dict["payload"])
                s.close_dict()
                continue

            if label == "begin-compute-solution":
                s.open_dict("compute-solution")
                s.update(line_dict["payload"])
                continue
            elif label == "end-compute-solution":
                s.update(line_dict["payload"])
                s.close_dict()
                continue

            if label == "begin-compute-embedding":
                s.open_dict("compute-embedding")
                s.update(line_dict["payload"])
                continue
            elif label == "end-compute-embedding":
                s.update(line_dict["payload"])
                s.close_dict()
                continue

            s.update(line_dict["payload"])

    return step_attempts


def print_log_list(a_list, indent=0):
    """Print list value in step attempt dictionary"""
    spaces = (indent + 2) * " "
    for entry in a_list:
        if type(entry) is list:
            print(f"{spaces}[")
            print_log_list(entry, indent + 2)
            print(f"{spaces}]")
        elif type(entry) is dict:
            print(f"{spaces}{{")
            print_log_dict(entry, indent + 2)
            print(f"{spaces}}}")
        else:
            print(f"{spaces}{entry}")


def print_log_dict(a_dict, indent=0):
    """Print dictionary value in step attempt dictionary"""
    spaces = (indent + 2) * " "
    for key in a_dict:
        if type(a_dict[key]) is list:
            print(f"{spaces}{key} :")
            print(f"{spaces}[")
            print_log_list(a_dict[key], indent=indent + 2)
            print(f"{spaces}]")
        elif type(a_dict[key]) is dict:
            print(f"{spaces}{key} :")
            print(f"{spaces}{{")
            print_log_dict(a_dict[key], indent=indent + 2)
            print(f"{spaces}}}")
        else:
            print(f"{spaces}{key} : {a_dict[key]}")


def print_log(log, indent=0):
    """
    Print the entries from a log file list of step attempts.
    """
    spaces = indent * " "
    subspaces = (indent + 2) * " "
    for entry in log:
        print(f"{spaces}{{")
        for key in entry:
            if type(entry[key]) is list:
                print(f"{subspaces}{key} :")
                print(f"{subspaces}[")
                print_log_list(entry[key], indent=indent + 2)
                print(f"{subspaces}]")
            elif type(entry[key]) is dict:
                print(f"{subspaces}{key} :")
                print(f"{subspaces}{{")
                print_log_dict(entry[key], indent=indent + 2)
                print(f"{subspaces}}}")
            else:
                print(f"{subspaces}{key} : {entry[key]}")
        print(f"{spaces}}}")


def get_history(log, key, step_status=None, time_range=None, step_range=None):
    """
    Extract the step/time series of the requested value.
    """

    steps = []
    times = []
    values = []

    for entry in log:

        step = np.longlong(entry["step"])
        time = np.double(entry["t_n"])

        if time_range is not None:
            if time < time_range[0] or time > time_range[1]:
                continue

        if step_range is not None:
            if step < step_range[0] or step > step_range[1]:
                continue

        if step_status is not None:
            if step_status not in entry["status"]:
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
                    sub_steps, sub_times, sub_values = get_history(
                        s[next_level_key], key
                    )
                    steps.extend(sub_steps)
                    times.extend(sub_times)
                    values.extend(sub_values)

    return steps, times, values

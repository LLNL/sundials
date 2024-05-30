#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Programmer(s): Cody Balos @ LLNL
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
    pattern = re.compile("\[(\w+)\]\[(rank \d+)\]\[(.*)\]\[(.*)\](.*)")
    matches = pattern.findall(line)
    line_dict = {}
    if matches:
        line_dict["loglvl"] = matches[0][0]
        line_dict["rank"] = matches[0][1]
        line_dict["scope"] = matches[0][2]
        line_dict["label"] = matches[0][3]
        line_dict["payload"] = parse_logfile_payload(
            matches[0][4], line_number, all_lines
        )
    return line_dict


def log_file_to_list(filename, step_scope_txt):
    """
    This function takes a debug log file from CVODE and creates a list where
    each list entry is a step attempt.

    E.g.,
    [
      [
          {
              "loglvl": "DEBUG",
              "rank": "rank 0",
              "scope": "<step_scope_txt>",
              "label": "enter-step-attempt-loop",
              "payload": {"step": "0", "h": "1e-06", "q": "1", "t_n": "0"},
          }, ...
      ], ...
    ]
    """
    with open(filename, "r") as logfile:
        log = []
        lines_for_this_step = None
        all_lines = logfile.readlines()
        for line_number, line in enumerate(all_lines):
            line_dict = parse_logfile_line(line.rstrip(), line_number, all_lines)
            if not line_dict:
                continue
            if (
                line_dict["scope"] == step_scope_txt
                and line_dict["label"] == "enter-step-attempt-loop"
            ):
                if lines_for_this_step is None:
                    lines_for_this_step = [line_dict]
                else:
                    log.append(lines_for_this_step)
                    lines_for_this_step = [line_dict]
            else:
                lines_for_this_step.append(line_dict)
    return log


def cvode_debug_file_to_list(filename):
    """
    This function takes a debug log file from CVODE and creates a list where
    each list entry is a step attempt. See log_file_to_list.
    """
    return log_file_to_list(filename, "CVODE::cvStep")

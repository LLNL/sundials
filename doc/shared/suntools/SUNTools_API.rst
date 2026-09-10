..
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2025-2026, Lawrence Livermore National Security,
   University of Maryland Baltimore County, and the SUNDIALS contributors.
   Copyright (c) 2013-2025, Lawrence Livermore National Security
   and Southern Methodist University.
   Copyright (c) 2002-2013, Lawrence Livermore National Security.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _SUNTools.API:

suntools API
============

The API documentation below is extracted from the RST-formatted docstrings in
the ``suntools`` source.  The parser functions return ordinary Python
dictionaries and lists, so their output can be passed directly to analysis or
plotting packages.

Log parsing
-----------

.. autofunction:: suntools.logs.log_file_to_list

.. autofunction:: suntools.logs.print_log

.. autofunction:: suntools.logs.get_history

.. autoclass:: suntools.logs.StepData
   :members:

CSV and table parsing
---------------------

.. autofunction:: suntools.csv.keys

.. autofunction:: suntools.csv.read

.. autofunction:: suntools.csv.write

.. autofunction:: suntools.table.parse_table

.. autofunction:: suntools.utils.str2num

Command-line API
----------------

.. autofunction:: suntools.cli.build_parser

.. autofunction:: suntools.cli.main

.. autofunction:: suntools.tune.cli.run_from_args

Tuning API
----------

The ``suntools.tune`` package provides a backend-neutral configuration and
trial runner for tuning SUNDIALS ``SetOptions`` parameters.  DeepHyper,
GPTune, and Ytopt are supported through separate backend modules.

Configuration models
--------------------

.. autoclass:: suntools.tune.ParameterSpec
   :members:

.. autoclass:: suntools.tune.BackendConfig
   :members:

.. autoclass:: suntools.tune.SearchConfig
   :members:

.. autoclass:: suntools.tune.ExecutableConfig
   :members:

.. autoclass:: suntools.tune.MetricConfig
   :members:

.. autoclass:: suntools.tune.ObjectiveConfig
   :members:

.. autoclass:: suntools.tune.ConstraintConfig
   :members:

.. autoclass:: suntools.tune.TuneConfig
   :members:

Configuration helpers
---------------------

.. autofunction:: suntools.tune.config.parse_parameter_spec

.. autofunction:: suntools.tune.config.parse_key_value

.. autofunction:: suntools.tune.config.load_config

.. autofunction:: suntools.tune.config.parse_regex_group

.. autofunction:: suntools.tune.config.config_from_args

Trial runner
------------

.. autoclass:: suntools.tune.runner.TrialResult
   :members:

.. autofunction:: suntools.tune.runner.expand_environment_variables

.. autofunction:: suntools.tune.runner.build_trial_argv

.. autofunction:: suntools.tune.runner.run_trial

.. autofunction:: suntools.tune.runner.run_baseline

.. autofunction:: suntools.tune.runner.run_trial_async

.. autofunction:: suntools.tune.runner.extract_objective

.. autofunction:: suntools.tune.runner.extract_constraint

.. autofunction:: suntools.tune.runner.objective_to_score

.. autofunction:: suntools.tune.runner.select_best

.. autofunction:: suntools.tune.runner.select_worst

.. autofunction:: suntools.tune.runner.write_results

.. autofunction:: suntools.tune.runner.format_command

Backend adapters
----------------

.. autofunction:: suntools.tune.deephyper_backend.to_deephyper_problem

.. autoclass:: suntools.tune.deephyper_backend.DeepHyperBackend
   :members:

.. autofunction:: suntools.tune.gptune_backend.to_gptune_problem

.. autoclass:: suntools.tune.gptune_backend.GPTuneBackend
   :members:

.. autofunction:: suntools.tune.ytopt_backend.to_ytopt_problem

.. autoclass:: suntools.tune.ytopt_backend.YtoptBackend
   :members:

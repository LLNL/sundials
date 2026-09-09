#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025-2026, Lawrence Livermore National Security,
# University of Maryland Baltimore County, and the SUNDIALS contributors.
# Copyright (c) 2013-2025, Lawrence Livermore National Security
# and Southern Methodist University.
# Copyright (c) 2002-2013, Lawrence Livermore National Security.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------------------

"""Parse command-line and YAML configurations for :mod:`suntools.tune`."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import yaml

from suntools.tune.models import (
    BackendConfig,
    ConstraintConfig,
    ExecutableConfig,
    ObjectiveConfig,
    ParameterSpec,
    SearchConfig,
    TuneConfig,
)


def parse_parameter_spec(name: str, spec: str) -> ParameterSpec:
    """Parse a command-line parameter specification.

    :param str name: ``SetOptions`` key for the parameter.
    :param str spec: Specification in one of the forms ``LOW:HIGH``,
                     ``LOW:HIGH:log``, ``int:LOW:HIGH``, or
                     ``choice:v1,v2,v3``.
    :returns: Validated backend-independent parameter model.
    :rtype: ParameterSpec
    :raises ValueError: If ``spec`` does not use a supported form.
    """

    if spec.startswith("choice:"):
        values = spec[len("choice:") :].split(",")
        return ParameterSpec(name=name, type="choice", values=values)

    if spec.startswith("int:"):
        parts = spec.split(":")
        if len(parts) != 3:
            raise ValueError("int parameter specs must use int:LOW:HIGH")
        lower = int(parts[1])
        upper = int(parts[2])
        return ParameterSpec(name=name, type="int", bounds=(lower, upper))

    parts = spec.split(":")
    if len(parts) not in (2, 3):
        raise ValueError("float parameter specs must use LOW:HIGH or LOW:HIGH:log")
    scale = "linear"
    if len(parts) == 3:
        if parts[2] != "log":
            raise ValueError("the only supported float scale suffix is :log")
        scale = "log"
    lower = float(parts[0])
    upper = float(parts[1])
    return ParameterSpec(name=name, type="float", bounds=(lower, upper), scale=scale)


def parse_key_value(items: Optional[Iterable[str]], option_name: str) -> Dict[str, str]:
    """Parse repeated ``KEY=VALUE`` command-line options.

    :param items: Option values to parse, or ``None``.
    :param str option_name: Name used in validation error messages.
    :returns: A dictionary containing the parsed key-value pairs.
    :rtype: dict[str, str]
    :raises ValueError: If an item does not contain a non-empty key and ``=``.
    """
    result: Dict[str, str] = {}
    if not items:
        return result
    for item in items:
        key, separator, value = item.partition("=")
        if not separator or not key:
            raise ValueError("%s entries must use KEY=VALUE" % option_name)
        result[key] = value
    return result


def parse_regex_group(value: Any) -> Any:
    """Convert a numeric regex group string to an integer.

    :param value: Group index or named group.
    :returns: An integer for digit-only strings; otherwise ``value`` unchanged.
    """
    if isinstance(value, str) and value.isdigit():
        return int(value)
    return value


def load_config(path: str) -> TuneConfig:
    """Load and validate a YAML tuning configuration.

    :param str path: YAML configuration path.
    :returns: Validated configuration with relative paths resolved against the
              YAML file's directory.
    :rtype: TuneConfig
    :raises ValueError: If the file is empty or fails model validation.
    """
    config_path = Path(path)
    with config_path.open("r") as fp:
        data = yaml.safe_load(fp)
    if data is None:
        raise ValueError("empty tune configuration")
    config = TuneConfig.model_validate(data)
    return _resolve_relative_paths(config, config_path.parent)


def _resolve_relative_paths(config: TuneConfig, base_dir: Path) -> TuneConfig:
    output_dir = config.search.output_dir
    executable_cwd = config.executable.cwd
    updates: Dict[str, Any] = {}
    if not output_dir.is_absolute():
        updates["search"] = config.search.model_copy(update={"output_dir": base_dir / output_dir})
    if not executable_cwd.is_absolute():
        updates["executable"] = config.executable.model_copy(
            update={"cwd": base_dir / executable_cwd}
        )
    if updates:
        config = config.model_copy(update=updates)
    return config


def config_from_args(args: Any) -> TuneConfig:
    """Build a :class:`TuneConfig` from parsed CLI arguments.

    When ``args.config`` is set, the YAML file is loaded and all other tuning
    fields are ignored. Otherwise the executable, parameters, objective, and
    optional constraint are assembled from the command-line namespace.

    :param args: Namespace containing the options created by
                 :func:`suntools.cli.build_parser`.
    :returns: Validated tuning configuration.
    :rtype: TuneConfig
    :raises ValueError: If required executable or parameter options are absent.
    """
    if getattr(args, "config", None):
        return load_config(args.config)

    executable: Sequence[str] = getattr(args, "executable", None) or []
    if executable and executable[0] == "--":
        executable = executable[1:]
    if not executable:
        raise ValueError("suntools tune requires an executable command")

    parameter_items: Optional[List[Tuple[str, str]]] = getattr(args, "params", None)
    if not parameter_items:
        raise ValueError("suntools tune requires at least one --params KEY SPEC")
    parameters = [parse_parameter_spec(name, spec) for name, spec in parameter_items]

    objective = ObjectiveConfig(
        metric=args.metric,
        direction=args.direction,
        source=args.objective_source,
        regex=args.objective_regex,
        group=parse_regex_group(args.objective_group),
    )

    constraint = None
    constraint_args = (
        args.constraint_metric,
        args.constraint_source,
        args.constraint_regex,
        args.constraint_upper_bound,
    )
    if any(value is not None for value in constraint_args):
        constraint = ConstraintConfig(
            metric=args.constraint_metric,
            source=args.constraint_source,
            regex=args.constraint_regex,
            group=parse_regex_group(args.constraint_group),
            upper_bound=args.constraint_upper_bound,
        )

    return TuneConfig(
        backend=BackendConfig(
            name=args.backend, options=parse_key_value(args.backend_option, "--backend-option")
        ),
        search=SearchConfig(
            max_evals=args.max_evals,
            workers=args.workers,
            repetitions=args.repetitions,
            output_dir=Path(args.output_dir),
        ),
        executable=ExecutableConfig(
            command=executable[0],
            args=list(executable[1:]),
            cwd=Path(args.cwd),
            env=parse_key_value(args.env, "--env"),
        ),
        parameters=parameters,
        objective=objective,
        constraint=constraint,
    )

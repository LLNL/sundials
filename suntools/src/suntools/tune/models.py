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

"""Data models used by the :mod:`suntools.tune` API.

The models provide a backend-independent description of an executable, its
tunable ``SetOptions`` parameters, and the metrics used to evaluate trials.
They validate both command-line and YAML configurations before a tuning
backend is started.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Tuple, Union

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

ParameterType = Literal["float", "int", "choice"]
ParameterScale = Literal["linear", "log"]
ObjectiveDirection = Literal["minimize", "maximize"]
MetricAggregation = Literal["sum", "mean"]


class ParameterSpec(BaseModel):
    """Backend-neutral representation of a tunable ``SetOptions`` parameter.

    :ivar str name: ``SetOptions`` key appended to the executable command.
    :ivar str type: Parameter kind: ``"float"``, ``"int"``, or ``"choice"``.
    :ivar tuple[float, float] bounds: Inclusive lower and upper bounds for
                                      numeric parameters.
    :ivar list[str] values: Allowed values for a choice parameter. A value may
                            contain whitespace-separated command-line arguments.
    :ivar str scale: Numeric sampling scale, either ``"linear"`` or ``"log"``.
    """

    model_config = ConfigDict(extra="forbid")

    name: str
    type: ParameterType
    bounds: Optional[Tuple[float, float]] = None
    values: Optional[List[str]] = None
    scale: ParameterScale = "linear"

    @field_validator("name")
    @classmethod
    def _name_must_not_be_empty(cls, value: str) -> str:
        if not value:
            raise ValueError("parameter name must not be empty")
        if value.startswith("-"):
            raise ValueError("parameter name must be a SetOptions key, not an option")
        return value

    @field_validator("values")
    @classmethod
    def _values_must_not_be_empty(cls, value: Optional[List[str]]) -> Optional[List[str]]:
        if value is not None and not value:
            raise ValueError("choice parameters require at least one value")
        if value is not None and any(not item.split() for item in value):
            raise ValueError("choice values must not be empty")
        return value

    @model_validator(mode="after")
    def _validate_by_type(self) -> "ParameterSpec":
        if self.type in ("float", "int"):
            if self.bounds is None:
                raise ValueError("float and int parameters require bounds")
            lower, upper = self.bounds
            if lower >= upper:
                raise ValueError("parameter lower bound must be less than upper bound")
            if self.scale == "log" and lower <= 0.0:
                raise ValueError("log-scaled parameters require a positive lower bound")
            if self.values is not None:
                raise ValueError("float and int parameters must not define values")
        elif self.type == "choice":
            if self.values is None:
                raise ValueError("choice parameters require values")
            if self.bounds is not None:
                raise ValueError("choice parameters must not define bounds")
            if self.scale != "linear":
                raise ValueError("choice parameters must use linear scale")
        return self

    def format_value(self, value: Any) -> str:
        """Return the command-line representation for a sampled value.

        :param value: Sampled value supplied by an optimization backend.
        :returns: The value converted to a command-line string. Integer
                  parameters are normalized through :class:`int` first.
        :rtype: str
        """

        if self.type == "int":
            return str(int(value))
        return str(value)

    def format_values(self, value: Any) -> List[str]:
        """Return one or more command-line values for a sampled value.

        :param value: Sampled value supplied by an optimization backend.
        :returns: One or more command-line tokens. Choice values are split on
                  whitespace; numeric parameters produce one token.
        :rtype: list[str]

        Choice values may contain whitespace-separated values when a
        SetOptions key accepts more than one argument, for example
        ``"ARKODE_DIRK_NONE ARKODE_ERK_NONE"`` for ``table_names``.
        """

        if self.type == "choice":
            return self.format_value(value).split()
        return [self.format_value(value)]


class BackendConfig(BaseModel):
    """Configuration for a tuning backend.

    :ivar str name: Backend name (``deephyper``, ``gptune``, or ``ytopt``).
    :ivar dict options: Backend-specific keyword options.
    """

    model_config = ConfigDict(extra="forbid")

    name: str = "deephyper"
    options: Dict[str, Any] = Field(default_factory=dict)


class SearchConfig(BaseModel):
    """Search budget and output settings.

    :ivar int max_evals: Maximum number of sampled configurations.
    :ivar int workers: Number of concurrent evaluator workers.
    :ivar int repetitions: Number of executions for each configuration.
    :ivar pathlib.Path output_dir: Directory where tuning results are written.
    """

    model_config = ConfigDict(extra="forbid")

    max_evals: int = Field(default=40, gt=0)
    workers: int = Field(default=1, gt=0)
    repetitions: int = Field(default=1, gt=0)
    output_dir: Path = Path("suntools-tune")


class ExecutableConfig(BaseModel):
    """Executable command and process environment for each trial.

    :ivar str command: Executable name or path.
    :ivar list[str] args: Arguments placed before sampled ``SetOptions`` pairs.
    :ivar pathlib.Path cwd: Working directory used to resolve the command.
    :ivar dict[str, str] env: Environment variable overrides.
    """

    model_config = ConfigDict(extra="forbid")

    command: str
    args: List[str] = Field(default_factory=list)
    cwd: Path = Path(".")
    env: Dict[str, str] = Field(default_factory=dict)

    @field_validator("command")
    @classmethod
    def _command_must_not_be_empty(cls, value: str) -> str:
        if not value:
            raise ValueError("executable command must not be empty")
        return value


class MetricConfig(BaseModel):
    """Metric extraction settings shared by objectives and constraints.

    :ivar str metric: Name of the metric in reports and result files.
    :ivar str source: ``stdout``, ``stderr``, or a file path.
    :ivar str or list[str] regex: Pattern(s) used to extract numeric values.
    :ivar int or str group: Match group containing the numeric value.
    :ivar str aggregation: Aggregation for multiple patterns (``sum`` or
                           ``mean``).
    """

    model_config = ConfigDict(extra="forbid")

    metric: str
    source: Optional[str] = None
    regex: Optional[Union[str, List[str]]] = None
    group: Union[int, str] = 1
    aggregation: MetricAggregation = "sum"

    @field_validator("metric")
    @classmethod
    def _metric_must_not_be_empty(cls, value: str) -> str:
        if not value:
            raise ValueError("metric must not be empty")
        return value

    @field_validator("regex")
    @classmethod
    def _regex_must_not_be_empty(
        cls, value: Optional[Union[str, List[str]]]
    ) -> Optional[Union[str, List[str]]]:
        if isinstance(value, list) and (not value or any(not item for item in value)):
            raise ValueError("metric regex values must not be empty")
        return value

    @model_validator(mode="after")
    def _set_default_source(self) -> "MetricConfig":
        if self.regex and self.source is None:
            self.source = "stdout"
        if not self.regex and self.metric != "wall_time":
            raise ValueError("non-wall_time metrics require regex")
        return self


class ObjectiveConfig(MetricConfig):
    """Objective metric and optimization direction."""

    metric: str = "wall_time"
    direction: ObjectiveDirection = "minimize"


class ConstraintConfig(MetricConfig):
    """An upper-bound constraint on a metric extracted from trial output."""

    upper_bound: float

    @model_validator(mode="after")
    def _validate_constraint(self) -> "ConstraintConfig":
        if not self.regex:
            raise ValueError("constraints require a metric regex")
        return self


class TuneConfig(BaseModel):
    """Complete backend-independent tuning configuration.

    :ivar BackendConfig backend: Optimization backend and backend options.
    :ivar SearchConfig search: Evaluation budget and output directory.
    :ivar ExecutableConfig executable: Command run for each trial.
    :ivar list[ParameterSpec] parameters: Tunable ``SetOptions`` parameters.
    :ivar ObjectiveConfig objective: Metric to optimize.
    :ivar ConstraintConfig constraint: Optional upper-bound metric constraint.
    """

    model_config = ConfigDict(extra="forbid")

    backend: BackendConfig = Field(default_factory=BackendConfig)
    search: SearchConfig = Field(default_factory=SearchConfig)
    executable: ExecutableConfig
    parameters: List[ParameterSpec]
    objective: ObjectiveConfig = Field(default_factory=ObjectiveConfig)
    constraint: Optional[ConstraintConfig] = None

    @field_validator("parameters")
    @classmethod
    def _parameters_must_be_unique(cls, value: List[ParameterSpec]) -> List[ParameterSpec]:
        if not value:
            raise ValueError("at least one parameter is required")
        names = [parameter.name for parameter in value]
        duplicates = sorted({name for name in names if names.count(name) > 1})
        if duplicates:
            raise ValueError("parameter names must be unique: " + ", ".join(duplicates))
        return value

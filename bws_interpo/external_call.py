import os
import sys
import re
from dataclasses import dataclass
from subprocess import Popen, PIPE

import ast
import inspect
from typing import Dict, Any, Callable, Set, List

import itertools
import functools
from functools import partial

from .func_params import extract_func_params
from .cmd_exec import ShellExecProvider


@dataclass
class Domain:
    name: str
    command: List[str]


def _func_is_empty(func):

    source = inspect.getsource(func)

    tree = ast.parse(source)

    # Find the function definition in the AST
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef):
            body = node.body

            return len(body) == 1 and isinstance(body[0], ast.Pass)

    return False  # In case no FunctionDef was found


def _param_to_cmd_args(
    param: inspect.Parameter, value: Any, args_map, boolean_params: Set[str]
):

    k = param.name
    param = param.replace(name=args_map[k]) if k in args_map else param
    boolean_params = {args_map[p]
                      if p in args_map else p for p in boolean_params}

    if param.name in boolean_params:
        return (
            [f"--{param.name}"] if value else []
        )  # Use `if value` to support implicit booleaness of Python
    else:
        is_keyword = param.kind == inspect.Parameter.KEYWORD_ONLY
        prefix = "--" if is_keyword else "-"

        return [prefix + param.name, str(value)]


def foreign_function(
    f,
    domain: Domain,
    func_name=None,
    args_map=None,
    args_validation=None,
    postprocess_fn=None,
    **run_args,
):
    """
    Use this decorator to expose XMIPP programs to Python.

    XMIPP programs are called in the command line using scipion, for example
    `scipion run xmipp_xmipp_volume_from_pdb -i ... -o ...`. The
    `foreign_function` decorator sets up wiring to transform a Python function
    call into command line arguments.

    When exposing an XMIPP program you can rename input arguments to make the
    function more Pythonic and add input validation. See Exposing XMIPP Programs
    for more details.
    """

    is_empty = _func_is_empty(f)
    if not is_empty:
        raise RuntimeError(
            f"Forward declared external scipion function {f.__name__} must be only contain a single pass statement."
        )

    if args_map is None:
        args_map = {}
    if args_validation is None:
        args_validation = {}

    func_name = func_name if func_name is not None else f.__name__

    params = inspect.signature(f).parameters
    boolean_params = {k for k, v in f.__annotations__.items() if v is bool}

    pos_args = {
        k
        for k, v in params.items()
        if v.kind == inspect.Parameter.POSITIONAL_OR_KEYWORD
    }

    if boolean_params.intersection(pos_args):
        raise RuntimeError(
            "Positional arguments cannot be declared as boolean flags")

    run_args.setdefault("shell", True)
    # run_args["stdout"]=PIPE
    run_args["stderr"] = PIPE

    args_validation = {k: re.compile(v) for k, v in args_validation.items()}

    @functools.wraps(f)
    # @inject
    def wrapper(
        *args,
        # __scipion_bridge_runner__: ShellExecProvider = Provide[Container.shell_exec],
        **kwargs,
    ):
        _ = f(
            *args, **kwargs
        )  # Call function for Python to throw error if args and kwargs aren't passed correctly

        merged_args = extract_func_params(args, kwargs, params)

        # Filter args that are None to support optional arguments
        merged_args = {k: v for k, v in merged_args.items() if v is not None}

        # Validate inputs before calling external program
        # arg_names = {k.name for k in merged_args}
        for arg, value in merged_args.items():
            if not arg.name in args_validation:
                continue

            pattern = args_validation[arg.name]
            if not re.fullmatch(pattern, value):
                raise ValueError(
                    f"Value '{value}' for does not have the required format for '{arg.name}'"
                )

        raw_args = [
            _param_to_cmd_args(p, v, args_map, boolean_params)
            for p, v in merged_args.items()
        ]
        if postprocess_fn is not None:
            raw_args = postprocess_fn(raw_args)

        raw_args = itertools.chain.from_iterable(raw_args)
        raw_args = domain.command + [func_name, *raw_args]

        runner = ShellExecProvider()
        return runner(func_name, domain, raw_args, run_args)

    return wrapper

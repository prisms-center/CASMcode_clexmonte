import inspect
import typing

import numpy as np


def to_dict(
    value: typing.Any, data: dict, option: str, write_null: bool = False, **kwargs
):
    if value is not None:
        methods = inspect.getmembers(value.__class__, inspect.isfunction)
        if "to_dict" in methods:
            data[option] = value.to_dict(data.get(option), **kwargs)
        else:
            data[option] = value
    elif write_null:
        data[option] = None


def required_from_dict(required_type: typing.Any, data: dict, option: str, **kwargs):
    if option not in data:
        raise Exception(
            f"Error parsing dict: missing required '{option}'"
            f"of type '{required_type.__name__}'"
        )
    try:
        methods = inspect.getmembers(required_type, inspect.isfunction)
        if "from_dict" in methods:
            value = required_type.from_dict(data.get(option), **kwargs)
        else:
            value = required_type(data.get(option))
    except Exception as e:
        print("what:", e)
        raise Exception(
            f"Error parsing dict: failed converting required '{option}'"
            f"to type '{required_type.__name__}'"
        )
    return value


def required_int_array_from_dict(data: dict, option: str):
    if option not in data:
        raise Exception(
            f"Error parsing dict: missing required '{option}'"
            f"of integer array-like type"
        )
    try:
        value = np.array(data.get(option), dtype="int")
    except Exception as e:
        print("what:", e)
        raise Exception(
            f"Error parsing dict: failed converting required '{option}'"
            f"to integer array"
        )
    return value


def required_array_from_dict(data: dict, option: str):
    if option not in data:
        raise Exception(
            f"Error parsing dict: missing required '{option}'" f"of array-like type"
        )
    try:
        value = np.array(data.get(option))
    except Exception as e:
        print("what:", e)
        raise Exception(
            f"Error parsing dict: failed converting required '{option}' to array"
        )
    return value


def optional_from_dict(
    required_type: typing.Any,
    data: dict,
    option: str,
    default_value: typing.Any = None,
    **kwargs,
):
    value = data.get(option)
    if value is not None:
        methods = inspect.getmembers(required_type, inspect.isfunction)
        if "from_dict" in methods:
            return required_type.from_dict(value, **kwargs)
        else:
            return required_type(value)
    return default_value

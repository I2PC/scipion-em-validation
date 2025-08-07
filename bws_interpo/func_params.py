from collections import OrderedDict


def extract_func_params(args, kwargs, params):
    defaults = {
        k: param for k, param in params.items() if param.default is not param.empty
    }

    args_dict = dict(zip(params.values(), args))
    kwargs_dict = {params[k]: v for k, v in kwargs.items()}

    merged_args = {
        **args_dict,
        **kwargs_dict,
    }

    for v in defaults.values():
        merged_args.setdefault(v, v.default)

    return OrderedDict((k, merged_args[k]) for k in params.values())

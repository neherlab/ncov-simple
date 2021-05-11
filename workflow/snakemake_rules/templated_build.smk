'''
This file generated the build configurations for the templated builds
'''

from itertools import product

patterns = config["templated-builds"]["build_patterns"]
subsamples = config["templated-builds"]["subsamples"]

for build_vars in product(*[x.items() for x in patterns.values()]):
    build_name_params = {k:v[0] for k,v in zip(patterns.keys(), build_vars)}
    build_params = {k:v[1] for k,v in zip(patterns.keys(), build_vars)}
    build_params.update({k:eval(v.format(**build_params)) for k,v in config['templated-builds'].get('subsampling_parameters',{}).items()})

    build_name = config["templated-builds"]["build_name"].format(**build_name_params)

    tmp = {}
    for subsample in subsamples:
        tmp[subsample] = {}
        tmp[subsample]["filters"] = subsamples[subsample]["filters"].format(**build_params)
        if "priorities" in subsamples[subsample]:
            tmp[subsample]["priorities"] = subsamples[subsample]["priorities"].format(**build_params)

    config['builds'][build_name] = {'subsamples': tmp}

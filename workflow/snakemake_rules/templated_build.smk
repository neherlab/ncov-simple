'''
This file generated the build configurations for the templated builds
'''

from itertools import product

for templated_build in config["templated-builds"].values():
    patterns = templated_build["build_patterns"]
    subsamples = templated_build["subsamples"]
    metadata_adjustments = templated_build.get("metadata_adjustments",{})
    auspice_config = templated_build.get("auspice_config",{})
    description = templated_build.get("description",{})

    for build_vars in product(*[x.items() for x in patterns.values()]):
        build_name_params = {k:v[0] for k,v in zip(patterns.keys(), build_vars)}
        build_params = {k:v[1] for k,v in zip(patterns.keys(), build_vars)}
        build_params.update({k:eval(v.format(**build_params)) for k,v in templated_build.get('subsampling_parameters',{}).items()})

        build_name = templated_build["build_name"].format(**build_name_params)

        tmp = {}
        for subsample in subsamples:
            tmp[subsample] = {}
            tmp[subsample]["filters"] = subsamples[subsample]["filters"].format(**build_params)
            if "priorities" in subsamples[subsample]:
                tmp[subsample]["priorities"] = subsamples[subsample]["priorities"].format(**build_params)
        config['builds'][build_name] = {'subsamples': tmp}

        tmp = []
        for adjustment in metadata_adjustments:
            tmp.append({"query": adjustment["query"].format(**build_params),
                        "src": adjustment["src"], "dst": adjustment["dst"]})
        config['builds'][build_name]['metadata_adjustments'] = tmp

        if(auspice_config != {}):
            config['builds'][build_name]['auspice_config'] = auspice_config

        if(description != {}):
            config['builds'][build_name]['description'] = description

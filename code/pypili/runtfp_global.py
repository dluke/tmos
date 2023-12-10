

# universal default setup
import parameters
import range_tfp

def setup():
    """
    """
    params = parameters.parse_paramdata(parameters.paramdata)
    args = parameters.ParameterList(params)
    # add whatever global configuration we want here
    # ...
    args.Pili.force_max_length = True
    args.Pili.max_length = 19.99
    return args


def exit_condition(cell):
    return False

def single_run(args):
    #import ctfp3d
    #ctfp3d.pyrun(args, None)

    singlepool = range_tfp.buildsinglepool(args, exit_condition)
    cores = 1
    range_tfp.runpool(singlepool, args, cores)




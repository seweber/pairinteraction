"""pipy - Useful scripts for calculating energy spectrum over some parameter, can also use parallelization

This should be called from the gui calls, and should in some way accept a list,... where to save the results to,
such that the gui can start plotting asynchronously.

The main function is calc_list, which will do all the simulations and save the results.
Should get a dict settings, which it will use to do the simulations.
settings should have: config (config for Atom) dict or Atom object or path to pickled Atom object
scriptoptions : listoptions : dict with the parameters for the simulations
    mainly should give a list of parameters (including distance, fields and symmetries) for the simulations
TODO symmetries should be handled extra

runtimeoptions : how many parallel processes,...

TODO distance change is easy (no basis change)
but field changes might change basis? probably should brutforce fix a "good" basis in advance applied for all fields

symmetries do for sure change basis, but can fix finite set of basis in advance
"""
import concurrent.futures
import copy
import logging
import os
import pickle
import time
from functools import partial

import numpy as np

from pipy import atom_from_config

logger = logging.getLogger(__name__)


def calc_list(settings):
    start_time = time.perf_counter()
    executor = settings.pop("ProcessPoolExecutor", None)
    if executor is None:
        num_pr = settings.get("runtimeoptions", {}).get("NUM_PROCESSES", 1)
        num_pr = os.cpu_count() if num_pr in [0, -1] else num_pr
        if num_pr > 1:
            executor = concurrent.futures.ProcessPoolExecutor(num_pr)

    param_list = get_param_list(settings)
    ip_list = list(range(len(param_list)))

    p_one_run = partial(one_run, settings, param_list)

    if executor is not None:
        with executor:
            results_list = list(executor.map(p_one_run, ip_list))
    else:
        results_list = []
        for i in ip_list:
            results_list.append(p_one_run(i))

    data = {"settings": settings, "results_list": results_list}
    time_needed = time.perf_counter() - start_time
    logger.info("Time needed for calc_list: %ss", time_needed)
    return data


def get_param_list(settings):
    scriptoptions = settings.setdefault("scriptoptions", {})
    if "param_list" in scriptoptions:
        return scriptoptions["param_list"]
    elif "listoptions" in scriptoptions and "steps" in scriptoptions["listoptions"]:
        listoptions = scriptoptions["listoptions"]
        steps = listoptions["steps"]
        k_lists = {}
        for k in ["Bx", "By", "Bz", "Ex", "Ey", "Ez", "distance"]:
            if listoptions.get("min" + k) is not None and listoptions.get("max" + k, None) is not None:
                k_lists[k] = np.linspace(listoptions["min" + k], listoptions["max" + k], steps)
        param_list = [{k: v[i] for k, v in k_lists.items()} for i in range(steps)]
        scriptoptions["param_list"] = param_list
    else:
        raise NotImplementedError("TODO: create param list from individual lists, including symmetries")
    return param_list


def one_run(settings, param_list, ip):
    config = settings["config"]

    # get default Atom from config (either given, load it or newly construct it)
    if "atom" in config:
        atom = config["atom"]
    elif "atom_path" in config:
        atom = load_atom(config["atom_path"])
    else:
        config = copy.deepcopy(config)
        atom = atom_from_config(config)
    atom.system  # make sure system is initialized and updateFromParams will have no effect on the basisStates
    # TODO it actually should have been initialized already before parallelization

    atom.updateFromParams(param_list[ip])
    atom.calcEnergies()

    if "save_results_list" in settings:
        settings["save_results_list"][ip] = atom
    if "save_to_queue" in settings:
        settings["save_to_queue"].put((ip, atom))

    return atom


def load_atom(path):
    if not os.path.splitext(path)[1] == ".pkl":
        raise NotImplementedError("TODO: load atom from other formats")
    if not os.path.isfile(path):
        raise OSError("File not found: " + path)
    with open(path, "rb") as f:
        atom = pickle.load(f)
    return atom

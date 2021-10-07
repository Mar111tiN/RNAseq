import os
from subprocess import check_call as run
from datetime import datetime as dt
import pandas as pd

ansii_colors = {
    "magenta": "[1;35;2m",
    "green": "[1;9;2m",
    "red": "[1;31;1m",
    "cyan": "[1;36;1m",
    "gray": "[1;30;1m",
    "black": "[0m",
}

colors = {
    "process": ansii_colors["green"],
    "time": ansii_colors["magenta"],
    "normal": ansii_colors["gray"],
    "warning": ansii_colors["red"],
    "success": ansii_colors["cyan"],
}


def show_output(text, color="normal", multi=False, time=True, **kwargs):
    """
    get colored output to the terminal
    """
    time = (
        f"\033{colors['time']}{dt.now().strftime('%H:%M:%S')}\033[0m : " if time else ""
    )
    proc = f"\033{colors['process']}Process {os.getpid()}\033[0m : " if multi else ""
    text = f"\033{colors[color]}{text}\033[0m"
    print(time + proc + text, **kwargs)


def show_command(command, list=False, multi=True, **kwargs):
    """
    prints the command line if debugging is active
    """

    if list:
        command = f"\033[1m$ {' '.join(command)}\033[0m"
    else:
        command = f"\033[1m$ {command}\033[0m"
    print(proc + command, **kwargs)
    return


def run_cmd(cmd, multi=False):
    show_command(cmd, multi=multi)
    exit = run(cmd, shell=True)


def sort_df(df, cols={"Chr": True, "Start": True}):
    """
    helper for sorting dfs for chromosomes using Chr, Start + cols in cols
    """
    # make Chr column categorical for sorting .. and sort
    chrom_list = [f"chr{i}" for i in range(23)] + ["chrX", "chrY"]

    df["Chr"] = pd.Categorical(df["Chr"], chrom_list)
    return df.sort_values(list(cols.keys()), ascending=list(cols.values()))

def add_config(config, path="", config_name=""):
    '''
    update the config file with a config either 
    - from a path relative to snakedir
    - from a config name listed in the configs list
    '''
    if config_name:
        config_file = os.path.join(snakedir, config['configs'][config_name])
    else:
        config_file = os.path.join(snakedir, path)

    with open(config_file, "r") as stream:
        added_config = load(stream, Loader=Loader)
        config.update(added_config)
    return config
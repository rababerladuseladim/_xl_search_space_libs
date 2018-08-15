import re


def pretty_up_sample_name(smpl_str):
    """
    turn 'fr_03_to_06' in 'fraction 3 to 6'

    :param smpl_str:
    :return:
    """
    lst_nos = re.findall(r"(\d+)", smpl_str)
    tit = "fraction " + str(int(lst_nos.pop(0)))
    if lst_nos:
        tit += " to " + str(int(lst_nos.pop()))
    return tit

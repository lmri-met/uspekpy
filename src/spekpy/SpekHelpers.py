# For compatibility with Python2 #
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info[0] < 3:
    BYTES=str
    STR = unicode
else:
    BYTES = bytes
    STR = str
##################################
from collections import OrderedDict
from datetime import datetime as dt
import time as tm


def get_current_time_stamp(mode='str'):
    """
    A helper function to get the current time stamp

    :return str : The current time stamp
    """
    if mode == 'str':
        time_stamp = dt.fromtimestamp(tm.time()).strftime('%Y-%m-%d %H:%M:%S')
    elif mode == 'file':
        time_stamp = dt.fromtimestamp(tm.time()).strftime('%Y-%m-%d_%H-%M-%S')
    return time_stamp


def format_parameter_str(parameter_name, parameter_value, parameter_unit, 
                         decimal_places=1):
    """
    A helper function to help format spekpy parameters into a string

    :param str parameter_name: The name of the paramter
    :param int,flot,bool,None,str parameter_value: The value of the pararmeter
    :param str parameter_unit: The unit of the parameter
    :param int decimal_places: The number of decimal places to be used for 
        numeric parameter values (Optional)
    :return str : The formatted string representing a parameter
    """
    if isinstance(parameter_value, (int, float)) \
        and not isinstance(parameter_value, bool):
        the_value = ('%.*f' % (decimal_places, parameter_value))
    elif isinstance(parameter_value, bool) or parameter_value is None \
        or isinstance(parameter_value, list):
        the_value = str(parameter_value)
    elif isinstance(parameter_value, BYTES) \
        or isinstance(parameter_value, STR):
        the_value = parameter_value
    else:
        raise Exception('Could not format the parameter '+parameter_name+'!')

    return parameter_name + ': ' + the_value + ' [' + parameter_unit + ']; '


def ord_dct(dict2order):
    """
    Helper function to order dictionaries

    :param dict dict2order: The dictionary to order
    :return: The ordered dictionary
    """

    return OrderedDict(dict2order)




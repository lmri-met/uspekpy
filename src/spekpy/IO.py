# For compatibility with Python2 #
from __future__ import print_function, division, absolute_import
import sys
if sys.version_info[0] < 3:
    BYTES = str
    STR = unicode
else:
    BYTES = bytes
    STR = str
##################################
import sys, os
from collections import OrderedDict
import json
import numpy as np
import spekpy.SpekConstants as Const

PATH = os.path.dirname(os.path.abspath(__file__))
SLASH = os.sep


def full_file(*args):
    """
    A function to form an absolute path to a file from a variable number of 
    parts that represent a path relative to the root directory

    :param str args: Variable number of parts
    :return str : The full path to a file
    """
    if all(isinstance(arg, (STR, BYTES)) for arg in args):
        path_file = PATH + SLASH + SLASH.join(args)
        return path_file
    else:
        raise Exception('full_file: Not all arguments were strings.')


def find_file(file_name=None, extension=None, directories=None):
    """
    A function that locates a file of a given extension in a given directory

    :param str file_name: The name of the file with or without extension
    :param str extension: The extension of the file
    :param list directories: A list with the directories to search for the file
    :return str full_file_name: The location of the file if the file exists
    """
    if file_name and extension and directories:
        if os.path.isfile(file_name):
            return file_name
        else:
            if file_name.endswith(extension):
                file_name = file_name
            else:
                file_name = file_name + extension

            for d in directories:
                full_file_name = full_file(Const.dir_data, d, file_name)
                if os.path.isfile(full_file_name):
                    return full_file_name
                    break

    else:
        raise Exception('File name, extension, or directory not provided')

def path_file(file_name):
    """
    A function to break down a path name into parts 

    :param str file_name: The name of the file with or without extension.
    :return str path: The path
    :return str file_name: The file's name (without path)
    :return str SLASH: The OS dependent slash delimeter
    """
    path = os.path.dirname(file_name)
    file_name = os.path.basename(file_name)
    return path, file_name, SLASH


def write_json_to_disk(data, file_name, ordereddict=1):
    """
    A function to write data to disk as a json file

    :param data: Data to save 
    :param str file_name: Name of file to create
    :param int ordereddict: If 1, the data is converted to ordered dictionary 
        before saving
    """
    if ordereddict == 1:
        jsonFormat = json.dumps(OrderedDict(data), indent=4, 
                                separators=(',', ': '))
    else:
        jsonFormat = json.dumps(data, indent=4, separators=(',', ': '))
    thedir=os.path.dirname(file_name)
    if thedir: # True is string non-empty
        if not os.path.exists(os.path.dirname(file_name)):
            os.makedirs(os.path.dirname(file_name))
    with open(file_name, 'w+') as outfile:
        outfile.write(jsonFormat)


def read_json_from_disk(file_name):
    """
    A function to read data from a json file

    :param str file_name: Name of file to create
    :return dict or list data: Data from json file
    """
    with open(file_name, 'r') as file:
         try:
            data = json.load(file)
            return data
         except:
             raise Exception('\nError in loading state\nThe file "' 
                        + file_name + '" could not be interpreted with json')


def write_spectrum_to_disk(header, data, file_name, delimiter, fmt):
    """
    A function to write a spectrum to file.

    :param str header: The header of the export
    :param numpy.array data: An array containing energy and the spectrum 
    :param str file_name: The file name of the export
    :param str delim: The character to use to delimit the data
    :param str fmt: The format of the numerical values in the export
    :return:
    """
    np.savetxt(file_name, np.transpose(data), delimiter=delimiter, fmt=fmt,
               header=header)

def read_spectrum_from_disk(spekpy_obj):
    """
    A function to read a spectrum from a file

    This function returns: 
    k: Energy bins of the spectrum
    brem: Bremstrahlung contribution of the spectrum differential in energy
    char: Characteristic emission of the spectrum differential in energy
    dk: The bin width of the energy bins

    :param Spek spekpy_obj: A spekpy state.
    :return numpy.array k, numpy.array brem, numpy.array char, float dk: A 
        spectrum from file
    
    """
    file_name = spekpy_obj.state.external_spectrum.external_spectrum
    delimiter = spekpy_obj.state.external_spectrum.ext_delim
    try:
        # Characteristic contribution specified
        k, spk, spk_ch = np.loadtxt(file_name, delimiter=delimiter, 
                                    unpack=True)
        char = spk_ch[k>1] # Remove values that are lower than 1 keV in energy
        brem = spk[k>1] - char  # Remove values lower than 1 keV in energy
    except:
        # Characteristic contribution unspecified
        k, spk = np.loadtxt(file_name, delimiter=delimiter, unpack=True)
        brem = spk[k>1]  # Remove values that are lower than 1 keV in energy
        char = np.zeros(brem.shape) # Remove values lower than 1 keV in energy
    k = k[k>1]  # Remove values that are lower than 1 keV in energy
    dk = k[1] - k[0]
    return k, brem, char, dk

def get_script_path():
    """
    A function to get the path where the Spek class is being used

    :param:
    :return str script_path: The path of the script
    """
    script_path = os.getcwd()
    return script_path

def is_file(file_name, base_path):
    """
    A function to check whether a file exists.

    :param str file_name: The name of the file to check
    :param str base_path: The path of where the file to check is
    :return str file_name_path: The name of existing file with path added
    """
    path = os.path.dirname(file_name)
    if path == '':
        file_name_path = base_path + SLASH + file_name
    else:
        file_name_path = file_name

    if os.path.isfile(file_name_path):
        return file_name_path
    else:
        raise Exception('Could not find file: ' + file_name_path)

def get_matls():
    """
    A function to create a list of all of the spekpy defined and user defined 
    material compisitions

    :return list matls_usr_lst, list matls_def_lst: 2 lists conataining the 
        available material definitions
    """
    try:
        matls_usr_dir_contents = sorted(os.listdir(full_file(Const.dir_data, 
                                                  Const.dir_matl_usr)))
    except:
        matls_usr_dir_contents = []
    try:    
        matls_def_dir_contents = sorted(os.listdir(full_file(Const.dir_data, 
                                                   Const.dir_matl_def)))
    except:
        matls_def_dir_contents = []
        
    matls_usr_lst = []
    matls_def_lst = []
    for matl_usr in matls_usr_dir_contents:
        if matl_usr.endswith('.comp'):
            matls_usr_lst.append(matl_usr[:-5])
    
    for matl_def in matls_def_dir_contents:
        if matl_def.endswith('.comp'):
            matls_def_lst.append(matl_def[:-5])

    return matls_usr_lst, matls_def_lst

def get_states():
    """
    A function to create a list of all of the spekpy defined and user defined 
    states

    :return list states_usr_lst, list states_def_lst: 2 lists conataining the 
        available material definitions
    """
    try:
        states_usr_dir_contents = sorted(os.listdir(full_file(Const.dir_data, 
                                                   Const.dir_state_usr)))
    except:
        states_usr_dir_contents = []
    
    try:
        states_def_dir_contents = sorted(os.listdir(full_file(Const.dir_data, 
                                                   Const.dir_state_def)))
    except:
        states_def_dir_contents = []
        
    states_usr_lst = []
    states_def_lst = []
    for state_usr in states_usr_dir_contents:
        if state_usr.endswith('.state'):
            states_usr_lst.append(state_usr[:-6])
    
    for state_def in states_def_dir_contents:
        if state_def.endswith('.state'):
            states_def_lst.append(state_def[:-6])

    return states_usr_lst, states_def_lst

def delete_file(file_name_parts):
    """
    A function to delete a specified file

    :param str, file_name_parts: variable number of parts 
    :return 
    """
    file_name_with_path=full_file(*file_name_parts)
    if os.path.exists(file_name_with_path):
        os.remove(file_name_with_path)
    else:
        print("Unable to delete specified file (it does not exist)")
    return

def print_matls(matl_dir=None):
    """
    A method to print all of the existing materials to the terminal. This can 
    be used to see which materials are available in spekpy or have been created
    by the user

    :return:
    """
    matls_usr, matls_def = get_matls()
    if matl_dir == "usr" or matl_dir == None:
        print('User defined materials:')
        for mu in matls_usr:  # mu: a single material in matls_user
            print(mu)
    print()
    if matl_dir == "def" or matl_dir == None:
        print('spekpy defined materials:')
        for md in matls_def:  # md: a single material in matls_def
            print(md)

def print_states(state_dir=None):
    """
    A method to print all of the existing states to the terminal. This can be 
    used to see which states are available in spekpy or have been created by 
    the user

    :return:
    """
    states_usr, states_def = get_states()
    if state_dir == "usr" or state_dir == None:
        print('User defined states:')
        for su in states_usr:  # su: a single state in states_usr
            print(su)
    print()
    if state_dir == "def" or state_dir == None:
        print('spekpy defined states:')
        for sd in states_def:  # sd: a single states in states_def
            print(sd)

def print_matl_info(matl):
    """
    A method to print the definition of a material to the terminal

    :return:
    """
    matl_file = find_file(file_name=matl, 
                          extension=Const.extension_matl_composition, 
                          directories=[Const.dir_matl_usr, Const.dir_matl_def])
    if matl_file is not None: 
        w = 0
        matl_comp = read_json_from_disk(matl_file)
        matl_def_str = '\nMaterial Definition:\n'
        matl_def_str = matl_def_str + 'Name: ' + matl + '\n'
        try:
            matl_def_str = matl_def_str + 'Comment: ' + \
                matl_comp['file_information']['comment'] + '\n'
        except:
            pass
        matl_def_str = matl_def_str + 'Density: ' + \
            str(matl_comp['composition']['density']) + ' [g cm^-3]\n'
        matl_def_str = matl_def_str + 'Number of elements: ' + \
            str(matl_comp['composition']['number_of_elements']) + '\n'
        matl_def_str = matl_def_str + 'Z\tMass Fraction\n'
        elements = matl_comp['composition']['elements']
        for el in elements:
            w += el[1]
            matl_def_str = matl_def_str + str(el[0]) + '\t' + str(el[1]) + '\n'
        matl_def_str = matl_def_str + 'Total:\t' + str(w) + '\n'

    else:
        raise Exception('Could not find the material definition: ' + matl)
    
    print(matl_def_str)

def print_matls_in_group(matl_group):
    """ 
    A function to print all of the materials that belong to a group of 
    definitions. ICRU and ICRP are supported as material groups in this 
    release. 

    :param str matl_group: The name of the material definition group
    :return:
    """

    if matl_group != 'ICRU' and matl_group != 'ICRP':
        raise Exception(matl_group + 
                ' is not a supported material group. Use either ICRP or ICRU')

    matls_user, matls_def = get_matls()
    print('Available material definitions for the material group ' + matl_group)
    for matl_def in matls_def:
        if matl_group in matl_def:
            print(matl_def)

    
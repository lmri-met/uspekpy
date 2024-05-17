# For compatibility with Python2 #
from __future__ import print_function, division, absolute_import
##################################
from spekpy.SpekHelpers import get_current_time_stamp, format_parameter_str \
    as fmt_param, ord_dct
from spekpy.IO import write_json_to_disk, read_json_from_disk, \
    write_spectrum_to_disk, full_file
import spekpy.SpekConstants as Const


class State(object):
    """
    A class to track, import and save the state of a spekpy instance
    """

    def __init__(self):
        self.comment = None
        self.spekpy_version = None
        self.script_path = None
        self.external_spectrum = ExternalSpectrumDef()
        self.model_parameters = ModelParameterDef()
        self.spectrum_parameters = SpectrumParameterDef()
        self.filtration = FiltrationDef()
        self.flags = SpekpyFlagsDef()

    def export_spectrum_to_disk(self, file_name, delimiter, 
                                std_results, spk_char):
        """
        A method to export the current spectrum to a file on the disk

        :param str file_name: The name of the spectrum file
        :param str delimiter: The delimeter to delimit values in text file
        :return:  
        """      
        if not file_name:
            file_name = 'spekpy_export ' + get_current_time_stamp(
                mode='file') + Const.extension_export_spectrum_file

        header = 'Header Begin\n'
        header = header + self.get_basic_info_str(
            file_type='exported_spectrum') + '\n'
        header = header + self.get_current_state_str('full', std_results)
        header = header + 'Spectrum Units:\n[keV], [Photons keV^-1 cm^-2],' \
            ' [Photons keV^-1 cm^-2]' + '\n'
        header = header + 'Header End: ' + str(len(header.split('\n'))) \
            + ' lines'
        data = [std_results.k, std_results.spk, spk_char]
        print("spectrum name")
        print(file_name)
        write_spectrum_to_disk(header, data, file_name, delimiter, fmt='%.4f')

    def get_basic_info_str(self, file_type='-'):
        """
        A method to format basic info into a single string

        :param str file_type: The type of file
        :return str info_str: String containing spekpy version, filetype and 
            comment  
        """
        info_str = 'spekpy version: ' + self.spekpy_version + '\n'
        info_str = info_str + 'file type: ' + file_type + '\n'
        info_str = info_str + 'comment: ' + self.comment
        return info_str

    def get_model_parameters_str(self):
        """
        A method to format the model parameters into a single string

        :return str mp_str: A string with the model parameters
        """
        mp = self.model_parameters.__dict__
        mp_str = fmt_param('Tube Voltage', mp['kvp'], 'kVp') + \
                 fmt_param('Anode Angle', mp['th'], 'degrees') + \
                 fmt_param('Energy Bin', mp['dk'], 'keV') + \
                 fmt_param('Bin shift fraction', mp['shift'], '') + \
                 fmt_param('Physics Mode', mp['physics'], 'str') + \
                 fmt_param('Mu Data Source', mp['mu_data_source'], 'str') + \
                 fmt_param('Target', mp['targ'], 'str')

        return mp_str

    def get_spectrum_parameters_str(self):
        """
        A method to format the spectrum parameters into a single string

        :return str sp_str: A string with the spectrum parameters
        """
        sp = self.spectrum_parameters.__dict__
        sp_str = fmt_param('x', sp['x'], 'cm') + \
                 fmt_param('y', sp['y'], 'cm') + \
                 fmt_param('z', sp['z'], 'cm') + \
                 fmt_param('Tube Load', sp['mas'], 'mAs') + \
                 fmt_param('Bremsstrahlung Emission', sp['brem'], 'bool') + \
                 fmt_param('Characteristic Emission', sp['char'], 'bool') + \
                 fmt_param('Oblique', sp['obli'], 'bool') + \
                 fmt_param('Ref. air Kerma', sp['ref_kerma'], 'uGy') + \
                 fmt_param('Ref. fluence', sp['ref_flu'], 'Photons cm^-2')
        return sp_str

    def get_filtration_str(self):
        """
        A method to format the filtration details into a single string

        :return str filt_str: A string with the filtration details
        """
        filt = self.filtration.__dict__
        filt_str = fmt_param('Filtration', filt['filters'], 'mm')
        return filt_str

    def get_current_results_str(self, std_results):
        """
        A method to format the current standard results into a single string
        
        :return str result_str: A string with the standard results
        """
        res = std_results.__dict__
        result_str = fmt_param('Fluence', res['flu'], 'Photons cm^-2', 
                               decimal_places=4) + ' ' + \
                     fmt_param('Air Kerma', res['kerma'], 'uGy', 
                               decimal_places=4) + \
                     fmt_param('Mean Energy', res['emean'], 'keV', 
                               decimal_places=4) + '\n'
        result_str = result_str + fmt_param('First HVL Al', 
                               res['hvl_1_al'], 'mm Al', 
                               decimal_places=4) + ' ' + \
                     fmt_param('Second HVL Al', res['hvl_2_al'], 'mm Al', 
                               decimal_places=4) + ' ' + \
                     fmt_param('Homogeneity Coefficient Al', res['hc_al'], '-',
                               decimal_places=4) + \
                     fmt_param('Effective Energy Al', res['eeff_al'], 'keV',
                               decimal_places=4) + '\n'
        result_str = result_str + fmt_param('First HVL Cu', res['hvl_1_cu'], 
                               'mm Cu', decimal_places=4) + ' ' + \
                     fmt_param('Second HVL Cu', res['hvl_2_cu'], 'mm Cu', 
                               decimal_places=4) + ' ' + \
                     fmt_param('Homogeneity Coefficient Cu', res['hc_cu'], 
                               '-', decimal_places=4) + \
                     fmt_param('Effective Energy Cu', res['eeff_cu'], 'keV', 
                               decimal_places=4)
        return result_str

    def get_spectrum_str(self, std_results):
        """
        A method to format the spectrum energy-fluence values into a single 
        string

        :return str spectrum_str: A string with the spectrum details
        """
        k = std_results.k
        spk = std_results.spk
        spectrum_str = 'Spectrum units:' + '\n' + \
            '[keV]; [Photons keV cm^-2]' + '\n'

        spk_str = ''
        for idx, k_ in enumerate(k):
            spk_str = spk_str + '%.4f; %.4f\n' % (k_, spk[idx])

        spectrum_str = spectrum_str + spk_str
        return spectrum_str

    def get_current_state_str(self, mode, std_results=None):
        """
        A method to summarize the current state of the spekpy instance as a 
        single string
        
        :mode str mode: whether 'full' or 'minimal' output 
        :return str current_state_str: A string with spekpy's current state
        """

        current_state_str = 'Inputs\n' + '-' * 6 + '\n'
        current_state_str = current_state_str + \
            self.get_model_parameters_str() + '\n'
        current_state_str = current_state_str + \
            self.get_spectrum_parameters_str() + '\n'
        current_state_str = current_state_str + \
            self.get_filtration_str() + '\n'

        if mode == 'minimal':
            pass
        elif mode == 'full':
            current_state_str = current_state_str + 'Outputs\n' \
                + '-' * 7 + '\n'
            current_state_str = current_state_str + \
                self.get_current_results_str(std_results) + '\n'

        return current_state_str

    def prepare_save_state(self, file_name, file_type='state_usr', 
                           extra_param_name=None, extra_param=None):
        """
        A method to summarize the current state as an ordered dictionary, 
        ready for saving state

        :return dict state: An ordered dictionary summarizing current state
        """
        ts = get_current_time_stamp()

        if not file_name:
            file_name = 'spekpy_state ' + ts + Const.extension_state_file

        prepare_state = [('file_information', 
                          self.file_info_dct(file_name, file_type)),
                         ('external_spectrum', 
                          self.external_spectrum_dct()),
                         ('model_parameters', 
                          self.model_parameters_dct()),
                         ('spectrum_parameters', 
                          self.spectrum_parameters_dct()),
                         ('filters', 
                          self.filtration.filters),
                         ('flags', 
                          self.flags_dct())]

        if extra_param:
            prepare_state.append((extra_param_name, extra_param))

        state = ord_dct(prepare_state)

        return state

    def file_info_dct(self, file_name, file_type):
        """
        A method to summarize basic information for a file as an ordered 
        dictionary

        :param str file_name: Name of a file
        :param str file_type: The type of file
        :return dict file_information: An ordered dictionary containing basic
        information on file
        """
        file_information = ord_dct([('file_name', file_name),
                                    ('file_type', file_type),
                                    ('spekpy_version', self.spekpy_version),
                                    ('comment', self.comment)])

        return file_information

    def external_spectrum_dct(self):
        """
        A method to summarize an external spectrum as an ordered dictionary.

        :return dict external_spectrum: An ordered dictionary containing 
        information about an external spectrum.
        """
        external_spectrum = ord_dct([('external_spectrum', 
                                     self.external_spectrum.external_spectrum),
                                     ('ext_delim', 
                                     self.external_spectrum.ext_delim),
                                     ('ext_z', 
                                     self.external_spectrum.ext_z),
                                     ('ext_mas', 
                                     self.external_spectrum.ext_mas)])

        return external_spectrum

    def model_parameters_dct(self):
        """
        A method to summarize a spectrum model as an ordered dictionary

        :return dict model_parameters: An ordered dictionary containing 
            information about spectrum model
        """
        model_parameters = ord_dct([('kvp', 
                                     self.model_parameters.kvp),
                                    ('th', 
                                     self.model_parameters.th),
                                    ('dk', 
                                     self.model_parameters.dk),
                                    ('shift', 
                                     self.model_parameters.shift),
                                    ('physics', 
                                     self.model_parameters.physics),
                                    ('mu_data_source', 
                                     self.model_parameters.mu_data_source),
                                    ('targ', 
                                     self.model_parameters.targ)])

        return model_parameters

    def spectrum_parameters_dct(self):
        """
        A method to summarize spectrum paramaters as an ordered dictionary

        :return dict spectrum_parameters: An ordered dictionary containing 
            information about spectrum parameters
        """
        spectrum_parameters = ord_dct([('x', 
                                        self.spectrum_parameters.x),
                                       ('y', 
                                        self.spectrum_parameters.y),
                                       ('z', 
                                        self.spectrum_parameters.z),
                                       ('mas', 
                                        self.spectrum_parameters.mas),
                                       ('ref_kerma', 
                                        self.spectrum_parameters.ref_kerma),
                                       ('ref_flu', 
                                        self.spectrum_parameters.ref_flu),
                                       ('brem', 
                                        self.spectrum_parameters.brem),
                                       ('char', 
                                        self.spectrum_parameters.char),
                                       ('obli', 
                                        self.spectrum_parameters.obli)])
        return spectrum_parameters


    def flags_dct(self):
        """
        A method to summarize spectrum flags as an ordered dictionary

        :return dict flags: An ordered dictionary containing information about
            the spectrum flags
        """
        flags = ord_dct([('external_spectrum', 
                          self.flags.external_spectrum),
                         ('kerma_normalized_spectrum', 
                          self.flags.kerma_normalized_spectrum),
                         ('mas_normalized_spectrum', 
                          self.flags.mas_normalized_spectrum),
                         ('flu_normalized_spectrum', 
                          self.flags.flu_normalized_spectrum)])

        return flags


    def save_state_as_json(self, file_name, comment):
        """
        A method to save a state as a json file

        :param str file_name: Name of state file to create
        :param str comment: Comment to add
        """
        if not file_name:
            file_name = 'spekpy_state_' + get_current_time_stamp(mode='file') \
                + Const.extension_state_file
        
        if not file_name[-4:] == Const.extension_state_file:
            file_name = file_name + Const.extension_state_file

        
        file_name_full = full_file(Const.dir_data, Const.dir_state_usr, 
                                   file_name)
        
        state = self.prepare_save_state(file_name, comment)
        write_json_to_disk(state, file_name_full)

    def load_state_data_json(self, file_name):
        """
        A method to load state data from a json file

        :param str file_name: Name of state file to read
        :param dict data: Data on state
        """
        data = read_json_from_disk(file_name)
        return data


class ExternalSpectrumDef(object):
    """
    A class to hold information on an external spectrum
    """
    def __init__(self):
        self.external_spectrum = None
        self.ext_delim = None
        self.ext_z = None
        self.ext_mas = None

class ModelParameterDef(object):
    """
    A class to hold information on a spectrum model
    """
    def __init__(self):
        self.kvp = None
        self.th = None
        self.dk = None
        self.shift = None
        self.physics = None
        self.mu_data_source = None
        self.targ = None


class SpekpyFlagsDef(object):
    """
    A class to hold information on spectrum flags
    """
    def __init__(self):
        self.external_spectrum = False
        self.kerma_normalized_spectrum = False
        self.mas_normalized_spectrum = False
        self.flu_normalized_spectrum = False


class SpectrumParameterDef(object):
    """
    A class to hold information on spectrum parameters
    """
    def __init__(self):
        self.x = None
        self.y = None
        self.z = None
        self.mas = None
        self.brem = None
        self.char = None
        self.obli = None
        self.ref_kerma = None
        self.ref_flu = None


class FiltrationDef(object):
    """
    A class to hold information on extrinsic filtration (inherent plus added,
    but not anode self-filtration)
    """
    def __init__(self):
        self.filters = []
        self.mut_fil = []

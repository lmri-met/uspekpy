# For compatibility with Python2 #
from __future__ import print_function, division, absolute_import
##################################
import numpy as np
import spekpy.SpekConstants as Const
from scipy import interpolate
from spekpy.IO import full_file, find_file, read_json_from_disk

def get_data():
    """
    A function to import data from disk.

    :return function _rg: electron csda range
    :return function _pe: electron probability function
    :return function _ne: electron frequency function
    :return dictionary _line_data: dictionary with L and K line data
    :return function _nist_brem_data: scaled brem xsection function
    """
    electron_data = ElectronData()
    brem_data = BremData()
    line_data = LineData()
    _rg = electron_data.get_rg()
    _pe, _ne = electron_data.get_pe_ne()
    _line_data = line_data.get_line_data()
    _nist_brem_data = brem_data.get_cross_section_interpolator() 
    return _rg, _pe, _ne, _line_data, _nist_brem_data

def logarithmic_interpolator(dependent_variables, independent_variables, 
                             interpolation_values):
    """
    A helper function to do logarithmic interpolations

    :param array dependent_variable: The dependent variable of the 
        interpolation
    :param array independent_variable: The independent variable of the 
        interpolation
    :param array interpolation_values: The values from where to make 
        interpolations
    :return array interpolated_values: The resulting interpolated values
    """
    ln_dependent_variables = np.log(dependent_variables)
    ln_independent_variables = np.log(independent_variables)
    interpolater = interpolate.interp1d(ln_independent_variables, 
                                    ln_dependent_variables, bounds_error=False,assume_sorted=True)
    ln_interpolation_values = np.log(interpolation_values)
    ln_interpolated_values = interpolater(ln_interpolation_values)
    interpolated_values = np.exp(ln_interpolated_values)
    return interpolated_values


def get_atomic_weight_data():
    """
    A function to get atomic weight data that is used when defining a material 
    using chemical formulae.
    The atomic weights are taken from the atwts.dat file in XCOM v3

    :return numpy.array atomic_weight_data: Atomic weight data for elements 
        1-92
    """
    atomic_weight_data_file_name = full_file(Const.dir_data, Const.dir_tables,
                                             Const.file_atomic_weight_data)
    data = read_json_from_disk(atomic_weight_data_file_name)
    return data['atwts']


class MuData:
    """
    A class to handle attenuation coefficients
    """
    def __init__(self, mu_data_source):
        self.mu_data_source = mu_data_source
        self.mu_over_rho_data = None
        self.__load_mu_data()

    def __load_mu_data(self):
        """
        A method to load the linear attenuation coefficient data from the 
        mu_data_source

        :return: The linear attenuation coefficient data is loaded into the 
        attribute mu_data in the instance
        """
        mu_data_source = self.mu_data_source
        mu_data_file_name = full_file(Const.dir_data, Const.dir_tables, 
                                      mu_data_source)
        data = read_json_from_disk(mu_data_file_name)
        self.mu_over_rho_coefficients_energies = \
            np.array(data['photon energy'],dtype=object)
        self.mu_over_rho_coefficients = \
            np.array(data['mu_over_rho'],dtype=object)
        return

    def get_mu_over_rho(self, atomic_number_element, energy_grid):
        """
        A method to get the mass attenuation coefficients for a specific 
        element and energy grid

        :param int atomic_number_element: The atomic number of an element
        :param array energy_grid: The energy grid for the mass attenuation 
            coefficients
        :return array mu: An array of mass attenuation coefficients [cm^2 g^-1]
        """
        try:
            mass_attenuation_coefficients = np.array(
                self.mu_over_rho_coefficients[atomic_number_element - 1]
                )
            mass_attenuation_coefficient_energies = (
                np.array(self.mu_over_rho_coefficients_energies
                         [atomic_number_element - 1]
                         ) * Const.conversion_MeV2keV)
            mu_over_rho = logarithmic_interpolator(
                mass_attenuation_coefficients, 
                mass_attenuation_coefficient_energies,
                energy_grid)
            return mu_over_rho
        except:
            raise Exception(
                'Could not interpolate mass attenuation coefficients!')

    def get_mu_over_rho_composition(self, composition_name, energy_grid):
        """
        A method to calculate the compositional mass attenuation coefficient. 
        The composition is defined in a material definition file located in the
        matls_def or matl_usr directory.

        :param str composition_name: The name of the material composition file
        :param array energy_grid: The energy grid for the compositional mass 
            attenuation coefficients
        :return array, float mu_over_rho_total, rho: An array with the mass 
            attenuation coefficients
        [cm^2 g^-1] for the material as well as the density of the composition 
            [g cm^-3]
        """
        # Load the definition of the material composition. First look in the 
        # ... directory with user-defined materials, if the definition file 
        # ... does not exist there, look in the directory with default defined
        # ... materials.
        try:
            material_composition_file = find_file(composition_name, 
                                    Const.extension_matl_composition,
                                    [Const.dir_matl_usr, Const.dir_matl_def])

            composition_data = read_json_from_disk(material_composition_file)
            rho_composition = composition_data['composition']['density']
            composition = [tuple(filt) for filt 
                           in composition_data['composition']['elements']]
            number_of_elements = \
                composition_data['composition']['number_of_elements']
            
            # Pre-allocate array
            mu_over_rho_composition = np.zeros([number_of_elements,
                                                len(energy_grid)]) 

            # Loop through array with elements and weights, get mass 
            # ... attenuation coefficients and append to array.
            for element_index, element in enumerate(composition):
                atomic_number_element = element[0]
                element_weight = element[1]
                mu_over_rho_element = self.get_mu_over_rho(
                    atomic_number_element, energy_grid)
                mu_over_rho_composition[element_index, :] = \
                    element_weight * mu_over_rho_element

            mu_total = np.sum(mu_over_rho_composition, axis=0)
            return mu_total, rho_composition

        except:
            raise Exception('Error when loading material composition file!')

    def get_mu_composition(self, composition_name, energy_grid):
        """
        A method to calculate the attenuation coefficient for a material 
        composition

        :param str composition_name: The name of the material composition file
        :param array energy_grid: The energy grid for the compositional 
            attenuation coefficients [keV]
        :return array mu_composition: An array with attenuation coefficients 
            for the material
        """
        mu_over_rho_composition, rho_composition = \
            self.get_mu_over_rho_composition(composition_name, energy_grid)
        mu_composition = mu_over_rho_composition * rho_composition
        return mu_composition

    def get_mu_t(self, composition_name, energy_grid, composition_thickness):
        """
        A method to get the product of attenuation coefficients (differential 
        in energy) of a composition and the thickness of a composition

        :param str composition_name:
        :param array energy_grid: The energy grid for the compositional 
            attenuation coefficients [keV]
        :param float composition_thickness: The thickness of the composition 
            [mm]
        :return float mu_material_thickness_product: The product of attenuation
            coefficients and the thickness of the composition
        """
        mu_composition = self.get_mu_composition(composition_name, energy_grid)
        mu_material_thickness_product = mu_composition * \
            composition_thickness * Const.conversion_mm2cm
        return mu_material_thickness_product


class MuEnAirData:
    """
    A class to handle mass absorption coefficients for air
    """
    def __init__(self, muen_data_source):
        """
        Init MuEnAir

        :param str muen_data_source: The muen_air data source, either 'nist' 
            or 'pene'
        """
        self.muen_over_rho_air_data_source = muen_data_source
        self.muen_over_rho_air_coefficients_energies = None
        self.muen_over_rho_air_coefficients = None
        self.__load_muen_over_rho_air_data()

    def __load_muen_over_rho_air_data(self):
        """
        An internal method that loads the muen_air data

        :return: Populates the class attributes muen_energy_grid and 
            muen_air_values with data
        """
        muen_over_rho_air_data_source = self.muen_over_rho_air_data_source
        muen_over_rho_air_data_file_name = full_file(Const.dir_data, 
                            Const.dir_tables, muen_over_rho_air_data_source)
        data = read_json_from_disk(muen_over_rho_air_data_file_name)
        self.muen_over_rho_air_coefficients_energies = \
            np.array(data['photon energy'])
        self.muen_over_rho_air_coefficients = \
            np.array(data['muen_over_rho_air'])
        return

    def get_muen_over_rho_air(self, energy_grid):
        """
        A method to interpolate the mass absorption coefficients for air for 
        specified energies from the data tables

        :param array energy_grid: The desired energy for the mass absorption 
            coefficients
        :return array muen_over_rho: The mass absorption coefficients
        """
        try:
            mass_absorption_coefficient_air = \
                self.muen_over_rho_air_coefficients
            mass_absorption_coefficient_air_energies = \
                self.muen_over_rho_air_coefficients_energies * \
                    Const.conversion_MeV2keV
            muen_over_rho_air = logarithmic_interpolator(
                mass_absorption_coefficient_air,
                mass_absorption_coefficient_air_energies, energy_grid)

            return muen_over_rho_air
        except:
            raise Exception(
                'Could not interpolate mass absorption coefficients!')


class BremData:
    """A class to handle the NIST bremsstrahlung cross-sections"""

    def __init__(self):
        self.brem_data = None
        self.__load_brem_data()

    def __load_brem_data(self):
        """
        An internal method to load the NIST bremsstrahlung cross-sections from
        a file

        :return: Populates the class attributes electron_kinetic_energy, 
        brem_photon_energy_fraction and brem_scaled_xsection with data
        """
        brem_data_file_name = full_file(Const.dir_data, Const.dir_tables, Const.file_brem_data)
        data = read_json_from_disk(brem_data_file_name)
        self.electron_kinetic_energy = np.array(data['ebr'])
        self.brem_photon_energy_fraction = np.array(data['ubr'])
        self.brem_scaled_xsection = np.array(data['xbr'])
        return

    def get_cross_section_interpolator(self):
        """
        An internal method to make an interpolation function to get 
        bremsstrahlung cross-sections
        
        :return cross_section_interpolator: function that interpolates brems 
            xsection
        """
        # Incident electron kinetic energy [MeV]
        # Fraction of incident energy carried away by photon [1]
        # Scaled cross-section differential in photon energy [mb]
        ebr = self.electron_kinetic_energy
        ubr = self.brem_photon_energy_fraction
        xbr = self.brem_scaled_xsection
        cross_section_interpolator = interpolate.RectBivariateSpline(
            ebr * Const.conversion_MeV2keV, ubr, xbr, kx=3, ky=3)
        return cross_section_interpolator


class ElectronData:
    """
    A class to handle the import of data tables that are used for the spekpy
    model
    """

    def __init__(self):
        self.__load_electron_data()

    def __load_electron_data(self):
        """
        An internal method to load electron penetration data from a file

        :return: Populates the class attributes pe_data, ne_data and range_data
        """
        pe_data_file_name = full_file(Const.dir_data, Const.dir_tables, 
                                      Const.file_pe_data)
        self.pe_data = read_json_from_disk(pe_data_file_name)
        ne_data_file_name = full_file(Const.dir_data, Const.dir_tables, 
                                      Const.file_ne_data)
        self.ne_data = read_json_from_disk(ne_data_file_name)
        range_data_file_name = full_file(Const.dir_data, Const.dir_tables, 
                                         Const.file_range_data)
        self.range_data = read_json_from_disk(range_data_file_name)

    def get_pe_ne(self):
        """
        An internal method to get the data for the conditional probability 
        function for electrons as well as the number of electrons of a certain
        energy at a certain depth.

        :return function pe: condition probability function
        :return function ne: electron frequency
        """
        # E0 is The electron initial kinetic energy [keV]
        # u is the fraction of the initial kinetic energy at depth [1]
        # t_scaled is The depth scaled by the CSDA range [1]       
        # Nearest-neighbour interpolation is used on a dense pre-interpolated 
        # ... 3D or 2D grid
        dE0 = self.pe_data['dE0']
        dt_scaled = self.pe_data['dt_scaled']
        du = self.pe_data['du']
        pe_data_grid = np.array(self.pe_data['pe'])

        def pe(E0, t_scaled, u):
            iu = np.rint(u / du).astype(int)
            iE0 = np.rint(E0 / dE0).astype(int)
            it_scaled = np.rint(t_scaled / dt_scaled).astype(int)
            v = pe_data_grid[iE0][np.ix_(it_scaled, iu)]
            return v

        ne_data_grid = np.array(self.ne_data['ne'])
        def ne(E0, t_scaled):
            iE0 = np.rint(E0 / dE0).astype(int)
            it_scaled = np.rint(t_scaled / dt_scaled).astype(int)
            v = ne_data_grid[iE0][np.ix_(it_scaled)]
            return v

        return pe, ne

    def get_rg(self):
        """
        An internal method to get the data for electron CSDA range

        :return: function: interpolatation function for range as function of 
            energy
        """
        erng = np.array(self.range_data['E0']) * Const.conversion_MeV2keV
        rng = np.array(self.range_data['range_csda'])
        return interpolate.interp1d(erng, rng, bounds_error=False,assume_sorted=True)

class LineData:
    """
    A class to handle the import of characteristic line data used for the 
    spekpy model
    """
   
    def __init__(self):
        """
        An internal method to load line data from a file

        :return: Populates the class attribute data
        """
        line_data_file_name = full_file(Const.dir_data, Const.dir_tables, 
                                        Const.file_line_data)
        self.data = read_json_from_disk(line_data_file_name)
    
    def get_line_data(self):
        """
        An internal method to get the data for characteristic lines

        :return: dictionary: dictionary containing line data
        """
        return self.data

# This next line means that Electron, Bremsstrahlung and Line data is read in 
# ... when DataTables is imported 
data = get_data()


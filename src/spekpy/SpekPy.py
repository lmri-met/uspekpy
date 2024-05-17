# For compatibility with Python2 #
from __future__ import print_function, division, absolute_import
##################################
from . import __version__
import spekpy.SpekConstants as Const
from spekpy.IO import find_file, path_file, read_spectrum_from_disk, \
    get_script_path, is_file, print_matl_info, print_matls_in_group, \
    print_states, print_matls, delete_file
from spekpy.SpekState import State, FiltrationDef
from spekpy.SpekModel import SpekModel
import copy
from spekpy.SpekTools import load_mu_data, change_filtration, \
    calculate_air_kerma_from_spectrum, calculate_mean_energy_from_spectrum, \
    calculate_effective_energy_from_spectrum, calculate_fluence_from_spectrum,\
    calculate_required_filter_thickness, \
    calculate_first_half_value_layer_from_spectrum, \
    calculate_second_half_value_layer_from_spectrum, \
    calculate_output_arrays, \
    calculate_homogeneity_coefficient_from_spectrum, \
    StandardResults, generate_spectrum, make_composition_def
from spekpy.DataTables import data
from numpy import isclose
from scipy.ndimage import gaussian_filter


class Spek:

    def __init__(self,
            kvp=None, th=None, dk=None, mu_data_source=None, physics=None,
            x=None, y=None, z=None, mas=None, brem=None, char=None, obli=None,
            comment=None, targ=None, shift=None, init_default=True):
        """Constructor method for the Spek class 

        SEE THE WIKI PAGES FOR MORE INFO:
        https://bitbucket.org/spekpy/spekpy_release/wiki/Further_information
        
        :param float kvp: tube potential [kV] (default: depends on target)
        :param float th: anode angle [degrees] (default: 12 degrees)
        :param float dk: energy bin width [keV] (default: 0.5 keV)
        :param string: mu_data_source (default: depends on physics model)
            options: ('pene' or 'nist')
        :param float physics: physics model (default: 'casim')
            options: ('casim', 'kqp', 'spekpy-v1', 'spekcalc',
                      'diff', 'uni', 'sim', 'classical')
        :param float x: displacement from central axis in anode-cathode 
            direction [cm] (default: 0 cm)
        :param float y: displacement from central axis in orthogonal 
            direction [cm] (default: 0 cm)
        :param float z: focus-to-detector distance [cm] (default: 100 cm)
        :param float mas: the tube current-time product [mAs] (default: 1 mAs)
        :param logical brem: whether bremsstrahlung x rays requested 
            (default: true)
        :param logical char: whether characteristic x rays requested
            (default: true)
        :param logical obli: whether increased oblique paths through filtration
            are assumed for off axis positions (default: true)
        :param string comment: any text annotation the user wishes to add
        :param string targ: the anode target material (default: 'W')
            options: ('W', 'Mo', or 'Rh')
        :param float shift: optional fraction of an energy bin to shift the
            energy bins (useful when matching to measurements) (default: 0.0)
        """
        self._rg, self._pe, self._ne, self._line_data, self._nist_brem_data \
            = data
        self.state = State()
        self.state.spekpy_version = __version__
        self.state.script_path = get_script_path()
        self.mu_data = None
        self.muen_air_data = None
        self.model=SpekModel()

        if init_default:
            # Assign defaults if not specified
            th = 12.0 if th is None else th
            dk = 0.5 if dk is None else dk
            targ = 'W' if targ is None else targ
            if kvp is None:
              kvp = 100.0 if targ == 'W' else 30.0  
            x = 0.0 if x is None else x
            y = 0.0 if y is None else y
            z = 100.0 if z is None else z
            mas = 1.0 if mas is None else mas
            brem = True if brem is None else brem
            char = True if char is None else char
            obli = True if obli is None else obli
            comment = None if comment is None else comment
            physics = Spek.alias('casim') if physics is None \
                else Spek.alias(physics.lower())
            if mu_data_source is None:
                if physics == 'spekcalc' or physics == 'spekpy-v1':
                  mu_data_source = 'nist'
                else:
                  mu_data_source = 'pene'
            # Check validity of inputs
            if physics == 'spekcalc' or physics == 'spekpy-v1':
                if targ != 'W':
                    raise Exception("Only tungsten anode available for "
                        "selected physics mode (targ = 'W')")
            if targ not in ['W', 'Mo', 'Rh']:
                raise Exception("Requested anode material is not available "
                        "(targ = 'W', 'Mo' or 'Rh')")
            if mu_data_source != 'nist' and mu_data_source != 'pene':
                raise Exception("Requested mu_data_source is unrecognized "
                                "(mu_data_source = 'nist' or 'pene')")
            if physics == 'spekcalc' or physics == 'spekpy-v1':
                if kvp<10. or kvp>1000.:
                    raise Exception("Requested kVp is out of range for "
                                    "selected physics model "
                                    "(kvp = 10 to 1000 kV)")
            else:
                if targ == 'W':
                    if kvp<10. or kvp>500.:
                        raise Exception("Requested kVp is out of range for "
                                        "selected physics model and target"
                                        "(kvp = 10 to 500 kV)")
                else:
                    if kvp<20. or kvp>50.:
                        raise Exception("Requested kVp is out of range for "
                                        "selected physics model and target"
                                        "(kvp = 20 to 50 kV)")
            if shift is not None:
                if shift > 0.5 or shift < -0.5:
                    raise Exception("Requested bin shift is outside allowed "
                                    "range (shift = -0.5 to +0.5)")
            # Assign parameters to state          
            self.set_state_parameters(kvp=kvp, th=th, dk=dk, physics=physics, 
                mu_data_source=mu_data_source, x=x, y=y, z=z, mas=mas, 
                brem=brem, char=char, obli=obli, targ=targ, shift=shift)

    @staticmethod
    def alias(name):
        """
        A static method to map user physics modes strings to the strings used
        internally in the software

        :param str name: The name of the physics mode
        :return str alias_name: The name used internally for the physics mode
        """
        aliases = {'legacy':'spekcalc','default':'spekpy-v1',
            'kqp':'spekpy-v2-kqp','sim':'spekpy-v2-sim','uni':'spekpy-v2-uni',
            'casim':'spekpy-v2-casim','diff':'spekpy-v2-diff',
            'classical':'spekpy-v2-classical'}
        aliases_list = list(aliases.keys())
        for alias_name in aliases_list:
            aliases[aliases[alias_name]] = aliases[alias_name]
        alias_name = aliases[name]
        return alias_name

    def set_state_parameters(self, **kwargs):
        """
        An internal method to set arbitrary attributes in the spekpy state. 
        If an attribute that affects the model (e.g., kvp, th, dk, physics, 
        mu_data_source) is in the argument, the spectrum parameters will be 
        updated.

        :param kwargs: Keyword arguments that are supported by the spekpy state
        :return :
        """
        if 'physics' in kwargs.keys():
            kwargs['physics'] = self.alias(kwargs['physics'])
        initial_state = copy.deepcopy(self.state)
        if kwargs:
            model_parameters_changed = False
            update_model = False
            update_external = False
            # Loop through the keywords, validate if they exist in the spekpy 
            # ... state and update if model parameters have changed
            for keyword in kwargs:
                if hasattr(self.state.external_spectrum, keyword):
                    setattr(self.state.external_spectrum, keyword,
                            kwargs[keyword])
                elif hasattr(self.state.model_parameters, keyword):
                    # Only update if the model parameters really have changed.
                    if getattr(self.state.model_parameters, keyword) is \
                        not kwargs[keyword]:
                        setattr(self.state.model_parameters, keyword, 
                                kwargs[keyword])
                        model_parameters_changed = True
                elif hasattr(self.state.spectrum_parameters, keyword):
                    setattr(self.state.spectrum_parameters, keyword, 
                            kwargs[keyword])
                else:
                    raise Exception("Keyword argument " + keyword + 
                                    " not recognized!")

            if  self.state.external_spectrum.external_spectrum is None:
                if model_parameters_changed:
                    update_model = True
            else:
                if model_parameters_changed:
                    update_external = True

            # For the two physics models below, the state needs updating if
            # ... the take-off angles changes
            if self.state.model_parameters.physics=='spekpy-v2-kqp' or \
                self.state.model_parameters.physics=='spekpy-v2-sim':
                x=self.state.spectrum_parameters.x
                y=self.state.spectrum_parameters.y
                z=self.state.spectrum_parameters.z
                x_init=initial_state.spectrum_parameters.x
                y_init=initial_state.spectrum_parameters.y
                z_init=initial_state.spectrum_parameters.z
                if z_init is not None:
                    if not isclose(x_init/z_init,x/z,atol=1e-5) or \
                        not isclose(y_init/z_init,y/z,atol=1e-5):
                        update_model = True
            
            if update_model:
                # Re-initialize the model parameters
                self.mu_data, self.muen_air_data = load_mu_data(
                    self.state.model_parameters.mu_data_source)
                self.spectrum_from_model()
                current_filtration = self.state.filtration.filters
                self.state.filtration = FiltrationDef()
                self.multi_filter(current_filtration)
            
            if update_external:
               self.mu_data, self.muen_air_data = load_mu_data(
                   self.state.model_parameters.mu_data_source)
            
            # Manage normalization by reference air kerma of fluence
            if 'ref_kerma' in kwargs or 'ref_flu' in kwargs:
                if 'ref_kerma' in kwargs and 'ref_flu' in kwargs:
                    raise Exception("A reference air kerma and reference" 
                                    "fluence cannot both be specified!")
                elif 'ref_kerma' in kwargs:
                    self.state.spectrum_parameters.ref_flu = None
                    kerma = self.get_kerma(ref_kerma=None)
                    self.model.norm = \
                        self.state.spectrum_parameters.ref_kerma / kerma
                elif 'ref_flu' in kwargs:
                    self.state.spectrum_parameters.ref_kerma = None 
                    flu = self.get_flu(ref_flu=None)
                    self.model.norm = \
                        self.state.spectrum_parameters.ref_flu / flu   
            else:
                self.model.norm = None
                self.state.spectrum_parameters.ref_kerma = None
                self.state.spectrum_parameters.ref_flu = None
        return self

    def parameters_for_calculation(self, **kwargs):
        """
        A function to handle parameters in spekpy calculations without setting
        the state of spekpy (transient) 

        This function recognized the following keyword arguments: 
        x, y, z, brem, char, obli

        :param kwargs: A set of keyword arguments for parameters that can 
            temporarly be used for calculations
        """
        calc_params = copy.deepcopy(self.state.spectrum_parameters)


        for keyword in kwargs:
            if hasattr(calc_params, keyword):
                setattr(calc_params, keyword, kwargs[keyword])
            else:
                raise Exception('Keyword argument '+ keyword +' not recognized')
 
    
        if getattr(calc_params, 'ref_kerma') is not None and \
            getattr(calc_params, 'ref_flu') is not None:
            raise Exception("A reference air Kerma and reference fluence "
                            "cannot both be specified!")
       
        if 'x' in kwargs or 'y' in kwargs or 'z' in kwargs:
            if self.state.model_parameters.physics == Spek.alias('kqp') or \
                self.state.model_parameters.physics == Spek.alias('sim'):
                x=self.state.spectrum_parameters.x
                y=self.state.spectrum_parameters.y
                z=self.state.spectrum_parameters.z
                x_calc=calc_params.x
                y_calc=calc_params.y
                z_calc=calc_params.z
                if not isclose(x_calc/z_calc,x/z,atol=1e-5) or \
                    not isclose(y_calc/z_calc,y/z,atol=1e-5):
                    raise Exception('Cannot update spatial position.\n' +
                    'Specified value of "physics" keyword prohibits a ' +
                    'change in take-off angles on-the-fly.\n' +
                    'Try using the "set()" method to change x, y or z.')
        
        return calc_params
    
    def spectrum_from_model(self):
        """
        Internal method to get spectra from spekpy's model of photon emission

        :return:
        """
        self.model = SpekModel().get_spectrum_parameters(self)

    def spectrum_from_external_source(self,z,mas):
        """
        Internal method used in loading an eternal spectrum from a file

        :param float z: The focus-to-detector distance for the file spectrum
        :param float mAs: The exposure setting for the file spectrum
        :return:
        """
        k, brem, char, dk = read_spectrum_from_disk(self)
        number_of_incident_electrons = mas * Const.conversion_A2mA \
            / Const.electron_charge
        # Convert external fluence spectrum to photons per solid angle per 
        # ... incident electron
        self.model.brem_k = brem * z**2 / number_of_incident_electrons
        self.model.char_k = char * z**2 / number_of_incident_electrons
        self.model.k = k
        self.state.model_parameters.dk = dk
        self.state.model_parameters.kvp = max(k) + 0.5*dk
        
    def set(self, **kwargs):
        """
        A method to set parameters in the spekpy state

        :param kwargs: A variable number of keyword arguments to change the 
            state of spekpy
        :return:
        """
        self.set_state_parameters(**kwargs)
        return self

    def summarize(self, mode='minimal'):
        """
        A method to print a summary of the current state and results thereof to
        the console

        :param str mode: The mode of summarization. Either 'full' or 'minimal'
            (default)
        :return:
        """
        if mode == 'full':
            std_results = self.get_std_results()
            summarization = self.state.get_current_state_str(mode, std_results)
        elif mode == 'minimal':
            summarization = self.state.get_current_state_str(mode)
        else:
            raise Exception("The mode must be either full or minimal "
                            "(default)")

        print(summarization)

    def filter(self, matl, t):
        """
        A method to alter the spectral filtration

        Example usage:
        spk.filter('Al', 1)

        :param str matl: The name of the desired filter material
        :param float t: The thickness of the desired filtration
        :return:
        """
        change_filtration(self, matl, t)
        return self

    def multi_filter(self, filter_list):
        """
        A method to apply multiple added filtrations at the same time

        :param filter_list:
        :return:
        """
        try:
            for filter_tuple in filter_list:
                change_filtration(self, filter_tuple[0], filter_tuple[1])
        except:
            raise Exception("Could not use multi-filter method. Check syntax "
                            "of your filter list!")
        return self

    def get_k(self, **kwargs):
        """
        A method to get the array of energies for the current spectrum

        :param kwargs: Keyword arguments to change parameters that are used 
            for the calculation
        :return array k: Array with photon energies (mid-bin values) [keV]
        """
        k = self.model.k
        return k

    def get_spk(self, **kwargs):
        """
        A method to get the spectrum for the parameters in the current spekpy 
        state

        :param kwargs: Keyword arguments to change parameters that are used for
            the calculation
        :return array spk: Array with photon distribution 
            [Photons cm^-2 keV^-1]
        """
        calc_params = self.parameters_for_calculation(**kwargs)
        spk = generate_spectrum(calc_params, self.state.model_parameters, 
                                self.state.external_spectrum, self.model, 
                                self.state.filtration)
        return spk

    def get_spectrum(self,edges=False, flu=True, diff=True, sig=None, 
                     addend=False,**kwargs):
        """
        A method to get the energy and spectrum for the parameters in the 
        current spekpy state

        :param bool edges: Keyword argument to determine whether midbin or edge
            of bins data are returned
        :param bool addend: Keyword argument to determine whether a zero
            end point is added to the spectrum
        :param bool flu: Whether to return fluence or energy-fluence
        :param bool diff: Whether to return spectrum differential in energy
        :param kwargs: Keyword arguments to change parameters that are used for
            the calculation
        :return array k: Array with photon energies (mid-bin or edge values)
            [keV]
        :return array spk: Array with corresponding photon fluences 
            [Photons cm^-2 keV^-1], [Photons cm^-2] or [Photons cm^-2 keV^1]
            depending of values of flu and diff inputs
        """
        k = self.get_k(**kwargs)
        dk = self.state.model_parameters.dk
        spk = self.get_spk(**kwargs)
        
        # Apply Gaussian filter to degrade energy resolution if requested
        if sig is not None:
            spk = gaussian_filter(spk, sig/dk, mode='constant')
        # Provide requested definition of spectrum    
        if flu and diff:
            spk = spk
        elif not flu and diff:
            spk = spk * k
        elif flu and not diff:
            spk = spk * dk
        else:
            spk = spk * k * dk
        # Format output as requested
        k_out, spk_out = calculate_output_arrays(k, spk, edges=edges, 
                                                 addend=addend)
        return k_out, spk_out
    
    def get_kerma(self, norm=False, **kwargs):
        
        """
        A method to get the air Kerma for the current spekpy state

        :param kwargs: Keyword arguments to change parameters that are used for
            the calculation
        :param bool norm: A boolean indicating if the air Kerma should be 
            normalized by the mAs
        :return float air_kerma: The air Kerma for the current spectrum at 
            location x,y,z from the focal spot
        """
        calc_params = self.parameters_for_calculation(**kwargs)
        air_kerma = \
            calculate_air_kerma_from_spectrum(
                self, calc_params, mas_normalized_air_kerma=norm)
        return air_kerma

    def get_hvl1(self, matl='Al', **kwargs):
        """
        Method to get the first half value layer for a desired material for the
        parameters in the current spekpy state

        :param str matl: The desired material name
        :param kwargs: Keyword arguments to change parameters that are used for
            the calculation
        :return float first_half_value_layer: The calculated first half value
            layer for the desired material [mm]
        """

        calc_params = self.parameters_for_calculation(**kwargs)
        first_half_value_layer = \
            calculate_first_half_value_layer_from_spectrum(
                self, calc_params, matl)
        return first_half_value_layer
    
    def get_hvl(self, matl='Al', **kwargs):
        """
        Method to get the first half value layer for a desired material for the
        parameters in the current spekpy state

        This is identical to "get_hvl1"

        :param str matl: The desired material name
        :param kwargs: Keyword arguments to change parameters that are used for
            the calculation
        :return float first_half_value_layer: The calculated first half value 
            layer for the desired material [mm]
        """

        calc_params = self.parameters_for_calculation(**kwargs)
        first_half_value_layer = \
            calculate_first_half_value_layer_from_spectrum(
                self, calc_params, matl)
        return first_half_value_layer

    def get_hvl2(self, matl='Al', **kwargs):
        """
        A method to get the second half value layer for a desired material for
        the parameters in the current spekpy state

        :param str matl: The desired material name
        :param kwargs: Keyword arguments to change parameters that are used for
            the calculation
        :return float second_half_value_layer: The calculated second half value
            layer for the desired material [mm]
        """

        calc_params = self.parameters_for_calculation(**kwargs)
        second_half_value_layer = \
            calculate_second_half_value_layer_from_spectrum(
                self, calc_params, matl)
        return second_half_value_layer

    def get_hc(self, matl='Al', **kwargs):
        """
        A method to get the homogeneity coefficient for a desired material for
        the parameters in the current spekpy state

        :param str matl: The desired material name
        :param kwargs: Keyword arguments to change parameters that are used for
            the calculation
        :return float homogeneity_coefficient: The calculated homogeneity 
            coefficient for the desired material
        """
        calc_params = self.parameters_for_calculation(**kwargs)
        homogeneity_coefficient = \
            calculate_homogeneity_coefficient_from_spectrum(
                self, calc_params, matl)
        return homogeneity_coefficient

    def get_matl(self, matl='Al', hvl_matl='Al',hvl=False, frac=False, 
                 **kwargs):
        """
        A method to calculate the thickness of a desired material that is 
        needed to obtain a first half value layer or fraction of air Kerma 
        (transmission). Take note, the default material is Al.

        :param str matl: The desired material name
        :param float hvl: The desired first half value layer of the desired 
            material
        :param float frac: The desired air Kerma fraction
        :param kwargs: Keyword arguments to change parameters that are used for
            the calculation
        :return float required_filter_thickness: The required thickness of the
            desired material [mm]
        """
        calc_params = self.parameters_for_calculation(**kwargs)
        required_filter_thickness = \
            calculate_required_filter_thickness(self, calc_params, matl, 
                                                hvl_matl, hvl, frac)
        return required_filter_thickness

    def get_emean(self, **kwargs):
        """
        A method to get a calculation of the mean energy of the spectrum for 
        the parameters in the current spekpy state

        :param kwargs: Keyword arguments to change parameters that are used for
            the calculation
        :return float mean_energy: The mean energy of the spectrum [keV]
        """
        calc_params = self.parameters_for_calculation(**kwargs)
        mean_energy = calculate_mean_energy_from_spectrum(self, calc_params)
        return mean_energy

    def get_eeff(self, matl='Al', **kwargs):
        """
        A method to get a calculation of the effective energy of the spectrum 
        for the parameters in the current spekpy state, for a desired material

        :param str matl: The desired filter material
        :param kwargs: Keyword arguments to change parameters that are used for
            the calculation
        :return float effective_energy: The calculated effective energy [keV]
        """
        calc_params = self.parameters_for_calculation(**kwargs)
        effective_energy = \
            calculate_effective_energy_from_spectrum(self, calc_params, matl)
        return effective_energy

    def get_flu(self, **kwargs):
        """
        A method to calculate the fluence of the spectrum for the parameters in
        the current spekpy state

        :param kwargs: Keyword arguments to change parameters that are used for
            the calculation
        :return float fluence: The fluence of the spectrum [Photons cm^-2]
        """
        calc_params = self.parameters_for_calculation(**kwargs)
        fluence = calculate_fluence_from_spectrum(self, calc_params)
        return fluence

    def get_eflu(self, **kwargs):
        """
        A method to calculate the energy fluence of the spectrum for the 
        parameters in the current spekpy state

        :param kwargs: Keyword arguments to change parameters that are used for
            the calculation
        :return float efluence: The fluence of the spectrum [Photons cm^-2 keV]
        """
        calc_params = self.parameters_for_calculation(**kwargs)
        efluence = calculate_fluence_from_spectrum(self, calc_params) * \
            self.get_emean(**kwargs)
        return efluence

    def get_norm_flu(self, kerma=1, **kwargs):
        """
        A method to calculate the kerma normalized fluence

        :param float kerma: Kerma used to normalize fluence with [uGy]
        :param kwargs: Keyword arguments to change parameters that are used for
            the calculation
        :return float norm_flu: The resulting Kerma normalized fluence 
            [Photons cm^-2 uGy^-1]
        
        """
        kerma_0 = self.get_kerma(**kwargs)
        flu_0 = self.get_flu(**kwargs)
        norm_flu = flu_0 * kerma / kerma_0
        return norm_flu
          

    def get_std_results(self, **kwargs):
        """
        A method to calculate standard results for a spectrum that is defined 
        with the current spekpy state

        Standard results include:
            k: Photon energies (mid-bin) [keV]
            spk: Photon energy spectrum [Photons cm^-2 keV^-1]
            flu: The fluence of the spectrum [Photons cm^-2]
            kerma: The air Kerma of the spectrum [uGy]
            emean: The mean energy of the spectrum [keV]
            hvl_1_al: First half value layer of aluminum [mmAl]
            hvl_2_al: Second half value layer of aluminum [mmAl]
            hc_al: Homogeneity coefficient of aluminum
            eeff_al: Effective energy of spectrum relative aluminum [keV]
            hvl_1_cu: First half value layer of copper [mmCu]
            hvl_2_cu: Second half value layer of copper [mmCu]
            hc_cu: Homogeneity coefficient of copper
            eeff_cu: Effective energy of spectrum relative copper [keV]

        :param kwargs: Keyword arguments to change the spekpy state
        :return StandardResults standard_results: A class containing the 
            standard results
        """
        calc_params = self.parameters_for_calculation(**kwargs)
        standard_results = \
            StandardResults().calculate_standard_results(self, calc_params)
        return standard_results

    def comment(self, comment=None):
        """
        An internal method that can be used to comment on the spekpy state. 
        This comment will be used for any exports
        or saved states.

        :param str comment: A comment for the current state of spekpy
        :return:
        """
        self.state.comment = comment
        return self

    def save_state(self, file_name=None, comment=None):
        """
        A method to save the state of spekpy to disk
        :param str comment: A comment that will be saved with the state
        :param str file_name: The name of the file
        :return:
        """

        if comment:
            self.state.comment = comment

        self.state.save_state_as_json(file_name, comment)

    @classmethod 
    def load_state(cls, state_name):
        """
        A class method to load a saved spectrum state
        
        :param str state_name: Name of the saved spectrum state
        :return Spek spek: Instance of spekpy Spek class
        """ 
        spek = Spek(init_default=False)
        state_file = \
            find_file(file_name=state_name, extension=Const.extension_state_file,
                      directories=[Const.dir_state_def, Const.dir_state_usr])
        
        state_data = spek.state.load_state_data_json(state_file)
        spek.state.spectrum_parameters.__dict__.update(
            **state_data['spectrum_parameters'])
        spek.state.model_parameters.__dict__.update(
            **state_data['model_parameters'])
        spek.state.external_spectrum.__dict__.update(
            **state_data['external_spectrum'])

        if spek.mu_data is None and spek.muen_air_data is None:
            spek.mu_data, spek.muen_air_data = \
                load_mu_data(spek.state.model_parameters.mu_data_source)
        
        spek.state.model_parameters.physics=spek.alias(
            spek.state.model_parameters.physics)
        spek.spectrum_from_model()
        filtration = [tuple(filt) for filt in state_data['filters']]
        spek.multi_filter(filtration)
        return spek

    @classmethod
    def load_from_file(cls, spectrum_name, spectrum_delimeter, z=100.0, 
                       mas=1.0, mu_data_source='nist'):
        """
        A class method to import an external spectrum from a file. The external
        spectrum is assumed to have been measured/estimated along the central 
        ray of the x-ray tube,
        i.e., that the x- and y- coordinates are 0 repsectively.

        :param str spectrum_name: The file name of the external spectrum
        :param str spectrum_delimeter: The delimeter used to separate energy 
            and fluence
        :param float ext_z: The Z-coordinate where the external spectrum was
            measured/estimated [cm]
        :param float ext_mas: The tube load that the external spectrum has been
            measured/estimated for
        :param str mu_data_source: The data source for attenuation coefficient,
            i.e., 'nist' or 'pene'
        :return Spek spek: An instance of spekpy where the external spectrum 
            has been loaded
        """       
        spek = Spek(init_default=False)
        external_spectrum = is_file(spectrum_name, spek.state.script_path)
        # Set parameters in self.state.external_spectrum
        spek.set_state_parameters(external_spectrum=external_spectrum, 
                        ext_delim=spectrum_delimeter, ext_z=z, ext_mas=mas)
        # Set parameters in self.state.spectrum_parameters
        spek.set_state_parameters(x=0, y=0, z=z, mas=mas, 
                                  obli=True, brem=True, char=True)
        spek.state.model_parameters.mu_data_source = mu_data_source
        spek.mu_data, spek.muen_air_data = load_mu_data(
            spek.state.model_parameters.mu_data_source)
        spek.spectrum_from_external_source(z,mas)
        return spek


    @staticmethod
    def make_matl(matl_name=None, matl_density=None, wt_matl_comp=None, 
                  chemical_formula=None, matl_comment=None):
        """
        A static method to create a new material definition

        :param str matl_name: The name of the material
        :param float matl_density: The density of the material
        :param list wt_matl_comp: A list of tuples containing information about
            the atomic number and mass fraction by weight of the material, 
            e.g., [(1, 0.111898), (8, 0.888102)] for water
        :param str chemical_formula: The chemical formula of the material, 
            e.g., 'H2O' for water
        :param str matl_comment: A comment to be included in the material 
            definition file
        """
        make_composition_def(matl_name, matl_density, wt_matl_comp, chemical_formula, matl_comment)


    @staticmethod
    def remove_matl(matl_name):
        """
        A static method to delete a material in the usr directory

        :param str matl_name: The name of a material in usr directory
        """
        if matl_name.endswith('.comp'):
            file_name_parts=[Const.dir_data, Const.dir_matl_usr, matl_name]
        else:
            file_name_parts=[Const.dir_data, Const.dir_matl_usr, 
                             matl_name + '.comp']
        delete_file(file_name_parts)

    @staticmethod
    def remove_state(state_name):
        """
        A static method to delete a state in the usr directory

        :param str state_name: The name of a state in usr directory
        """
        if state_name.endswith('.state'):
            file_name_parts=[Const.dir_data, Const.dir_state_usr, state_name]
        else:
            file_name_parts=[Const.dir_data, Const.dir_state_usr, 
                             state_name + '.state']
        delete_file(file_name_parts)
    
    @staticmethod
    def show_matls(matl_name=None, matl_group=None, matl_dir=None):
        """
        A static method to print information about a given material definition,
        alternatively, to print materials that belong to a logical group or 
        directory ("usr" ofr "def"). The following logical groupings include:
        materials defined by the ICRU or materials defined by the ICRP.

        :param str matl_name: The name of the given material
        :param str matl_group: The name of a material group
        :param str matl_dir: The name of a directory ("usr" or "def")
        """
        if matl_name and matl_group is None and matl_dir is None:
            print_matl_info(matl_name)
        elif matl_group and matl_name is None and matl_dir is None:
            print_matls_in_group(matl_group)
        elif matl_dir and matl_name is None and matl_group is None:
            print_matls(matl_dir)
        elif matl_dir is None and matl_name is None and matl_group is None:
            print_matls()
        else:
            raise Exception("Specify only one of matl_name, matl_group, and "
                            "matl_dir.")

    @staticmethod
    def show_states(state_dir=None):
        """
        A static method to print saved states in a specified directory ("usr" 
        or "def") or print all saved states

        :param str state_dir: The name of the given state directory ("usr" 
            or "def")
        """
        print_states(state_dir)
        

    def export_spectrum(self, file_name=None, comment=None, delim=';'):
        """
        An internal method that can be used to export a spectrum to disk

        :param str file_name: The file name of the exported spectrum
        :param str comment: A comment that will be included in the exported 
            spectrum
        :return:
        """
        path, file_name_nopath, slash = path_file(file_name)
        
        if path == '':
            file_name = get_script_path() + slash + file_name_nopath

        if comment:
            self.state.comment = comment
        if comment is None: 
            self.state.comment = ''

        std_results = self.get_std_results()
        spk_char = self.get_spk(brem=False)

        self.state.export_spectrum_to_disk(file_name, delim, 
                                           std_results, spk_char)


    @staticmethod
    def clone(spekpy_obj):
        """
        A static method to clone a Spek class instance

        :param Spek spekpy_obj: The Spek class instance
        :return Spek spekpy_clone: The cloned class instance
        """
        spekpy_clone = copy.deepcopy(spekpy_obj)
        return spekpy_clone



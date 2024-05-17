# For compatibility with Python2 #
from __future__ import print_function, division, absolute_import
##################################
import os
from numpy import sin, cos, arctan, arccos, deg2rad, degrees
from numpy import exp, sqrt, log, sum, abs, floor, round, divide, trapz
from numpy import pi, nan
from numpy import array, zeros, ones, insert, squeeze, meshgrid, arange
from numpy import load, finfo, argmax, where, linspace, matmul
import scipy.interpolate as interp
from scipy.ndimage import map_coordinates

#[1] Omar A, Andreo P and Poludniowski G. A model for the emission of K 
# ... and L x rays from an x-ray tube. NIM B 2018;437:36-47.
#[2] Omar A, Andreo P and Poludniowski G. A model for the energy and angular 
# ... distribution of x rays emitted from an x-ray tube. Part I. 
# ... Bremsstrahlung production. Accepted by Med Phys 2020.
#[3] Omar A, Andreo P and Poludniowski G. A model for the energy and angular 
# ... distribution of x rays emitted from an x-ray tube. Part II. 
# ... Validation of x-ray spectra from 20 to 300 kV. Accepted by Med Phys 2020.

class SpekAniso:

    
    def __init__(self,target,shape,diffuse,brsrc,x,y,z,k,th,kvp):
        """Constructor method for the SpekAniso class 

        :param string target: The abbreviation symbols for the target element
        :param string shape: The brems shape-function 
            options: ('kqp', 'sim', 'uni' or 'casim')
        :param bool diffuse: Whether instant diffusion is assumed for electrons
        :param string brsrc: The source of first differential brem xsection
            options: ('nist' or 'classical')
        :param float x: The lateral position in anode-cathode direction [cm]
        :param float y: The lateral position in other direction [cm]
        :param float z: The focus-to-detection distance [cm]
        :param float array k: The x-ray energy bins to evaluate at [keV]
        :param float th: The anode angle [degrees]
        :param float kvp: The tube potential [kV]
        """
        # The location of this file
        filedir = os.path.dirname(os.path.abspath(__file__))
        # Path to data file directory
        self.dirpath = os.path.join(filedir, 'data', 'tables', 'advanced')
        # Dictionary with targets modelled and corresponding atomic number
        target_Z = {'W':74, 'Mo':42, 'Rh':45}
        # Dictionary with targets modelled at corresponding atomic mass
        target_amu = {'W':183.84, 'Mo':95.94, 'Rh':102.906}
        # The atomic weight of the target
        self.amu = target_amu[target]
        # Electron mass [keV/c**2]
        self.me = 510.9989 
        # Classical electron radius [cm]
        self.re = 2.81794e-13
        #  Reduced Planck constant [kev s]
        self.hbar = 6.626e-34/(2*pi*1.602e-16)
        # Speed of light [cm/s]
        self.c = 2.998e10
        # X-ray energy bins
        self.k = k
        # Incident electron energy
        self.E0 = kvp
        # The model for bremsstrahlung shape function
        self.shape = shape
        # Whether instant diffusion is assumed
        self.diffuse = diffuse
        # The first-differential bremsstrahlung cross-section to use
        self.brsrc = brsrc
        # The atomic nu ber of the target
        self.Z = target_Z[target]
        # Angles for x and y axes w.r.t z
        bet = arctan(y / z)
        alp = arctan(x / z)
        # Change to takeoff angles
        self.varphi = deg2rad(th) + alp
        self.vartheta = bet
        # Get bremsstrahlung and characterstic contributions
        self.brem_k, self.brem_kt, \
            self.anode_self_filtration , self.t = self.__Brems()
        self.char_k, self.char_kt, \
            self.anode_self_filtration_char, self.t_char = self.__Char()
  
    
    def __get_range_fn(self):
        """Internal method to generate a function for electron CSDA range as a 
        function of energy. 
    
        :return function range_fn: Interpolation function for CSDA [g/cm**2] 
            as a fn of energy [keV] 
        """
        # Read in data file
        filepath = os.path.join(self.dirpath, 'CSDA.npz')
        data = load(filepath)
        # Extract data: electron kinetic energy grid [keV]
        E0_data = data['E0' + str(self.Z)].astype('float64')
        # Extract data: CSDA range grid [g/cm**2]
        rng_data = data['rng' + str(self.Z)].astype('float64')
        # Generate 1d interpolation function
        range_fn = interp.interp1d(E0_data, rng_data, kind='quadratic', 
                                         bounds_error = False)
        return range_fn


    def __get_logmac_fn(self, Rayleigh=True):
        """Internal method to generate log mass attenuation coefficient as a 
        function of log energy. 
    
        :param logical Rayleigh: If True then Rayleigh (coherent) contribution 
            included
        :return function logmac: Interpolation function for log(mac) as a fn of 
            log(photon energy [keV]) 
        :return float rho: Density of the element [g/cm3]
        """
        # Outputs dependent on whether Rayleigh/coherent scattering included
        if Rayleigh: # Include Rayleigh/coherent scattering
            filepath = os.path.join(self.dirpath, 
                                  'MassAttCoeff_Z' + str(self.Z) + '.npz')
        else: # Exclude Rayleigh/coherent scattering
            filepath = os.path.join(self.dirpath, 
                                  'MassAttCoeff_Z' + str(self.Z) + '_woRa.npz')
        # Read in data from file
        data = load(filepath)
        # Extract date grid for log energy of photon [keV]
        logk_data = data['logk'].astype('float64')
        # Extract data grid for log of mass attenuation coefficient [cm**2/g]
        logmac_data = data['logmac'].astype('float64')
        # Generate 1d interpolation function
        # assume_assorted set to True to prevent interp1d reordering arrays,
        # ... if precision loss creates ambiguities in the energy values at
        # ... K/L edges
        logmac_fn = interp.interp1d(logk_data, logmac_data, kind = 'linear', 
                            fill_value = 'extrapolate',assume_sorted=True)
        # Extract density [g/cm**3]
        rho = data['rho'].astype('float64')
        return logmac_fn, rho

    
    def __gen_epen_arrays(self, csda_range, E, theta_e, domega_e):
        """Internal method to interpolate electron penetration data arrays from 
        tables. 
    
        :param float csda_range: CSDA range of electrons [cm]
        :param array E: 1d array of electron energies at depth [keV]
        :param array theta_e: Array of electron orientations with respect to 
            normal incidence [degree]
        :param float domega_e: Element of solid angle [degree**2]
        :return array x_cm: 1d array of depths for electron 
            penetration arrays [cm]
        :return array NE: 2D array of electron frequency density with energy 
            E at depth x [keV**-1]
        :return array pomega_e: 3D array of electron prob density with angle 
            theta_e for electron energy E and depth x [degree**-1]
        """
        # List of targets with available electron penetration data 
        Zlist = array([42, 74])
        # The index of nearest element in list to actual target material
        iZ = abs( Zlist-self.Z ).argmin()
        # The Z of the nearest element
        Zdata=int( Zlist[iZ] )
        # Read and extract info for the electron incident energies available
        filepath = os.path.join( self.dirpath,'Z' + str(Zdata) + '_E0sim.npz' )
        data = load( filepath )
        E0_sim = data['E0_sim'].astype('float64')
        # Find the appropriate data file (fileNr depends on energy)
        fileNr = sum( insert( E0_sim[1:], 0, 0.0 ) / float(self.E0)<=1)
        # Ensure file number does not go out of range (primarily to avoid
        # ... floating point rounding errors at max and min kV values)
        fileNr = max(1,min(fileNr, E0_sim.size-1))
        # Read and extract the data file
        filepath = os.path.join( self.dirpath,'Z' + str(Zdata) + '_' + \
                                str(fileNr) + '.npz' )
        data = load( filepath )
        # 1D depth array for the data tables (as fraction of CSDA range)
        x_data = data['x'].astype('float64')
        # 1D array of electron energies for the data tables (as fraction of E0)
        u_data = data['u'].astype('float64')
        # 1D array of incident electron energies for the data tables [keV]
        E0_data = data['E0_table'].astype('float64')
        # 3D array electron frequency density for electron energy, depth and
        # ... incident electron energy [keV**-1]
        # Note that NE is differential in dE (or du*E0), not du
        NE_data = data['NE'].astype('float64')
        # Define meshfrids for u values, x indices and E0 values for interp.
        xind_data = array( range(x_data.size) ).astype('float64')
        [u_grid, xind_grid, E0_grid] = meshgrid( E / self.E0, xind_data, 
                                                   self.E0, indexing='ij')
        # Calculate the indices (floats) corresponding to u_grid values
        # Round to prevent floating point precision problems at boundaries
        imap = ( u_grid - u_data[0] ) / ( u_data[1] - u_data[0] )
        imap = imap.round(decimals=3)
        # Calculate the indices (floats) corresponding to xind_grid values
        # Note that we work with xind_grid instead of x_grid as x_data is not
        # ... uniformly spaced
        # Round to prevent floating point precision problems at boundaries
        jmap = ( xind_grid - xind_data[0] ) / ( xind_data[1] - xind_data[0] )
        jmap = jmap.round(decimals=3)
        # Calculate the indices (floats) corresponding to E0_grid values
        # Round to prevent floating point precision problems at boundaries
        kmap = ( E0_grid - E0_data[0] ) / ( E0_data[1] - E0_data[0] ) 
        kmap = kmap.round(decimals=3)
        # Outputs dependent on whether "instant diffusion" is assumed
        if self.diffuse: # Instant diffusion assumed
            theta_e_data = data['theta_e'].astype('float64')
            pomega_e_data = data['pomega_e'].astype('float64')
            # Calculate the diversion for each electron energy and depth
            # numpy where and divide used to avoid divide by zero warning
            diversion_num = sum( pomega_e_data * 
                    abs( sin( deg2rad( theta_e_data[:,None,None,None] ) ) ),
                    axis=0)
            diversion_den = sum( pomega_e_data * 
                    abs( cos( deg2rad( theta_e_data[:,None,None,None] ) ) ) * 
                    abs( sin( deg2rad( theta_e_data[:,None,None,None] ) ) ),
                    axis=0)
            diversion = where( diversion_den > 0, 
                      divide(diversion_num, diversion_den, 
                                where = diversion_den > 0),
                                2)
            # Correct NE from explicit diversion factor to diversion = 2 
            NE_data = NE_data * 2. / diversion
            # Determine expected normalization (NE integrated over
            # ... dE and dx) at self.E0 using linear interpolation
            NE_integrated = trapz(trapz(NE_data,u_data,axis=0),
            	x_data,axis=0)*E0_data
            weight = 1.0 - (self.E0 - E0_data[0]) / (E0_data[1] - E0_data[0])
            NE_integrated_expected = weight * NE_integrated[0] + \
            	(1.0 - weight) * NE_integrated[1]
            # 3D linear interpolation of NE onto the pre-defined grid
            NE = map_coordinates( NE_data, array( [imap, jmap, kmap] ), 
                                 order=1)
            # The 3rd dim is interpolated at a single value so we can remove
            # ... the singlet dim
            NE = squeeze(NE)
            # Correct actual to expected normalization
            NE_integrated_actual = (trapz( trapz(NE, E/self.E0, axis=0), 
            	x_data, axis=0)*self.E0)
            NE = NE * NE_integrated_expected / NE_integrated_actual
            # Define uniform electron direction distribution
            pomega_e = ones( [theta_e.size, E.size, x_data.size] )
            pomega_e = pomega_e / (sum( pomega_e, axis=0) * domega_e)
        else: # Explicit treatment
            # Determine expected normalization (NE integrated over
            # ... dE and dx) at self.E0 using linear interpolation
            NE_integrated = trapz(trapz(NE_data,u_data,axis=0),
            	x_data,axis=0)*E0_data
            weight = 1.0 - (self.E0 - E0_data[0]) / (E0_data[1] - E0_data[0])
            NE_integrated_expected = weight * NE_integrated[0] + \
            	(1.0 - weight) * NE_integrated[1]
            # 3D linear interpolation of NE onto the pre-defined grid
            NE = map_coordinates( NE_data, array( [imap, jmap, kmap] ), 
                                 order=1)
            # The 3rd dim is interpolated at a single value so we can remove
            # ... the singlet dim
            NE = squeeze(NE)
            # Correct actual to expected normalization
            NE_integrated_actual = (trapz( trapz(NE, E/self.E0, axis=0), 
            	x_data, axis=0)*self.E0)
            NE = NE * NE_integrated_expected / NE_integrated_actual
            # 1D array of electron directions for data tables
            theta_e_data = data['theta_e'].astype('float64')
            # 4D array of electron angular distribution for data tables
            # The four dims refer to theta_e, u, x and E0
            pomega_e_data = data['pomega_e'].astype('float64')
            # Define meshrids for interpolation
            [theta_e_grid, u_grid, xind_grid, E0_grid] = meshgrid( 
                theta_e, E / self.E0, xind_data, self.E0, indexing='ij')
            # Calculate indices (floats) corresponding to theta_e_grid values
            # Round to prevent floating point precision problems at boundaries
            imap = ( theta_e_grid - theta_e_data[0]) \
                / ( theta_e_data[1] - theta_e_data[0] )
            imap=imap.round(decimals=3)
            # Calculate indices (floats) corresponding to u_grid values
            # Round to prevent floating point precision problems at boundaries
            jmap = ( u_grid - u_data[0] ) \
                / ( u_data[1] - u_data[0] )
            jmap=jmap.round(decimals=3)
            # Calculate indices (floats) corresponding to xind_grid values
            # Round to prevent floating point precision problems at boundaries
            kmap = ( xind_grid - xind_data[0] ) \
                / ( xind_data[1] - xind_data[0] )
            kmap=kmap.round(decimals=3)
            # Calculate indices (floats) corresponding to E0_grid values
            # Round to prevent floating point precision problems at boundaries
            lmap = ( E0_grid - E0_data[0] ) / ( E0_data[1] - E0_data[0] )
            lmap=lmap.round(decimals=3)
            # 4D linear interpolation of NE onto the pre-defined grid
            # cval set to tiny nonzero so later normalization in pomega_e 
            # ... doesn't create NaN elements in pomega_e 
            pomega_e = map_coordinates( pomega_e_data, 
                            array([imap,jmap,kmap,lmap]), 
                            order=1, cval = finfo(dtype='float64').resolution)
            # The 4th dim is interpolated at a single value so we can remove 
            # ... the singlet dimension
            pomega_e = squeeze(pomega_e)
            # Ensure correct normalization
            pomega_e = pomega_e / ( sum( pomega_e, axis=0) * domega_e)
        # Return x array in units of cm rather than as a fraction of CSDA range
        x_cm = x_data*csda_range  
        return x_cm, NE, pomega_e
  
    
    def __gen_sig_arrays(self,E_grid,k_grid):
        """Internal method to interpolate bremsstrahlung cross-section 
        (dsig/dk) data arrays from tables. 
    
        :param array E_grid: 2D Electron kinetic energy mesh grids [keV]
        :param array k_grid: 2D Bremsstrahlung energy mesh grids [keV] 
        :return array sig: 2d array of bremsstrahlung cross-sections
            corresponding to E_grid/k_grid [cm**-2.keV**-1]
        :return array sig_tip: 1D array of bremsstrahlung cross-sections at tip
            corresponding to E_grid = k_grid [cm**-2.keV**-1]
        """ 
        if self.brsrc == 'nist': # NIST cross.section
            # Read the NIST bremsstrahlung cross-sections dsig/dx
            filepath = os.path.join(self.dirpath, 'NIST.npz')
            data = load(filepath)
            # Tabulated initial electron kinetic energies
            E_NIST = data['E'].astype('float64')
            # Tabulated ratios of photon emission to electron kinetic energy
            u_NIST = data['u'].astype('float64')
            # Tabulated scaled cross-section, chi
            chi_NIST = data[ 'chi' + str(self.Z) ].astype('float64')
            # Create interpolation function
            chi_fun = interp.RegularGridInterpolator( (log(E_NIST), u_NIST), 
                            chi_NIST, bounds_error = False, fill_value = 0)
            # Cross-section [cm**-2.keV**-1] calculated on E/k grid
            grid = zeros( [E_grid.shape[0], E_grid.shape[1], 2] )
            grid[:,:,0] = log(E_grid)
            grid[:,:,1] = k_grid / E_grid
            grid = squeeze(grid)
            E_in = (E_grid + self.me) / self.me
            p_in = sqrt( E_in**2 - 1. )
            sig = ( (self.Z / (p_in / E_in))**2 / k_grid) * chi_fun(grid) * 1e-27
            # Cross-section [cm**-2.keV**-1] calculated at tip (E = k)
            grid_tip = zeros([self.k.size, 2])
            grid_tip[:,0] = log(self.k)
            grid_tip[:,1] = ones( [self.k.size] )
            grid_tip = squeeze(grid_tip)
            E_in = (self.k + self.me) / self.me
            p_in = sqrt( E_in**2 - 1. )
            sig_tip = ( ( self.Z / (p_in / E_in) )**2 / self.k ) * \
                chi_fun(grid_tip) * 1e-27
        elif self.brsrc == 'classical': # Kramers classical cross-section
            # Note: we work in units of keV and cm NOT J, kg, m
            # where: me[keV/c^2] = me[kg]*c[m/s]**2/1.602e-16
            # and: hbar[keV s] = hbar[J s]/1.602e-16
            gamma_grid=(E_grid + self.me)/self.me
            betasq_grid=1-gamma_grid**(-2)
            gamma=(self.k+self.me)/self.me
            betasq=1-gamma**(-2)
            sig = (self.Z**2 / (k_grid * betasq_grid)) \
                * (16e0*pi*self.me*self.re**3/(3*sqrt(3)*self.hbar*self.c)) 
            sig_tip = (self.Z**2 / (self.k * betasq)) \
                * (16e0*pi*self.me*self.re**3/(3*sqrt(3)*self.hbar*self.c)) 
        return sig,sig_tip
 
    
    def __gen_shape_arrays(self, theta, domega_e, pomega_e, E, x):
        """Internal method to integrate electron direction distribution and 
        bremsstrahlung shape function over electron direction, for specified
        take-off angle (omega_x) 
    
        :param array theta: 1D array of angles between x-ray and electron,
            defined by omega_e and omega_x [degrees]
        :param float domega_e: element of solid angle [degrees**2] 
        :param array pomega_e: 3D array of electron prob density with angle 
            theta_e for electron energy E and depth x [degree**-1]
        :param array E: 1D array of electron kinetic energies
        :param array x: 1d array of depths for electron 
            penetration arrays [cm]
        :return array fomega_x: 3D array of bremsstrahlung integral 
            [solid angle**-1] as a function of k, E and x
        """
        # Output dependent on specified shape function
        if self.shape == 'uni':  # Uniform isotropic emission
            # Create output array
            fomega_x = zeros( [self.k.size, E.size, x.size] )
            # Generate isotropic distribution array, a 3D array for k, E and x
            fomega_x[:, :, :] = 1. / (4. * pi)
        elif self.shape == 'kqp': # Kissel-Quarles-Pratt shape function
            # Read and extract shape-function data
            filepath = os.path.join(self.dirpath, 'KQP.npz')
            data = load(filepath)
            # 1D array for angle of emission w.r.t. initial electron direction
            theta_KQP = data['theta'].astype('float64')
            # 1D array for x-ray energy, k, as fraction of electron energy, E
            u_KQP = data['u'].astype('float64')
            # 1D array for electron kinetic energy, E
            E_KQP = data['E'].astype('float64')
            # Shape-function [solid angle**-1]
            S_KQP = data['S' + str(self.Z)].astype('float64')
            # Define meshrids for interpolation
            [theta_grid, E_grid] = meshgrid( theta, E, indexing='ij' )
            # Calculate the indices (floats) corresponding to theta_grid values
            # Round to prevent floating point precision problems at boundaries
            imap = ( theta_grid - theta_KQP[0] ) \
                / ( theta_KQP[1] - theta_KQP[0] )
            imap=imap.round(decimals=3)
            # Calculate the indices (floats) corresponding to E_grid values
            # Round to prevent floating point precision problems at boundaries
            kmap = ( E_grid - E_KQP[0] ) \
                / ( E_KQP[1] - E_KQP[0] )
            kmap=kmap.round(decimals=3)
            # Create output array
            fomega_x = zeros( [self.k.size, E.size, x.size] )
            # Iterate through x-ray energy values
            for i in range(self.k.size):
                # Calculate u grid
                u_grid = self.k[i] / E_grid
                # Calculate the indices (floats) corresponding to u_grid values
                # Round to prevent float point precision problems at boundaries
                jmap = ( u_grid - u_KQP[0] ) \
                    / ( u_KQP[1] - u_KQP[0] )
                jmap=jmap.round(decimals=3)
                # 3D linear interpolation of S onto the pre-defined grid
                S = map_coordinates(S_KQP, array([imap,jmap,kmap]), order=1,
                                    cval=finfo(dtype='float64').resolution)
                # Ensure correct normalization
                S = S / ( ( sum(S, axis=0) / theta.size ) * 4. * pi )
                # Do the integral of electron distribution and shape function
                fomega_x[i, :, :] = sum( pomega_e * S[:, :, None], axis=0) \
                    * domega_e
        else: # SIM or simple shape-function
            # Define meshrids for interpolation
            [theta_grid, E_grid] = meshgrid(theta, E, indexing='ij')
            fomega_x = zeros( [self.k.size, E.size, x.size] )
            # Calculate relativistic gamma and beta
            gam = (E_grid + self.me) / self.me
            bet = sqrt( 1. - 1. / gam**2 )
            # Calculate the "SIM" shape-function
            S = (1. - bet * cos( deg2rad( theta_grid ) ) )**-2
            # Ensure correct normalization
            S = S / ( ( sum(S, axis=0) / theta.size ) * 4. * pi )
            #  Do the integral of electron distribution and shape function, 
            # ... a 3D array for k, E and x (actually independent of k and x)
            fomega_x[:,:,:] = ( sum( pomega_e * S[:,:,None], 
                                    axis=0) * domega_e )[None, :, :]
        return fomega_x    
  
    # Ref [2-3]
    def __Brems(self):
        """Internal method to calculate bremsstrahlung spectrum 
        
        :return array brem_k: 1D array of bremsstrahlung spectrum 
            [photons.keV**-1.solid-angle] 
        :return array brem_kx: 2D array of bremsstrahlung spectrum at depths
            corresponding to the x array [photons.keV**-1.cm**-1.solid-angle]
        :return array anode_self_filtration: 1D array for the product of MAC,
            density and depth for x-ray energies, k [dimensionless]
        :return array x: 1D array of depths sampled [cm]
        """ 
        # Define spectrum bin spacing in keV 
        dk = self.k[1] - self.k[0]
        
        # Get an interpolation function for log of the mass attenuation coef 
        # ... in the target as a function of log of photon energy.
        # Also gets the target density 
        logmac_target, rho_target = self.__get_logmac_fn()
        
        # Get an interpolation function for CSDA range of electrons in target
        # ... (scaled by density) as a function of electron energy
        range_fn = self.__get_range_fn() # Range function for target (g cm**-2)
        # Define the CSDA range in units of cm
        csda_range = range_fn(self.E0) / rho_target
        
        # Atomic density [atoms cm**-3]
        n = (rho_target / self.amu) * 6.0221*1e23 
      
        # Electron kinetic energy increment for sampling
        dE = 0.01 * self.E0
        # Number of energies to sample
        nE = int( round( (self.E0-dE)/dE ) ) + 1
        # The electron kinetic energies to sample
        E = linspace(dE, self.E0, nE)
      
        # Electron kinetic energy and bremsstrahlung energy mesh grids
        [E_grid, k_grid] = meshgrid(E, self.k, indexing='ij')
        # A logical array defining which combinations kinematically possible
        emission_possibility = where(k_grid/E_grid>1, 0, 1)
        # Index of the minimum possible electron energy for each brems energy
        ind_Emin = argmax(emission_possibility, axis=0)
        # Index range for bremss tip calculation
        ind_Erng = array([ind_Emin-1,ind_Emin]);
      
        # The Golden Angle in degrees
        golden_angle = 360. * ( 1. - 1. / ( ( 1.+sqrt(5.) ) / 2.) );
        # Define number of electron orientations to sample (ntheta_e+1)
        ntheta_e = 500
        # Sample points on linear scale (equiv. to equal solid angle increment)
        linear_sampling = arange(0.,ntheta_e + 0.5)
        # Convert sampling points to theta_e values
        theta_e  = degrees(arccos(1.0 - 2.0 * linear_sampling / ntheta_e))
        # Calculate a unique phi_e angle for each theta_e, using golden angle
        phi_e    = linear_sampling * golden_angle
        # Calculate direction-cosine vector for all electron orientations
        omega_e  = array( [ cos( deg2rad(theta_e) ), 
                    sin( deg2rad(theta_e) ) * sin( deg2rad(phi_e) ),
                    -sin( deg2rad(theta_e)) * cos( deg2rad(phi_e) ) ]).T
        # Solid angle element for each electron orientation
        domega_e = 360. * (2. * 180. / pi) / float(ntheta_e + 1);
        
        # X-ray emission direction vector defined by take-off angles
        omega_x = array( [ sin(-self.varphi)*cos(self.vartheta), 
                              sin(self.vartheta), 
                              cos(self.vartheta)*cos(-self.varphi) ] )
        
        # Emission angle relative to initial electron direction (degrees)
        theta = degrees( arccos( matmul(omega_e, omega_x) ).real ) 
        
        # Generate arrays for:
        # * The depths, x (NB: not uniformly spaced) [1d]
        # * The electron frequency distribution at the depths, x, differential 
        # ... in electron energies, E [2d]
        # * The electron ang. dist. for electron energies, E, at depths, x [3d]
        x, NE, pomega_e = self.__gen_epen_arrays(csda_range, E, theta_e, 
                                                 domega_e)
        
        # Generate arrays for:
        # * The brems cross-section dsig/dk for sampled E/k grid [2d]
        # * The brems cross-section dsig/dk for the tip E = k [1d]
        sig, sig_tip = self.__gen_sig_arrays(E_grid, k_grid)
        
        # Generate arrays for:
        # * The integral of (pomega_e * brems shape fn) over omega_e 
        # This returns a value corresponding to directio omega_x, sampled at
        # ... values of k, E and x [3d]
        fomega_x = self.__gen_shape_arrays(theta, domega_e, pomega_e, E, x) 
        
        mac_sub = zeros( [11] )
        N_brems = zeros( self.k.shape )
        anode_filtration = zeros( [self.k.size, x.size] )
        brem_kx = zeros( [self.k.size ,x.size] )
        # Iterate over x-ray emission energy
        for i in range(self.k.size):
            # Iterate over subsamples of x-ray energy in the bin
            for j in range(mac_sub.size):
                # mac values at subsampled energy
                mac_sub[j] = exp( logmac_target( 
                    log( self.k[i]-dk*0.5+(dk*0.1)*float(j) )
                    ) )
            # Average attenuation for the x-ray energy bin for each depth
            attn_anode = exp( -mac_sub.mean() * rho_target * x * 
                            sin(self.varphi)**(-1) * cos(self.vartheta)**(-1) )
            # Average filtration for each x-ray energy bins and depth
            anode_filtration[i,:] = mac_sub.mean() * rho_target * x
            # Generate interpolation fn for electron freq. distribution at the 
            # ... tip (the tip region is where k ~ E)
            # The function returns a vector of values for x, for each E value
            NE_tip = interp.interp1d( E[ ind_Erng[:,i] ], NE[ ind_Erng[:,i] ].T
                        , kind='linear', bounds_error=False, fill_value=0.0)
            # Generate interpolation fn for angular distribution at the 
            # ... tip (the tip region is where k ~ E)
            # The function returns a vector of values for x, for each E value
            fomega_x_tip = interp.interp1d( E[ ind_Erng[:,i] ], 
                            squeeze( fomega_x[i,ind_Erng[:,i],:] ).T,
                            kind='linear', bounds_error=False, fill_value=0.0)
            # The number of emissions per keV per solid angle at depths, x,
            # ... differential in electron energy, E
            N_brems_xE = n * NE * (sig[:,i])[:,None] * \
                squeeze( fomega_x[i,:,:] )
            # The number of emissions per keV per solid angle at depths, x
            # Sum of (1) integral contribution from lowest physical sampling 
            # ... point (E[ind_Emin[i]]) up to max (E[-1]) and 
            # (2) tip region from k to lowest sampled point (E[ind_Emin[i]])
            N_brems_x = trapz( N_brems_xE[ ind_Emin[i]: ], E[ ind_Emin[i]: ]
                            , axis=0 ) \
                      - trapz( [ N_brems_xE[ ind_Emin[i],: ], 
                            n * NE_tip( self.k[i] ) * sig_tip[i] * 
                            fomega_x_tip( self.k[i] ) ], 
                                [ E[ ind_Emin[i] ], self.k[i] ] 
                            , axis=0 )
            # Number of emissions per keV per solid angle including self-attn       
            N_brems[i] = trapz( N_brems_x * attn_anode, x )
            # Store the emissions for each depth
            brem_kx[i,:] = N_brems_x
      
        # Assign return variables based on the physical model
        # 'kqp and 'sim' mode return anode-filtered data summed over depths
        # Other models return un-filtered data un-summed over depths and the 
        # ... filtration for the depths and energies
        if self.shape == 'kqp' or self.shape == 'sim':
            brem_k = N_brems
            brem_kx = None
            anode_self_filtration = None
        elif self.shape == 'casim' or self.shape == 'uni':
            brem_k = None
            brem_kx = brem_kx
            anode_self_filtration = anode_filtration
        else:
            raise Exception('Physics model not recognized!')
          
        return brem_k, brem_kx, anode_self_filtration, x


    def __gen_char_dists(self,x):
        """Internal method to calculate interpolated depth-distributions for 
        characteristic spectrum 
    
        :param x: 1D array of depths to sample (as fraction of CSDA range) 
        :return array Kdist: 1D array of K characteristic emissions at 
            depths x, per incident electron
        :return array L1dist: 1D array of L1 characteristic emissions at 
            depths x, per incident electron
        :return array L2dist: 1D array of L2 characteristic emissions at 
            depths x, per incident electron
        :return array L3dist: 1D array of L3 characteristic emissions at 
            depths x, per incident electron
        """
        # Read and extract the tabulated data
        filepath = os.path.join(self.dirpath, 'FluorescenceZ' 
                                + str(self.Z) + '.npz')
        data = load(filepath)
        # K depths tabulated
        x_K_data = data['x_K'].astype('float64')
        # K incident energies tabulated
        E0_K_data = data['E0_K'].astype('float64')
        # K emission data
        K_data =  data['K'].astype('float64')
        # Assign the L shell data to arrays if it exists
        try:
            x_L1_data = data['x_L1'].astype('float64')
            E0_L1_data = data['E0_L1'].astype('float64')
            x_L2_data = data['x_L2'].astype('float64') 
            E0_L2_data = data['E0_L2'].astype('float64') 
            x_L3_data = data['x_L3'].astype('float64') 
            E0_L3_data = data['E0_L3'].astype('float64') 
            L1_data =  data['L1'].astype('float64')
            L2_data =  data['L2'].astype('float64')
            L3_data =  data['L3'].astype('float64')
            # L data present
            L_included = True
        except:
            # L data absent
            L_included = False
        # Create interpolation function for K depth-data depending on x and E0
        Kfun = interp.RegularGridInterpolator( (x_K_data, E0_K_data), K_data, 
            bounds_error=False, fill_value=finfo(dtype='float64').resolution)
        # Create interpolation function for L depth-data (if available) 
        # ... depending on x and E0
        if L_included:
            L1fun = interp.RegularGridInterpolator( (x_L1_data, E0_L1_data), 
                L1_data, bounds_error=False, 
                fill_value=finfo(dtype='float64').resolution)
            L2fun = interp.RegularGridInterpolator( (x_L2_data, E0_L2_data), 
                L2_data, bounds_error=False, 
                fill_value=finfo(dtype='float64').resolution)
            L3fun = interp.RegularGridInterpolator( (x_L3_data, E0_L3_data), 
                L3_data, bounds_error=False,
                fill_value=finfo(dtype='float64').resolution)
        # Create meshgrids for interpolation function
        [x_grid,E0_grid] = meshgrid(x, self.E0, indexing='ij')
        # Create variable array to hold the grid data
        var = zeros( [x_grid.shape[0], x_grid.shape[1], 2] )
        var[:,:,0] = x_grid
        var[:,:,1] = E0_grid
        # Interpolate for depths, x, and the E0 value. Remove singlet dimension
        # ... as interpolation is only done for 1 value of E0
        Kdist = squeeze( Kfun(var) )
        # Do interpolations for L shells, if data available
        if L_included:
            L1dist = squeeze( L1fun(var) )
            L2dist = squeeze( L2fun(var) )
            L3dist = squeeze( L3fun(var) )
        else:
            L1dist = zeros(Kdist.shape)
            L2dist = zeros(Kdist.shape)
            L3dist = zeros(Kdist.shape)
        return Kdist,L1dist,L2dist,L3dist


    def __get_char_data(self):
        """Internal method to get K and L shell data from data file
    
        :return float UK: K-edge energy [keV]
        :return float UL1: L1-edge energy [keV]
        :return float UL2: L2-edge energy [keV]
        :return float UL3: L3-edge energy [keV]
        :return array kK: K-line energies [keV]
        :return float kL1: L1-line energies [keV]
        :return float kL2: L2-line energies [keV]
        :return float kL3: L3-line energies [keV]
        :return array PK: K-line probabilities per K characteristic emission
        :return float PL1: L1-line probabilities per L1 characteristic emission
        :return float PL2: L2-line probabilities per L2 characteristic emission
        :return float PL3: L3-line probabilities per L3 characteristic emission
        """
        # Read an extract data from file
        filepath = os.path.join(self.dirpath, 'Z' + str(self.Z) + 'ch.npz')
        data = load(filepath)
        UK = data['UK'].astype('float64')
        # Normalize per K-line emission rather than K relaxation
        # ... (the remainder in un-normlized data is Auger electron contrib.)
        PK = data['PK'].astype('float64') / sum(data['PK'].astype('float64'))
        kK = data['kK'].astype('float64')
        # Extract L shell data if present
        try:
            UL1 = data['UL1'].astype('float64')
            UL2 = data['UL2'].astype('float64')
            UL3 = data['UL3'].astype('float64')
            kL1 = data['kL1'].astype('float64')
            kL2 = data['kL2'].astype('float64')
            kL3 = data['kL3'].astype('float64')
            # Normalize per L-line emission rather than L relaxation (the 
            # ... remainder in un-normlized data is Auger electron contrib.)
            PL1 = data['PL1'].astype('float64') \
                / sum(data['PL1']).astype('float64')
            PL2 = data['PL2'].astype('float64') \
                / sum(data['PL2'].astype('float64'))
            PL3 = data['PL3'].astype('float64') \
                / sum(data['PL3'].astype('float64'))
            return UK, kK, PK, UL1, UL2, UL3, kL1, kL2, kL3, PL1, PL2, PL3
        except:
            return UK, kK, PK, nan, nan, nan, nan, nan, nan, nan, nan, nan

    # Ref. [1]
    def __Char(self):
        """Internal method to calculate characteristic spectrum 
    
        :return array char_k: 1D array of characteristic spectrum 
            [photons.keV**-1.solid-angle] 
        :return array char_kx: 2D array of characteristic spectrum at depths
            corresponding to the x array [photons.keV**-1.cm**-1.solid-angle]
        :return array anode_self_filtration_char: 1D array for the product of 
            MAC, density and depth for x-ray energies, k [dimensionless]
        :return array x: 1D array of depths sampled [cm]
        """ 
        # Define spectrum bin spacing in keV 
        dk = self.k[1] - self.k[0]    
      
        # Get characteristic data for K and L edges 
        # (edge energies [U], line energies [k], probabilities [P])
        UK, kK, PK, UL1, UL2, UL3, kL1, kL2, kL3, PL1, PL2, PL3 = \
            self.__get_char_data()
    
        # Get an interpolation function for log of the mass attn coef in the 
        # ... target as a fn of log of photon energy. Also gets target density       
        mac_target, rho_target = self.__get_logmac_fn(Rayleigh=False)
      
        # Get an interpolation function for CSDA range of electrons in target
        # ... (scaled by density) as a function of electron energy
        range_target = self.__get_range_fn()
        # Define the CSDA range in units of cm
        csda_range = range_target(self.E0) / rho_target
        # Define the maximum electron penetration depth considered
        dmax = 20.*csda_range 
        
        # Define number of depths to sample
        if self.shape == 'kqp' or self.shape == 'sim':
            N = 1500 # High number as still fast runtime compared to brems part
        else:
            N = 50 # Lower as otherwise execution time becomes similar to brems
          
        # Depth sampling spacing  
        dx = csda_range / float(N)
        # Array of depths to sample
        x = arange(0,dmax+0.5*dx,dx)
      
        # Get the frequency distribution of characteristic emissions for target
        # ... for the incident kinetic electron energy, E0 and depths, x
        Kdist, L1dist, L2dist, L3dist = self.__gen_char_dists(x/csda_range)
        
        # Intitialize arrays
        char_k = zeros(self.k.size)
        mu_times_x = zeros([self.k.size,x.size])
        char_kx = zeros([self.k.size,x.size])
        ## Self-filtration of x-rays produced in the target

        if self.E0 > UL1: # If L1-shell ionization is possible, do calculation
            NL1 = zeros([kL1.size])
            for i in range(kL1.size):
                # Attenuation factor in escape target, considering photon 
                # ... energy and the take-off angle
                attn_anode = exp( -exp( mac_target(log(kL1[i]) ) )
                        * rho_target * x 
                        * sin(self.varphi)**-1 * cos(self.vartheta)**-1 )
                # Number of characteristic photons in the line per solid angle
                NL1[i] = (1./(4.*pi)) * trapz(L1dist*attn_anode,x/csda_range)
                # The bin index of the line
                ind = floor( ( kL1[i] - (self.k[0] - dk*0.5) ) 
                            / dk ).astype(int)
                if ind<0:
                    continue
                # Adds characteristic photons per solid angle per keV to bin
                char_k[ind] = char_k[ind] + NL1[i] * PL1[i] / dk
                # Adds characteristic photons depth distrubution per solid 
                # ... angle per keV to the bin
                char_kx[ind,:] = char_kx[ind,:] + \
                    (1./(4.*pi)) * (L1dist * attn_anode / csda_range) \
                        * PL1[i] / dk
    
        if self.E0 > UL2: # If L2-shell ionization is possible, do calculation
            NL2 = zeros([kL2.size])
            for i in range(kL2.size):
                # Attenuation factor in escape target, considering photon 
                # ... energy and the take-off angle
                attn_anode = exp( -exp( mac_target(log(kL2[i]) ) ) 
                            * rho_target * x 
                            * sin(self.varphi)**-1 * cos(self.vartheta)**-1 )
                # Number of characteristic photons in the line per solid angle
                NL2[i]= (1./(4.*pi)) * trapz(L2dist*attn_anode,x/csda_range)
                # The bin index of the line
                ind = floor( ( kL2[i] - (self.k[0] - dk*0.5) ) 
                            / dk ).astype(int)
                if ind<0:
                    continue
                # Adds characteristic photons per solid angle per keV to bin
                char_k[ind] = char_k[ind] + NL2[i] * PL2[i] / dk
                # Adds the characteristic photons depth distrubution per solid 
                # ... angle per keV to the bin
                char_kx[ind,:] = char_kx[ind,:] + \
                    (1./(4.*pi)) * (L2dist * attn_anode / csda_range) \
                        * PL2[i] / dk
    
        if self.E0 > UL3: # If L3-shell ionization is possible, do calculation
            NL3 = zeros([kL3.size])
            for i in range(kL3.size):
                # Attenuation factor in escape target, considering photon 
                # ... energy and the take-off angle
                attn_anode = exp( -exp( mac_target(log(kL3[i]) ) ) 
                            * rho_target * x
                            * sin(self.varphi)**-1 * cos(self.vartheta)**-1 )
                # Number of characteristic photons in the line per solid angle
                NL3[i]= (1./(4.*pi)) * trapz(L3dist*attn_anode,x/csda_range)
                # The bin index of the line
                ind = floor( ( kL3[i] - (self.k[0] - dk*0.5) ) \
                            / dk ).astype(int)
                if ind<0:
                    continue
                # Adds characteristic photons per solid angle per keV to bin
                char_k[ind] = char_k[ind] + NL3[i] * PL3[i] / dk
                # Adds the characteristic photons depth distrubution per solid
                # ... angle per keV to the bin
                char_kx[ind,:] = char_kx[ind,:] + \
                    (1./(4.*pi)) * (L3dist * attn_anode / csda_range) \
                        * PL3[i] / dk
                
        if self.E0 > UK: # If K-shell ionization is possible, do calculation
            NK = zeros([kK.size])
            for i in range(kK.size):
                # Attenuation factor in escape target, considering photon 
                # ... energy and the take-off angle
                attn_anode = exp( -exp( mac_target(log(kK[i]) ) ) 
                            * rho_target * x 
                            * sin(self.varphi)**-1 * cos(self.vartheta)**-1 )
                # Number of characteristic photons in the line per solid angle
                NK[i]= (1./(4.*pi)) * trapz(Kdist*attn_anode,x/csda_range)
                # The bin index of the line
                ind = floor((kK[i] - (self.k[0] - dk * 0.5)) \
                            / dk).astype(int)
                if ind<0:
                    continue
                # Adds characteristic photons per solid angle per keV to bin
                char_k[ind] = char_k[ind] + NK[i]*PK[i]/dk
                # Adds the characteristic photons depth distrubution per solid
                # ... angle per keV to the bin
                char_kx[ind,:] = char_kx[ind,:] + \
                    (1./(4.*pi)) * (Kdist * attn_anode / csda_range) \
                        * PK[i] / dk

        for i in range(self.k.size):
            # linear attenuation coef. times depth for energies and depths
            mu_times_x[i,:] = exp( mac_target( log(self.k[i]) ) ) \
                * rho_target * x
            # Attenuation factor to escape target, for the photon energy and 
            # ... the depth bins for the take-off angle
            attn_factor =  exp( -mu_times_x[i,:] * 
                               sin(self.varphi)**-1 * cos(self.vartheta)**-1 )
            # Characteristic freq. dist. per solid angle per keV, with the
            # ... self-filtration removed (applied again outside this class)
            char_kx[i,:] = divide(char_kx[i,:],attn_factor,where=attn_factor>0)
    
        if self.shape == 'kqp' or self.shape == 'sim':
            char_k = char_k
            char_kx = None
            anode_self_filtration_char = None
        elif self.shape == 'casim' or self.shape == 'uni':
            char_k = None
            char_kx = char_kx
            anode_self_filtration_char = mu_times_x
        else:
            raise Exception('Physics model not recognized!')
        
        return char_k, char_kx, anode_self_filtration_char, x



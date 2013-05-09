#!/usr/bin/env python

"""
The AMPX module contains a class for reading AMPX cross section libraries, in
*working* format.  The AMPX library format is widely used within
the SCALE package of ORNL.  The format is specified in the SCALE documentation,
but also in Appendix A of
 
  Ford, W. E. et al "ANSL-V: ENDF/B-V Based Multigroup Cross-Section Libraries
    for Advanced Neutron Source (ANS) Reactor Studies", Technical Report
    ORNL-6618 (1990)

which is available on the web.  AMPX comes in two flavors: master and 
working.  The master format is complex and holds a lot of data specific to
self-shielding.  The working format is more straightforward and is what one
gets during SCALE analysis sequences.  For example, one can use the SCALE
modules nitawl and ice to create isotopic and macroscopic libraries, 
respectively.
"""

import struct
import math

from binaryreader import _BinaryReader, _FortranRecord

class AMPX(_BinaryReader):
    """ An AMPX object represents a binary AMPX working library 

    :Attributes:
      **header** : dict
        Dictionary of all header information
      **number_groups_n** : int
        Number of neutron energy groups
      **number_groups_p** : int
        Number of photon energy groups
      **number_groups** : int
        Total number of neutron and photon energy groups
      **number_nuclides** : int
        Number of nuclides (or mixtures) in the library
      **n_eb** : list of doubles
        Neutron energy bounds
      **p_eb** : list of doubles
        Photon energy bounds
      **data** : data
        List of nuclide data dictionaries in the order read
      
    Parameters
    ----------
    filename : str
        Path of the AMPX file to load.
    endian : str
        Endian marker (may be needed for cross-platform compatibility)
    """

    def __init__(self, filename, endian = "@"):
        super(AMPX, self).__init__(filename, 'rb', endian)

        self.header = {}            
        self.number_groups_n = 0    
        self.number_groups_p = 0    
        self.number_groups = 0
        self.data = []           
        self.master = False # flag for eventual handling of master libraries
                
    def read(self) :
        """ Read the binary file.
        """

        # Read the header.  Keep some common parameters as attributes.
        self.header = self._read_record_type_1()
        self.number_groups_n = self.header['igm'] 
        self.number_groups_p = self.header['ipm']
        self.number_groups   = self.number_groups_n + self.number_groups_p
        self.number_nuclides = self.header['nnuc']
        
        # Read energy and lethargy boundaries
        if self.number_groups_n > 0:
          self.n_eb, self.n_ub = self._read_record_type_2(self.number_groups_n)
        if self.number_groups_p > 0:
          self.p_eb, self.p_ub = self._read_record_type_2(self.number_groups_p)
        
        # Read cross section records.  Note, this data is actually repeated
        # in the library.  The first pass gets the required bounds, etc.
        for n in range(0, self.number_nuclides):
            # Initialize dictionary for this nuclide
            data = {}
            # All library nuclides start with record type 3
            self._read_record_type_3(data)
            self.data.append(data)

        # Read cross section data
        for n in range(0, self.number_nuclides):
           
            # Re-read the meta-data.  Opting not to check against the former.
            tmp = {} 
            self._read_record_type_3(tmp)
            
            # Read the cross sections
            data = self.data[n]
            if self.master :
                pass
            else :
                # 1-D cross sections
                #   Neutrons
                if self.number_groups_n > 0:
                    data['sigma_n'] =                                          \
                        self._read_record_type_9(data['num_1d_n_reactions'],   \
                                                 self.number_groups_n)
                #   Photons
                if self.number_groups_p > 0:
                    data['sigma_p'] =                                          \
                        self._read_record_type_9(data['num_1d_p_reactions'],   \
                                                 self.number_groups_p)
                # 2-D cross sections
                #   Note, both neutrons and photons read together.  This 
                #   produces a (ng_n + ng_p)^2 matrix.  
                order = data['scatter_expansion_order']
                # One of both of these is nonzero if scattering matrix present
                length = max(data['max_length_p0_array'], \
                             data['length_total_scatter'])
                if length > 0 :
                    data['scatter'] = \
                        self._read_record_type_12(order, self.number_groups)
                            
    def _read_record_type_1(self):
        """ Read the header.
  
        This defines the number of isotopes, energy groups, etc.  Note, the
        somewhat terse keys are straight from the documentation.
        """
        # Get file control record
        fc = self.get_fortran_record()

        header = {}
        # Identification number
        header['idtape'] = fc.get_int()[0]       
        # Number of sets of data
        header['nnuc']   = fc.get_int()[0]       
        # Number of neutron energy groups
        header['igm']    = fc.get_int()[0]       
        # First thermal neutron group (i.e. first that gets upscatter source)
        header['iftg']   = fc.get_int()[0]       
        # Master library version type (only for working)
        header['msn']    = fc.get_int()[0]       
        # Number of photon energy groups 
        header['ipm']    = fc.get_int()[0]       
        # Not used
        header['I1']     = fc.get_int(1)[0]      
        # Produced by weighting a working library in XSDRNPM [0/1 = no/yes]
        header['I2']     = fc.get_int(1)[0]
        # Not used      
        header['I3']     = fc.get_int(1)[0]
        # Not used
        header['I4']     = fc.get_int(1)[0]  
        # Library title
        header['title']  = fc.get_string(400)[0]

        return header
    
    def _read_record_type_2(self, ng):
        """ Read the energy and lethargy boundaries.  
       
        Note, the lethargy u = ln(E_0/E) is defined using E_0 = 10 MeV.
        
        Parameters
        ----------
        ng : int
            Number of groups for which boundaries are to be read
        """
        fc = self.get_fortran_record()
        eb = fc.get_float(ng+1)
        ub = fc.get_float(ng+1)
        return eb, ub

    def _read_record_type_3(self, data):
        """ Read a cross section set directory 
        
        This contains essentially meta-data for the actual cross section 
        data read from records 4 and up.
        
        Parameters
        ----------
        data : dict
            Dictionary of data for an isotope or mixture
        """
        
        # Read the directory record.  While the record size is the same for 
        # master and working libraries, the contents have different meanings.
        # Here, we'll use generic names where meaning differs.
        fc = self.get_fortran_record()
        data['description'] = fc.get_string(18*4)[0]            #  1-18
        data['id'] = fc.get_int()[0]                            # 19
        if self.master:
            data['num_resonance_sets'] = fc.get_int()[0]        # 20
            data['num_unresolved_energies'] = fc.get_int()[0]   # 21
            data['num_1d_n_reactions'] = fc.get_int()[0]        # 22
            data['num_2d_n_reactions'] = fc.get_int()[0]        # 23
            data['i24'] = fc.get_int()[0]                       # 24   
            data['num_1d_p_reactions'] = fc.get_int()[0]        # 25
            data['num_2d_p_reactions'] = fc.get_int()[0]        # 26 
            data['num_2d_np_reactions'] = fc.get_int()[0]       # 27   
            data['i28'] = fc.get_int()[0]                       # 28
        else:
            data['id_parent_set'] = fc.get_int()[0]             # 20
            data['zone_number'] = fc.get_int()[0]               # 21
            data['num_zones'] = fc.get_int()[0]                 # 22
            data['length_total_scatter'] = fc.get_int()[0]      # 23
            data['scatter_expansion_order'] = fc.get_int()[0]   # 24   
            data['sequence'] = fc.get_int()[0]                  # 25
            data['num_sets_nuclide'] = fc.get_int()[0]          # 26 
            data['max_length_p0_array'] = fc.get_int()[0]       # 27   
            data['num_1d_n_reactions'] = fc.get_int()[0]        # 28            
        data['mass_number'] = fc.get_float()[0]                 # 29
        data['zaid'] = fc.get_float()[0]                        # 30
        for i in range(31, 34):                         
            data['i'+str(i)] = fc.get_int()[0]                  # 31-33
        data['power_per_fission'] = fc.get_float()[0]           # 34 [W-s/fiss]
        data['energy_per_fission'] = fc.get_float()[0]          # 35 [W-s/cap]
        if self.master:
            data['max_length_2d'] = fc.get_int()[0]             # 36            
            data['num_bondarenko_sets'] = fc.get_int()[0]       # 37
            data['num_bondarenko_sigma_0'] = fc.get_int()[0]    # 38
            data['num_bondarenko_T'] = fc.get_int()[0]          # 39    
            data['num_bondarenko_groups'] = fc.get_int()[0]     # 40  
            data['i41'] = fc.get_int()[0]                       # 41
        else:
            for i in range(36, 41): 
                data['i'+str(i)] = fc.get_int()[0]              # 36-40
            data['num_1d_p_reactions'] = fc.get_int()[0]        # 41
        data['i42'] = fc.get_int()[0]                           # 42
        if self.master:
            data['sigma_p'] = fc.get_float()[0]                 # 43
        else:
            data['i43'] = fc.get_float()[0]                     # 43
        data['i44'] =  fc.get_int()[0]                          # 44
        data['endf_mat_fast_n'] = fc.get_int()[0]               # 45
        data['endf_mat_thermal_n'] = fc.get_int()[0]            # 46
        data['endf_mat_p'] = fc.get_int()[0]                    # 47
        data['endf_mat_p_production'] = fc.get_int()[0]         # 48
        data['nuclide'] = fc.get_string(4)[0]                   # 49
        data['number_records'] = fc.get_int()[0]                # 50    
            
    def _read_record_type_4(self, data):
        """ Read resonance parameters
        """
        pass
    
    def _read_record_type_5(self, data):
        """ Read first record of Bondarenko block 
        """
        pass
    
    def _read_record_type_6(self, data):
        """ Read directory for Bondarenko data
        """
        pass
    
    def _read_record_type_7(self, data):
        """ Read infinite dilution values for Bondarenko data
        """
        pass   
    
    def _read_record_type_8(self, data):
        """ Read Bondarenko factors
        """
        pass     
    
    def _read_record_type_9(self, num_mt, num_g):
        """ Read temperature-independent average cross sections
        
        These are sometimes called the one-dimensional cross sections.  The
        structure is 
          MT1 (sigma_1(I), I=1, IGM)
          MT2 (...)
            ...
        The MTx parameters identify the reaction and are stored as *floats*
        
        Parameters
        ----------
        num_mt : int
            Number of reaction types 
        num_g : int
            Number of groups 
        """
        fc = self.get_fortran_record()
        sigma = {}
        for m in range(0, num_mt):
            # Get the endf reaction identifier 
            mt = int(fc.get_float()[0])        
            # Read in the associated cross section data
            sigma[mt] = fc.get_float(num_g)
        return sigma

    def _read_record_type_12(self, order, num_g):
        """ Read scattering matrix (and other 2-d data)
        

        Parameters
        ----------
        order : int
            Order of Legendre expansion of the scattering matrix
        num_g : int
            Number of groups to read in the record.  For working libraries,
            this can be the total of neutron and photon groups.
        """

        # The AMPX format keeps scattering in a sparse matrix format 
        # based on the "magic" number of the form
        #
        # JJJKKKIII where
        #    III = sink group
        #    JJJ = lower bound of groups entering group
        #    KKK = upper bound of groups entering group 
        #
        # The magic number is stored, but the full scattering matrix is 
        # kept in memory for easier manipulation.  The astute reader will
        # recognize this limits AMPX libraries to less than 1000 groups.
        
        # List of scattering matrices and their magic numbers
        scatter = []
        
        # Loop over all the Legendre orders
        for l in range(0, order + 1):       
            
            fc = self.get_fortran_record()  # Each order is a new record
            L  = fc.get_int()[0]            # Total size of the record
            S = []                          # Scattering matrix for this order
            magic = []                      # Magic numbers for this order

            for g in range(0, num_g):
               
                magic.append(fc.get_int()[0])
                lower = (magic[g] / 1000000) - 1
                upper = (magic[g] - 1000000 * (lower + 1)) / 1000
                assert(magic[g] % 100 == g + 1)
              
                # Get the matrix row, reversing to get high E -> low E
                S.append([0] * num_g)
                S[g][lower:upper] = fc.get_float(upper-lower)[::-1]
            
            scatter.append([S, magic])      
                  
        return scatter

            

#
# Copyright 2015 Olli Tapaninen, VTT Technical Research Center of Finland
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

import numpy as np
from scipy.interpolate import griddata
import warnings
from pkg_resources import resource_string
import io


class Illuminants(object):

    def __init__(self):
        '''
        Constructor
        '''
        # Some physical constants
        # Boltzmann constant
        self.kb = 1.3806488e-23
        # Speed of light
        self.c = 299792458.0
        # Planck constant
        self.h = 6.62606957e-34

        # Load S_n(lambda) for D-illuminant
        data = resource_string(__name__,  'data/SforD-illuminant.txt')

        dtypes = {'names':   ['wavelen',   'S0',  'S1',  'S2'],
                  'formats': [float,      float, float, float]
                  }
        self.D_S_Data = np.loadtxt(io.BytesIO(data),
                                   dtype=dtypes)

        # Load Illuminant F Spectrums
        data = resource_string(__name__,  'data/IlluminantF.txt')

        dtypes = {'names':   ['wavelen',  'F01', 'F02',
                              'F03', 'F04', 'F05', 'F06', 'F07',
                              'F08', 'F09', 'F10', 'F11', 'F12'
                              ],
                  'formats': [float, float, float, float, float,
                              float, float, float, float, float,
                              float, float, float]
                  }
        self.illumF = np.loadtxt(io.BytesIO(data),
                                 dtype=dtypes)

    def blackbodySpectrum(self, wavelen, T_K):
        '''
        Returns the spectrum of a black body at a temperature of T_K
        for given wavelengths [nm]
        '''
        T_K = float(T_K)
        return (2.0 * np.pi * self.h * self.c ** 2) /\
            ((wavelen * 10 ** -9) ** 5) *\
            (1.0 / (np.exp((self.h * self.c) /
                           (wavelen * 10 ** -9 * self.kb * T_K)) - 1.0))

    def blackbodySpectrum560(self, wavelen, T_K):
        '''
        Returns the spectrum of a black body at a temperature of T_K
        for given wavelengths [nm]. Normalized to 1 at 560 nm.
        '''
        T_K = float(T_K)

        le = (wavelen * 10 ** -9)**-5 * \
            (np.exp(1.4388 * 10**-2 / ((wavelen * 10 ** -9) * T_K)) - 1)**-1

        le560 = (560.0 * 10 ** -9)**-5 * \
            (np.exp(1.4388 * 10**-2 / ((560.0 * 10 ** -9) * T_K)) - 1)**-1
        return le / le560

    def illuminantA(self, wavelen):
        '''
        Returns CIE standard illuminant A for given vavelengths.
        SPD 100 at 560nm
        '''
        T_K = 2848.0
        return (100.0 * (560.0 / wavelen) ** 5 *
                (np.exp((1.435e7) / (T_K * 560.0)) - 1.0) /
                (np.exp((1.435e7) / (T_K * wavelen)) - 1.0)
                )

    def illuminantB(self, wavelen):
        '''
        B served as a representative of noon sunlight at T = 4874K

        This is not the official standard illuminant B, but a
        illuminant D at 4874K, of which the B approximates.
        '''
        return self.illuminantD(wavelen, 4874.0)

    def illuminantC(self, wavelen):
        '''
        C represented average day light with a CCT of 6774 K

        This is not the official standard illuminant C, but a
        illuminant D at 6774K, of which the C approximates.
        '''
        return self.illuminantD(wavelen, 6774.0)

    def illuminantD(self, wavelen, T_K):
        '''
        Returns the spectral power distribution of CIE standard illuminant D
        for given wavelengths and T_K.

        wavelen in nm

        See, http://en.wikipedia.org/wiki/Standard_illuminant
        '''
        T_K = float(T_K)

        # Calculate chromaticity coordinates
        if T_K >= 4000 and T_K <= 7000:
            xd = (0.244063
                  + 0.09911 * 10 ** 3 / T_K
                  + 2.9678 * 10 ** 6 / T_K ** 2
                  - 4.6070 * 10 ** 9 / T_K ** 3
                  )
        elif T_K >= 7000 and T_K <= 25000:
            xd = (0.237040
                  + 0.24748 * 10 ** 3 / T_K
                  + 1.9018 * 10 ** 6 / T_K ** 2
                  - 2.0064 * 10 ** 9 / T_K ** 3
                  )
        else:
            warnings.warn(
                'Illuminant D is not defined for T = {} K'.format(T_K))
            xd = (0.237040
                  + 0.24748 * 10 ** 3 / T_K
                  + 1.9018 * 10 ** 6 / T_K ** 2
                  - 2.0064 * 10 ** 9 / T_K ** 3
                  )

        yd = -3.000 * xd ** 2 + 2.870 * xd - 0.275

        M = 0.0241 + 0.2562 * xd - 0.7341 * yd
        M1 = (-1.3515 - 1.7703 * xd + 5.9114 * yd) / M
        M2 = (0.03000 - 31.4424 * xd + 30.0717 * yd) / M

        S0 = griddata(self.D_S_Data['wavelen'],
                      self.D_S_Data['S0'], wavelen, fill_value=0.0)
        S1 = griddata(self.D_S_Data['wavelen'],
                      self.D_S_Data['S1'], wavelen, fill_value=0.0)
        S2 = griddata(self.D_S_Data['wavelen'],
                      self.D_S_Data['S2'], wavelen, fill_value=0.0)

        S = S0 + M1 * S1 + M2 * S2

        return S

    def illuminantE(self, wavelen):
        '''
        Returns the spectrum of CIE Illuminant E.
        Constant value 100 at all wavelengths.


        wavelen in nm
        '''
        return np.ones(len(wavelen), dtype=float)

    def illuminantF(self, wavelen, F=1):
        '''
        Returns the spectrum of CIE Illuminant F 1-12

        wavelen in nm
        '''
        return griddata(self.illumF['wavelen'],
                        self.illumF['F{:{fill}2d}'.format(F, fill=0)],
                        wavelen)

    def illuminantD_M(self, wavelen, T_K):
        '''
        Returns the spectral power distribution of CIE standard illuminant D
        for given wavelengths and T_K proportionaly mixed with
        blackbody spectrum.

        Available for 4500K < T_K < 5500K

        wavelen in nm

        See, TM-30-15
        '''
        if not (T_K > 4500 and T_K < 5500):
            raise RuntimeError(
                'Illuminant D_M only available for 4500K < T_K < 5500K, not %d'
                % T_K)

        S_P = self.blackbodySpectrum560(wavelen, T_K)
        print('----')
        print(np.trapz(y=S_P, x=wavelen))
        S_D = self.illuminantD(wavelen, T_K)
        print(np.trapz(y=S_D, x=wavelen))
        print('----')

        S_M = (5500.0 - T_K) / 1000.0 * S_P + \
              (1.0 - (5500.0 - T_K) / 1000.0) * S_D
        return(S_M)

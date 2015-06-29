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

import io
import warnings

from pkg_resources import resource_string
from scipy.interpolate import griddata

import numpy as np

from . import ColorTransforms, Illuminants


warnings.simplefilter("always")


class CIEXYZ(object):

    def __init__(self):
        '''
        Constructor
        '''
        self.ct = ColorTransforms.ColorTransforms()
        self.illuminants = Illuminants.Illuminants()

        ######################################
        data = resource_string(__name__, "data/CIE_1931_STDcolorObs.txt")

        # Load CIE 1931 2deg standard observer data
        dtypes = {'names':   ['wavelen',   'x',  'y',  'z'],
                  'formats': [float,      float, float, float]
                  }
        self.observer2Data1931 = np.genfromtxt(io.BytesIO(data),
                                               dtype=dtypes)
        # print(self.observer2Data1931)

        # Load CIE 1964 10deg standard observer data
        data = resource_string(__name__, "data/CIE_1964_STDcolorObs.txt")

        dtypes = {'names':   ['wavelen',   'x',  'y',  'z'],
                  'formats': [float,      float, float, float]
                  }
        self.observer10Data1964 = np.loadtxt(io.BytesIO(data), dtype=dtypes)

        ######################################
        data = resource_string(__name__, "data/CIE_2006_STDcolorObs_2deg.txt")

        # Load CIE 2006 2deg standard observer data
        dtypes = {'names':   ['wavelen',   'x',  'y',  'z'],
                  'formats': [float,      float, float, float]
                  }
        self.observer2Data2006 = np.loadtxt(io.BytesIO(data), dtype=dtypes)
        self.observer2Data = self.observer2Data1931
        # Load CIE 2006 10deg standard observer data
        data = resource_string(__name__, "data/CIE_2006_STDcolorObs_10deg.txt")

        dtypes = {'names':   ['wavelen',   'x',  'y',  'z'],
                  'formats': [float,      float, float, float]
                  }
        self.observer10Data2006 = np.loadtxt(io.BytesIO(data), dtype=dtypes)
        self.observer10Data = self.observer10Data1964
        ######################################

        # Load Test Sample Spectrums
        data = resource_string(__name__, 'data/TCS1-14Spectrum.txt')

        dtypes = {'names':   ['wavelen',  'TCS01', 'TCS02',
                              'TCS03', 'TCS04', 'TCS05', 'TCS06', 'TCS07',
                              'TCS08', 'TCS09', 'TCS10', 'TCS11', 'TCS12',
                              'TCS13', 'TCS14', 'TCS15'],
                  'formats': [float, float, float, float, float,
                              float, float, float, float, float,
                              float, float, float, float, float,
                              float]
                  }
        self.tscSpect = np.loadtxt(io.BytesIO(data),
                                   dtype=dtypes, delimiter=',')

    def calculateLumFluxFromSpectrum(self, wavelen, intens, observer=2):
        '''
        Calculates luminous flux for given spectra.
        wavelen in nm
        '''
        if observer == 2:
            y = griddata(self.observer2Data['wavelen'],
                         self.observer2Data['y'], wavelen, fill_value=0.0)
        else:
            y = griddata(self.observer10Data['wavelen'],
                         self.observer10Data['y'], wavelen, fill_value=0.0)

        return 683.002 * np.trapz(y * intens, x=wavelen)

    def calculateCCT(self, xyz):
        '''
        Calculates the correlated color temperature.

        Approximation for range 3 000K - 800 000K
        '''
        xe = 0.3366
        ye = 0.1735
        A0 = -949.86315
        A1 = 6253.80338
        t1 = 0.92159
        A2 = 28.70599
        t2 = 0.20039
        A3 = 0.00004
        t3 = 0.07125

        n = (xyz[0] - xe) / (xyz[1] - ye)

        # Approximation for range 3 000 - 50 000 K
        CCT = A0 + A1 * np.exp(-n / t1) + A2 * \
            np.exp(-n / t2) + A3 * np.exp(-n / t3)

        # Approximation for range 2856K - 6504K
        if CCT < 3000:
            n = (xyz[0] - 0.3320) / (xyz[1] - 0.1858)
            CCT = -449.0 * n ** 3 + 3525.0 * n ** 2 - 6823.3 * n + 5520.33
            if CCT < 2856.0:
                warnings.warn('CCT out of approximation range')

        # Range 50 000 - 800 000 K
        elif CCT > 50000.0:
            xe = 0.3356
            ye = 0.1691
            A0 = 36284.48953
            A1 = 0.00228
            t1 = 0.07861
            A2 = 5.4535e-36
            t2 = 0.01543

            n = (xyz[0] - xe) / (xyz[1] - ye)

            CCT = A0 + A1 * np.exp(-n / t1) + A2 * \
                np.exp(-n / t2) + A3 * np.exp(-n / t3)

        if CCT > 800000.0:
            warnings.warn('CCT out of approximation range')
            CCT = np.NAN

        return CCT

    def cieXYZFromSpectrum(self, wavelen, intens,
                           reflect=None,
                           observer=2,
                           normalize=100.0):
        '''
        Will calculate the tristimulus values for given spectra.
        Y is normalized to 100 by default.

        wavelen in nm

        Use normalize = None for no normalization
        observer = 2 for CIE 1931 2 deg std observer
        observer = 10 for CIE 1964 10 deg std observer

        reflect is an array of reflectance or transmittance coefficients [0,1],
        if None, they all will be ones.
        '''
        if observer == 2:
            x = griddata(self.observer2Data['wavelen'],
                         self.observer2Data['x'], wavelen, fill_value=0.0)
            y = griddata(self.observer2Data['wavelen'],
                         self.observer2Data['y'], wavelen, fill_value=0.0)
            z = griddata(self.observer2Data['wavelen'],
                         self.observer2Data['z'], wavelen, fill_value=0.0)
        else:
            x = griddata(self.observer10Data['wavelen'],
                         self.observer10Data['x'], wavelen, fill_value=0.0)
            y = griddata(self.observer10Data['wavelen'],
                         self.observer10Data['y'], wavelen, fill_value=0.0)
            z = griddata(self.observer10Data['wavelen'],
                         self.observer10Data['z'], wavelen, fill_value=0.0)

        if reflect is None:
            reflect = np.ones(len(x), dtype=float)

        X = np.trapz(x * reflect * intens, x=wavelen)
        Y = np.trapz(y * reflect * intens, x=wavelen)
        Yn = np.trapz(y * intens, x=wavelen)
        Z = np.trapz(z * reflect * intens, x=wavelen)

        XYZ = np.array([X, Y, Z], dtype=float)

        # Normalize
        if normalize is not None:
            XYZ /= Yn
            XYZ *= normalize

        return XYZ

    def ciexyzFromXYZ(self, XYZ):
        L = np.sum(XYZ)
        return XYZ / L

    def ciexyzFromSpectrum(self, wavelen, intens, reflect=None):
        return self.ciexyzFromXYZ(self.cieXYZFromSpectrum(wavelen,
                                                          intens,
                                                          reflect=reflect
                                                          )
                                  )

    def calculateCRI(self, wavelen, intens):
        '''
        Calculates CIE 1995 CRI values
        wavelen in nm
        '''
        # Get color coordinates
        xyz = self.ciexyzFromSpectrum(wavelen, intens)
        uv = self.ct.transformToCIEuv(xyz)

        T = self.calculateCCT(xyz)

        if T > 1.0e100:
            return np.zeros(15)

        # Select reference source
        if T < 5000:
            refIntens = self.illuminants.blackbodySpectrum(wavelen, T)
        else:
            refIntens = self.illuminants.illuminantD(wavelen, T)

        ref_xyz = self.ciexyzFromSpectrum(wavelen, refIntens)
        ref_uv = self.ct.transformToCIEuv(ref_xyz)

        if np.sqrt(np.sum((uv - ref_uv) ** 2)) > 5.4e-3:
            warnings.warn('Test light not white enough! CRI has no meaning!')

        CRI = np.zeros(15)
        for tsc in np.arange(1, 15):

            sampleSp = self.testSampleSpectrum(wavelen, tsc)

            # Illuminate Test samples under reference light
            refTest_XYZ = self.cieXYZFromSpectrum(wavelen,
                                                  refIntens,
                                                  reflect=sampleSp)
            # Illuminate Test samples under test light
            test_XYZ = self.cieXYZFromSpectrum(wavelen,
                                               intens,
                                               reflect=sampleSp)

            refTest_xyz = self.ct.transformToxyz(refTest_XYZ)
            test_xyz = self.ct.transformToxyz(test_XYZ)

            refTest_uv = self.ct.transformToCIEuv(refTest_xyz)
            test_uv = self.ct.transformToCIEuv(test_xyz)

            test_uv_adapt = self.chromaticAdaptation(uv, ref_uv, test_uv)

            refTest_UVW = self.ct.transformToCIEUVW(refTest_uv,
                                                    ref_uv,
                                                    refTest_XYZ[1])
            test_UVW = self.ct.transformToCIEUVW(test_uv_adapt,
                                                 ref_uv,
                                                 test_XYZ[1])

            CRI[tsc] = 100.0 - 4.6 * \
                np.sqrt(np.sum((test_UVW - refTest_UVW) ** 2))

        # 0-100 Scaling
        CRI[1:] = 10.0 * np.log(np.exp(CRI[1:] / 10.0) + 1.0)

        CRI[0] = np.average(CRI[1:9])

        return CRI

    def testSampleSpectrum(self, wavelen, TCS=1):
        '''
        Returns the spectrum of CIE Standard test color sample 1-14

        wavelen in nm
        '''
        return griddata(self.tscSpect['wavelen'],
                        self.tscSpect['TCS{:{fill}2d}'.format(TCS, fill=0)],
                        wavelen, fill_value=0.0)

    def chromaticAdaptation(self, uvt, uvr, uvi):
        '''
        CIE (1995) uses this von Kries chromatic transform equation.

        Arguments:
        uvt: Test light chromaticity values u, v  (numpy array or list)
        uvr: Reference lamp chromaticity values u, v
        uvi: TCS test sample chromaticity values u, v

        Returns numpy array of chromaticly adapted [u,v]
        '''
        ct = (4.0 - uvt[0] - 10.0 * uvt[1]) / uvt[1]
        dt = (1.708 * uvt[1] - 1.481 * uvt[0] + 0.404) / uvt[1]
        cti = (4.0 - uvi[0] - 10.0 * uvi[1]) / uvi[1]
        dti = (1.708 * uvi[1] - 1.481 * uvi[0] + 0.404) / uvi[1]
        cr = (4.0 - uvr[0] - 10.0 * uvr[1]) / uvr[1]
        dr = (1.708 * uvr[1] - 1.481 * uvr[0] + 0.404) / uvr[1]
        return np.array([(10.872 + 0.404 * cr / ct * cti - 4.0 * dr /
                          dt * dti) /
                         (16.518 + 1.481 * cr / ct * cti - dr / dt * dti),
                         (5.520) /
                         (16.518 + 1.481 * cr / ct * cti - dr / dt * dti)
                         ])

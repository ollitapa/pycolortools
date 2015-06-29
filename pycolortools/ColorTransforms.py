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


class ColorTransforms(object):

    def __init__(self):
        '''
        Constructor
        '''
        pass

    def rgbMatrix(self, space='sRGB', refWhite='D65'):

        if (space == 'AdobeRGB' and refWhite == 'D65'):

            a = np.array([[0.5767309,    0.1855540,    0.1881852, ],
                          [0.2973769,    0.6273491,    0.0752741, ],
                          [0.0270343,    0.0706872,    0.9911085]]
                         )
        if (space == 'sRGB' and refWhite == 'D65'):
            a = np.array([[0.4124564,    0.3575761,    0.1804375],
                          [0.2126729,    0.7151522,    0.0721750],
                          [0.0193339,    0.1191920,    0.9503041]
                          ])
        if (space == 'SMPTE-C' and refWhite == 'D65'):
            a = np.array([[0.3935891,    0.3652497,    0.1916313],
                          [0.2124132,    0.7010437,    0.0865432],
                          [0.0187423,    0.1119313,    0.9581563]
                          ])

        return a

    def invrgbMatrix(self, space='sRGB', refWhite='D65'):
        a = self.rgbMatrix(space, refWhite)
        return np.linalg.inv(a)

    def sRGBgammaCorrection(self, C):
        if C > 1.0:
            C = 1.0
        if C < 0.0:
            C = 0.0

        if C > 0.0031308:
            return (1 + 0.055) * C ** (1 / 2.4) - 0.055
        else:
            return 12.92 * C

    def XYZToRGB(self, xyz, space='sRGB', refWhite='D65'):
        return np.dot(self.invrgbMatrix(space, refWhite), xyz)

    def RGBToXYZ(self, rgb, space='sRGB', refWhite='D65'):
        return np.dot(self.rgbMatrix(space, refWhite), rgb)

    def RGBToirgb(self, RGB):
        return np.array([self.sRGBgammaCorrection(RGB[0]),
                         self.sRGBgammaCorrection(RGB[1]),
                         self.sRGBgammaCorrection(RGB[2])])

    def xyzTorgb(self, xyz):
        RGB = self.XYZToRGB(xyz)
        return self.RGBToirgb(RGB)

    def transformToLab(self, XYZ, XYZn):
        def f(t):
            if t > (6.0 / 29.0) ** 3:
                return t ** (1.0 / 3.0)
            else:
                return 1.0 / 3.0 * (29.0 / 6.0) ** 2 * t + 4.0 / 29.0

        L = 116.0 * f(XYZ[1] / XYZn[1]) - 16.0
        a = 500.0 * (f(XYZ[0] / XYZn[0]) - f(XYZ[1] / XYZn[1]))
        b = 200.0 * (f(XYZ[1] / XYZn[1]) - f(XYZ[2] / XYZn[2]))

        return np.array([L, a, b])

    def transformToxyz(self, XYZ):
        L = np.sum(XYZ)
        return XYZ / L

    def transformToXYZ(self, xy, Y=100.0):
        return np.array((Y * xy[0] / xy[1], Y, Y * (1 - xy[0] - xy[1]) / xy[1]))

    def transformToRGB(self, XYZ):
        M = np.array(
            [[3.1956, 2.4478, -0.1434], [-2.5455, 7.0492, 0.9963], [0, 0, 1.0]])
        return np.dot(M, XYZ)

    def transformToCIEuv(self, xyz):
        '''
        Transforms CIE 1960 xyz colors to uvw.
        '''
        return np.array([(4.0 * xyz[0]) / (12.0 * xyz[1] - 2.0 * xyz[0] + 3),
                         (6.0 * xyz[1]) / (12.0 * xyz[1] - 2.0 * xyz[0] + 3)])

    def transformToCIEUVW(self, uv, ref_uv, Y=100.0):
        '''
        Transforms CIE 1964 xyz colors to U* V* W*.
        '''
        W = 25.0 * Y ** (1 / 3.0) - 17.0
        U = 13.0 * W * (uv[0] - ref_uv[0])
        V = 13.0 * W * (uv[1] - ref_uv[1])
        return np.array((U, V, W))

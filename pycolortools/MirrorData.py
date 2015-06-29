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
from pkg_resources import resource_string
import io


class MirrorData(object):

    def __init__(self):
        '''
        Constructor
        '''

        # Directory for datafiles
        data = resource_string(
            __name__, "data/ThorlabsMirrorPF20-03-G01_20deg.txt")

        # Load A Thorlabs mirror PF20-03-G01 regular reflectance data
        dtypes = {'names':   ['wavelen',   'p-pol',  's-pol',  'ave'],
                  'formats': [float,      float, float, float]
                  }

        self.PF20_03_G01_deg_20 = np.loadtxt(io.BytesIO(data), dtype=dtypes)

        data = resource_string(
            __name__, "data/ThorlabsMirrorPF20-03-G01_8deg.txt")

        # Load A Thorlabs mirror PF20-03-G01 regular reflectance data
        dtypes = {'names':   ['wavelen',   'p-pol',  's-pol',  'ave'],
                  'formats': [float,      float, float, float]
                  }

        self.PF20_03_G01_deg_8 = np.loadtxt(io.BytesIO(data), dtype=dtypes)

    def regReflect(self, wavelen, lens='PF20-03-G01', pol='ave', deg='8deg'):
        '''
        Returns regular reflectance data for mirrors.

        lens = 'PF20-03-G01' Thorlabs PF20-03-G01 
        pol  = p-pol, s-pol, ave
        deg  = 8deg, 20deg
        '''
        if lens == 'PF20-03-G01':

            if deg == '8deg':
                y = griddata(self.PF20_03_G01_deg_8['wavelen'],
                             self.PF20_03_G01_deg_8[pol],
                             wavelen, fill_value=0.0)
            else:
                y = griddata(self.PF20_03_G01_deg_20['wavelen'],
                             self.PF20_03_G01_deg_20[pol],
                             wavelen, fill_value=0.0)
        else:
            y = 0

        return y

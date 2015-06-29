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


class PlotColors(object):

    def __init__(self):
        '''
        Constructor
        '''
        self.colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
                       "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]
        self.extendedColors = [
            "#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78",
            "#2ca02c", "#98df8a", "#d62728", "#ff9896",
            "#9467bd", "#c5b0d5",  "#8c564b", "#c49c94",
            "#e377c2", "#f7b6d2", "#7f7f7f", "#c7c7c7",
            "#bcbd22", "#dbdb8d", "#17becf", "#9edae5"]

        self.lines = ['-', '--', '-.', ':']
        self.markers = []

    def plotColors(self, size, extended=False):
        '''
        Returns a list of nice colors for plotting.
        @param size: Number of colors desired 
        @param extended: Will extend the color list by adding lighter versions. 
        '''
        if size > len(self.colors) or extended:

            colors = self.extendedColors
            while size > len(colors):
                colors.extend(colors)

            return colors[:size]
        else:
            return self.colors[:size]

    def plotStyles(self, size, extended=False):

        if size > len(self.colors) or extended:

            ss = len(self.extendedColors)
            linegroup = [self.lines[0]] * ss
            lines = linegroup

            i = 1
            while size > len(lines):
                lines.extend([self.lines[int(i % 4)]] * ss)
                ++i

            colors = self.extendedColors
            while size > len(colors):
                colors.extend(colors)

            return colors[:size], lines[:size]
        else:
            lines = [self.lines[0]] * size
            return self.colors[:size], lines[:size]

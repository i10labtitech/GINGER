# Copyright (C) 2018 Itoh Laboratory, Tokyo Institute of Technology
# 
# This file is part of GINGER.
# 
# GINGER is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# GINGER is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with GINGER; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

import sys
import numpy as np
import matplotlib.pyplot as plt

args = sys.argv

data = np.loadtxt(args[1], delimiter="\t")
data_log = np.log10(data)
data_histo = np.histogram(data_log, bins=[i * 0.05 for i in range(120)])

#print(data_histo)

predicted_average = 0
for i in range(len(data_histo[0])):
    if i >= 4 and i < (len(data_histo[0]) - 4):
        point = data_histo[0][i]
        tmp_max = max(data_histo[0][i-4:i+5])
        if tmp_max == point and tmp_max > len(data)/100:
            print(str(data_histo[0][i]) + "\t" + "{:.2f}".format(data_histo[1][i]))
            predicted_average = data_histo[1][i] + 0.025

data_trim = data_log[data_log > predicted_average]
#print(data_trim)

diff_sum = 0
for i in data_trim:
    diff_sum += abs(i -predicted_average) ** 2
predicted_SD = (diff_sum / (len(data_trim)-1) ) ** 0.5
print("average" + "\t" + "{:.3f}".format(predicted_average))
print("SD" + "\t" + "{:.3f}".format(predicted_SD))
print("3sigma" + "\t" + "{:.3f}".format(predicted_SD*3))
print("average + 3sigma" + "\t" + "{:.3f}".format(predicted_SD*3 + predicted_average))
print("predicted limit" + "\t" + "{:.3f}".format(10**(predicted_SD*3 + predicted_average)))

limit = 10**(predicted_SD*3 + predicted_average)

data_over = data[data > limit]
print(len(data_over))

# plot
'''
fig =  plt.figure()
ax = fig.add_subplot(111)
ax.hist(data_trim, bins = [i * 0.05 for i in range(120)])
fig.show()
plt.show()
'''

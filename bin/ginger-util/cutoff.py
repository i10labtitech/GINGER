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
#import matplotlib.pyplot as plt

args = sys.argv

data = np.loadtxt(args[1], delimiter="\t")
data_histo = np.histogram(data, bins=[i * 0.05 for i in range(201)])

#print(data_histo)

zone_range = 5
z = zone_range

data_smooth = [np.mean(data_histo[0][i-z:i+z+1]) for i in range(len(data_histo[1]))[z:-z]]
data_delta = [data_smooth[i+1] - data_smooth[i] for i in range(len(data_smooth)-1)]
data_delta_x = [i+0.025 for i in data_histo[1][z:-z-1]]

#print(data_smooth)

min_data = []
max_data = []

predicted_average = 0
for i in range(len(data_smooth)):
    if i >= 1 and i < (len(data_smooth) - 1):
        point = data_smooth[i]
        tmp_min = min(data_smooth[i-1:i+2])
        if tmp_min == point and tmp_min > 10:
            #print("min: " + str(data_smooth[i]) + "\t" + "{:.2f}".format(data_histo[1][i+z]))
            min_data.append(data_histo[1][i+z])

highest = max([i for i in data_histo[0]])
right_max = 0
for i, j in zip(data_histo[0], data_histo[1]):
    if i > highest*0.5:
        right_max = j
#print(right_max)

threshold = 0
for i in min_data:
    if i > right_max:
        threshold = i
        print("{:.2f}".format(i))
        break

#print(data_trim)

#diff_sum = 0
#for i in data_trim:
#    diff_sum += abs(i -predicted_average) ** 2
#predicted_SD = (diff_sum / (len(data_trim)-1) ) ** 0.5
#print("average" + "\t" + "{:.3f}".format(predicted_average))
#print("SD" + "\t" + "{:.3f}".format(predicted_SD))
#print("3sigma" + "\t" + "{:.3f}".format(predicted_SD*3))
#print("average + 3sigma" + "\t" + "{:.3f}".format(predicted_SD*3 + predicted_average))
#print("predicted limit" + "\t" + "{:.3f}".format(10**(predicted_SD*3 + predicted_average)))

#limit = 10**(predicted_SD*3 + predicted_average)
#
#data_over = data[data > limit]
#print(len(data_over))

# plot
y_lim = 2000
'''
fig =  plt.figure()
ax = fig.add_subplot(111)
ax.hist(data, bins = [i * 0.1 for i in range(101)])
ax.vlines([threshold], 0, y_lim, colors="k", linestyle="solid")
ax.set_ylim(0,y_lim)
fig.show()
'''
#plt.show()


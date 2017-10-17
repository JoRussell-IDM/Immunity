import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import datetime
import json
from pprint import pprint
from malaria.createSimDirectoryMap import createSimDirectoryMap
from operator import itemgetter
import itertools
import scipy.stats as stats
def hex_to_RGB(hex):
  ''' "#FFFFFF" -> [255,255,255] '''
  # Pass 16 to the integer function for change of base
  return [int(hex[i:i+2], 16) for i in range(1,6,2)]


def RGB_to_hex(RGB):
  ''' [255,255,255] -> "#FFFFFF" '''
  # Components need to be integers for hex to make sense
  RGB = [int(x) for x in RGB]
  return "#"+"".join(["0{0:x}".format(v) if v < 16 else
            "{0:x}".format(v) for v in RGB])

def color_dict(gradient):
  ''' Takes in a list of RGB sub-lists and returns dictionary of
    colors in RGB and hex form for use in a graphing function
    defined later on '''
  return {"hex":[RGB_to_hex(RGB) for RGB in gradient],
      "r":[RGB[0] for RGB in gradient],
      "g":[RGB[1] for RGB in gradient],
      "b":[RGB[2] for RGB in gradient]}

def linear_gradient(start_hex, finish_hex="#FFFFFF", n=10):
  ''' returns a gradient list of (n) colors between
    two hex colors. start_hex and finish_hex
    should be the full six-digit color string,
    inlcuding the number sign ("#FFFFFF") '''
  # Starting and ending colors in RGB form
  s = hex_to_RGB(start_hex)
  f = hex_to_RGB(finish_hex)
  # Initilize a list of the output colors with the starting color
  RGB_list = [s]
  # Calcuate a color at each evenly spaced value of t from 1 to n
  for t in range(1, n):
    # Interpolate RGB vector for color at the current value of t
    curr_vector = [
      int(s[j] + (float(t)/(n-1))*(f[j]-s[j]))
      for j in range(3)
    ]
    # Add it to our list of output colors
    RGB_list.append(curr_vector)

  return color_dict(RGB_list)


exp_list = ['Outbreak_healthseek_51yr_filtered','Outbreak_nohealthseek_51yr_filtered']
list_exp_channels = []

for expname in exp_list:
    exp_dir = createSimDirectoryMap(expname)
    exp_dir.sort_values('outbreak_start_day', inplace=True)
    exp_dir.reset_index(drop= True, inplace=True)

    interval_list = exp_dir.outbreak_start_day.unique()
    exp_channels = []
    for interval,exp in exp_dir.groupby(by = 'outbreak_start_day'):
        infected_channel_by_startday = []
        for p in exp.index:
            IC_path = os.path.join(exp.outpath[p], 'output')
            files = [x for x in os.listdir(IC_path) if 'Filtered' in x]
            interval = interval
            try:
                with open(os.path.join(IC_path,files[0])) as IC:
                    data = json.load(IC)
                    infected = data['Channels']['Infected']['Data']
                    infected_channel_by_startday.append(infected)
            except:
                print('IC not available for start day' + str(interval))
        exp_channels.append(infected_channel_by_startday)
    list_exp_channels.append(exp_channels)
#colors = ['#B21F35','#FF7435','#FFCB35','#009E47','#0052A5','#06A9FC','#681E7E','#000000']
colors = linear_gradient('#0000cc','#cc0000',n=51)
fig,axarr = plt.subplots(nrows = 1, ncols = 2,sharey = True)


for i in range(len(exp_list)):#range(len(infected_channel)):
    for j in range(0,len(interval_list),10):
        for p in range(len(list_exp_channels[i][j])):
            axarr[i].plot(list_exp_channels[i][j][p],color = colors['hex'][j],alpha = 0.2)
        axarr[i].plot(np.mean(list_exp_channels[i][j],axis = 0),color = colors['hex'][j], alpha = 1)
   # _ = plt.plot(infected_channel[i][(interval_list[i]):],color = colors[i],alpha = 1) #for plotting from InsetChart


_  = plt.show()


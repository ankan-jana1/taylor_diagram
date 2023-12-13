# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 15:48:21 2022

@author: ankan_jana
"""


import statistics
import pandas as pd
import numpy as np
import logging
import matplotlib.pyplot as plt
from matplotlib.projections import PolarAxes
import mpl_toolkits.axisartist.grid_finder as gf
import mpl_toolkits.axisartist.floating_axes as fa


class TaylorDiagram(object):
  def __init__(self, STD ,fig=None, rect=111, label='_'):
    self.STD = STD
    tr = PolarAxes.PolarTransform()
    # Correlation labels
    rlocs = np.concatenate(((np.arange(11.0) / 10.0), [0.95, 0.99]))
    tlocs = np.arccos(rlocs) # Conversion to polar angles
    gl1 = gf.FixedLocator(tlocs) # Positions
    tf1 = gf.DictFormatter(dict(zip(tlocs, map(str, rlocs))))
    # Standard deviation axis extent
    self.smin = 0
    self.smax = 1.6 * self.STD
    gh = fa.GridHelperCurveLinear(tr,extremes=(0,(np.pi/2),self.smin,self.smax),grid_locator1=gl1,tick_formatter1=tf1,)
    if fig is None:
      fig = plt.figure()
    ax = fa.FloatingSubplot(fig, rect, grid_helper=gh)
    fig.add_subplot(ax)
    # Angle axis
    ax.axis['top'].set_axis_direction('bottom')
    ax.axis['top'].label.set_text("Correlation coefficient")
    ax.axis['top'].toggle(ticklabels=True, label=True)
    ax.axis['top'].major_ticklabels.set_axis_direction('top')
    ax.axis['top'].label.set_axis_direction('top')
    # X axis
    ax.axis['left'].set_axis_direction('bottom')
    ax.axis['left'].label.set_text("Standard deviation")
    ax.axis['left'].toggle(ticklabels=True, label=True)
    ax.axis['left'].major_ticklabels.set_axis_direction('bottom')
    ax.axis['left'].label.set_axis_direction('bottom')
    # Y axis
    ax.axis['right'].set_axis_direction('top')
    ax.axis['right'].label.set_text("Standard deviation")
    ax.axis['right'].toggle(ticklabels=True, label=True)
    ax.axis['right'].major_ticklabels.set_axis_direction('left')
    ax.axis['right'].label.set_axis_direction('top')
    # Useless
    ax.axis['bottom'].set_visible(False)
    # Contours along standard deviations
    ax.grid()
    self._ax = ax # Graphical axes
    self.ax = ax.get_aux_axes(tr) # Polar coordinates
    # Add reference point and STD contour
    l , = self.ax.plot([0], self.STD, 'k*', ls='', ms=8, label=label)
    t = np.linspace(0, (np.pi / 2.0))
    r = np.zeros_like(t) + self.STD
    self.ax.plot(t, r, 'k--', label='_')
    # Collect sample points for latter use (e.g. legend)
    self.samplePoints = [l]
  def add_sample(self,STD,r,*args,**kwargs):
    l,= self.ax.plot(np.arccos(r), STD, *args, **kwargs) # (theta, radius)
    self.samplePoints.append(l)
    return l
  def add_contours(self,levels=5,**kwargs):
    rs, ts = np.meshgrid(np.linspace(self.smin, self.smax), np.linspace(0, (np.pi / 2.0)))
    RMSE=np.sqrt(np.power(self.STD, 2) + np.power(rs, 2) - (2.0 * self.STD * rs  *np.cos(ts)))
    contours = self.ax.contour(ts, rs, RMSE, levels, **kwargs)
    return contours

def srl(obsSTD, s, r, l, fname):
  fig=plt.figure(figsize=(6,6))
  dia=TaylorDiagram(obsSTD, fig=fig, rect=111, label='Observed')
  plt.clabel(dia.add_contours(colors='0.5'), inline=1, fontsize=12)
  cs = plt.matplotlib.cm.Set1(np.linspace(0, 1, len(l)))
  srlc = zip(s, r, l, cs)
  for i in srlc:
    dia.add_sample(i[0], i[1], label=i[2], c=i[3], marker='s', ls='')
  spl = [p.get_label() for p in dia.samplePoints]
  fig.legend(dia.samplePoints, spl, numpoints=1, prop=dict(size='small'), loc=[0.75, 0.75])
  plt.title('Taylor diagram')
  plt.savefig(fname, dpi=300, bbox_inches='tight')
  plt.clf()
  plt.close(fig)
  
def correlationcoefficient(evaluation, simulation):
    """
    Correlation Coefficient
        .. math::
         r = \\frac{\\sum ^n _{i=1}(e_i - \\bar{e})(s_i - \\bar{s})}{\\sqrt{\\sum ^n _{i=1}(e_i - \\bar{e})^2} \\sqrt{\\sum ^n _{i=1}(s_i - \\bar{s})^2}}
    :evaluation: Observed data to compared with simulation data.
    :type: list
    :simulation: simulation data to compared with evaluation data
    :type: list
    :return: Corelation Coefficient
    :rtype: float
    """
    if len(evaluation) == len(simulation):
        correlation_coefficient = np.corrcoef(evaluation, simulation)[0, 1]
        return correlation_coefficient
    else:
        logging.warning("evaluation and simulation lists does not have the same length.")
        return np.nan


data=pd.read_excel('combinedresults.xlsx')


data_imd = list(data['IMD'])
data_sebal= list(data['SEBAL'])
data_sebs = list(data['SEBS'])
data_ssebi = list(data['SSEBI'])
#data_ssebop = list(data['SSEBop'])
data_metric = list(data['METRIC'])

sd_imd= statistics.stdev(data_imd)
sd_sebal= statistics.stdev(data_sebal)
sd_sebs= statistics.stdev(data_sebs)
sd_ssebi= statistics.stdev(data_ssebi)
#sd_ssebop= statistics.stdev(data_ssebop)
sd_metric= statistics.stdev(data_metric)

coref_sebal=correlationcoefficient(data_imd, data_sebal)
coref_sebs=correlationcoefficient(data_imd, data_sebs)
coref_ssebi=correlationcoefficient(data_imd, data_ssebi)
#coref_ssebop=correlationcoefficient(data_imd, data_ssebop)
coref_metric=correlationcoefficient(data_imd, data_metric)

obsSTD = sd_imd


s = [sd_sebal,sd_sebs,sd_ssebi,sd_metric]
#s1 = [0.980035327,0.997244197, 1.003002031]
r = [coref_sebal,coref_sebs,coref_ssebi,coref_metric]

#r1 = [0.82,0.72,0.8]

l = ['SEBAL','SEBS','S-SEBI', 'METRIC']
#l1 = ['A', 'B','C']

fname = 'TaylorDiagram.png'
srl(obsSTD, s,  r,  l, fname)

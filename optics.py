#!/usr/bin/env python

"""
To install pyclustering: 

git clone https://github.com/annoviko/pyclustering
cd pyclustering
pip install . --user
"""

import numpy as np
import random as rd
import colorlover as cl

from pyclustering.cluster.optics import (optics, ordering_analyser,
                                         ordering_visualizer)
from pyclustering.samples.definitions import FCPS_SAMPLES
from pyclustering.utils import read_sample

import plotly.offline as py
from plotly.tools import make_subplots
import plotly.graph_objs as go



def gaussian_clusters(num_clusters):
  sample = []
  for i in range(1, num_clusters + 1):
    for _ in range(40):
      sample.append([i * 1.5 + rd.random(), i * 1.5 + rd.random()])

  sample = np.array(sample)
  np.random.shuffle(sample)
  sample = sample.tolist()
  return sample



def analyze(sample, radius, neighbors):
  # Run cluster analysis where connectivity radius is bigger than real

  # Create OPTICS algorithm for cluster analysis
  optics_instance = optics(sample, radius, neighbors)

  # Run cluster analysis
  optics_instance.process()

  # Obtain results of clustering
  clusters = optics_instance.get_clusters()
  noise = optics_instance.get_noise()

  # Obtain reachability-distances
  ordering = ordering_analyser(optics_instance.get_ordering())

  return ordering, clusters, noise



def plot_results(sample, ordering, clusters, noise, name, epsilon_cutoff):
  # Find valleys boundaries
  cluster_ordering = np.array(ordering.cluster_ordering)
  if epsilon_cutoff > np.max(cluster_ordering):
    epsilon_cutoff = np.max(cluster_ordering)
  peaks = np.where(cluster_ordering > epsilon_cutoff)[0]
  val_boundaries = [0] + list(peaks) + [len(ordering.cluster_ordering)]
  print 'valleys boundaries: %s' % val_boundaries
  num_valleys = len(val_boundaries) - 1
  print 'num_valleys: %s' % num_valleys

  # Color palette (max 9 colors)
  num_colors = min(len(cl.scales), 9)
  color_palette = cl.scales[str(num_colors)]
  palette_file = 'palette.html'
  # Save and look at the color palette.
  with open(palette_file, 'w+') as f:
    f.write(cl.to_html(color_palette))
    print 'color palette saved to: %s' % palette_file
  colors = color_palette['qual']['Set1']
  # Make the set bigger if need be
  while len(colors) < num_valleys:
    colors += colors

  # Data traces to plot
  traces1 = [go.Scatter(
    name='Cluster %s data' % i,
    x=[sample[sample_id][0] for sample_id in
       clusters[0][val_boundaries[i]:val_boundaries[i + 1]]],
    y=[sample[sample_id][1] for sample_id in
       clusters[0][val_boundaries[i]:val_boundaries[i + 1]]],
    mode='markers',
    marker=dict(color=colors[i])
  ) for i in range(num_valleys)]

  traces2 = [go.Scatter(
    name='Cluster %s reachability distances' % i,
    x=range(val_boundaries[i], val_boundaries[i + 1] + 1),
    y=ordering.cluster_ordering[val_boundaries[i]:val_boundaries[i + 1] + 1],
    mode='lines',
    line=dict(color=colors[i])
  ) for i in range(num_valleys)]

  traces2.append(go.Scatter(
    name='Reachability threshold',
    y=[epsilon_cutoff for _ in ordering.cluster_ordering],
    mode='lines',
    line=dict(color='grey', dash='dash')
  ))

  # Traces layouts
  layout1 = go.Layout(title='Input data',
                      xaxis=go.XAxis(title='x'),
                      yaxis=go.YAxis(title='y'))
  layout2 = go.Layout(title='Clustering structure',
                      xaxis=go.XAxis(title='Index (cluster order of the '
                                           'objects)'),
                      yaxis=go.YAxis(title='Reachability distance'),
                      annotations=[
                        dict(
                          x=val_boundaries[i],
                          y=ordering.cluster_ordering[val_boundaries[i]],
                          xref='x',
                          yref='y',
                          text='Reachability distance over threshold',
                          showarrow=True,
                          arrowhead=7,
                          ax=0,
                          ay=-40
                        ) for i in range(1, num_valleys)]
                      )

  fig1 = {'data': traces1, 'layout': layout1}
  url1 = py.plot(fig1, auto_open=False, filename='%s_data.html' % name)
  print url1

  fig2 = go.Figure(data=traces2, layout=layout2)
  url2 = py.plot(fig2, auto_open=False, filename='%s_reachability.html' % name)
  print url2

  # With subplots:
  fig3 = make_subplots(rows=1, cols=2, subplot_titles=['Input Data',
                                                       'Clustering Structure'])
  for trace1 in traces1:
    fig3.append_trace(trace1, row=1, col=1)
  for trace2 in traces2:
    fig3.append_trace(trace2, row=1, col=2)

  # Subplot layout example:
  fig3['layout'].update(xaxis1=go.XAxis(title='x'),
                        yaxis1=go.YAxis(title='y'))
  # Or:
  fig3['layout']['xaxis1'].update(title='X values')  # , range=[0, 6])
  fig3['layout']['yaxis1'].update(title='Y values')  # , range=[0, 6])
  fig3['layout']['xaxis2'].update(title='Index (cluster order of the objects)')
  fig3['layout']['yaxis2'].update(
    title='Reachability distance')  # , range=[0,1])

  description = ('<b><a href="https://en.wikipedia.org/wiki/OPTICS_algorithm">'
                 'OPTICS</a> %s Example</b>' % name)

  fig3['layout'].update(width=1200, height=600, title=description)

  url3 = py.plot(fig3, filename='%s_clustering-structure.html' % name
                 , auto_open=False)
  print url3



def main():
  # OPTICS demo with 2D vectors
  radius = 2.0
  neighbors = 2
  epsilon_cutoff = 0.5
  sample_2d = read_sample(FCPS_SAMPLES.SAMPLE_LSUN)  # or: gaussian_clusters(3)
  ordering_2d, clusters_2d, noise_2d = analyze(sample_2d, radius, neighbors)
  plot_results(sample_2d, ordering_2d, clusters_2d, noise_2d, '2D',
               epsilon_cutoff)

  # Plot input data and clustering structure
  ordering_visualizer.show_ordering_diagram(ordering_2d)



if __name__ == '__main__':
  main()

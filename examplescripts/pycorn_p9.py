#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
PyCORN - script to extract data from files generated
by UNICORN Chromatography software supplied with ÄKTA Systems
and to plot the data with plotnine.
'''

import os
import pandas as pd
import plotnine as p9
import numpy as np
import glob
import re
import multiprocessing

# set working directory to script directory
os.chdir(os.path.dirname(os.path.realpath(__file__)))

# import sys
# basepath = os.path.dirname(os.path.abspath(__file__))
# modpath = os.path.join(basepath , "..")
# print(modpath)
# sys.path.append(modpath)
from pycorn import pc_uni6

# example for reading out start times from selected files
# for path in glob.glob('../../export/*.zip'):
#     print(path)
#     data = pc_uni6(path)
#     print(re.search('[0-9]{1,2}:[0-9]{2}:[0-9]{2} [AP]M',data['Run Log']['data'][0][1]).group())
#     filename = os.path.basename(path)[:-4]
# sys.exit()

# example for plotting data over time
# data = pc_uni6('../../export/20230119_mito_overlay_SCX_m13_FOR_frac_StepAndLongerGradient_30tubes_FTfrac3 001 001.zip')
# df = pd.DataFrame.from_dict(data.volume_to_time('UV 1_214 (1)'))
# df[['time','value']] = pd.DataFrame(df['data'].tolist(), index=df.index)
# df = df.drop(columns=['data'])
# df['time'] = pd.to_numeric(df['time'])
# df['value'] = pd.to_numeric(df['value'])
# plot = (p9.ggplot(df)
#     + p9.geom_point(p9.aes(x='time', y='value'),size=0.3)
# )
# plot.draw(show=True)
# sys.exit()

def make_limits(range):
    '''Make limits for y-axis'''
    return (-0.1*range[1], np.max(range))

def plot_data(path):
    '''Plot data from UNICORN .zip-file'''
    # read data from file
    data = pc_uni6(path)
    filename = os.path.basename(path)[:-4]

    # get chromatography type from filename (based on personal naming convention)
    chrom_type_lst = re.findall(r'SEC|SAX|SCX', path)
    if len(chrom_type_lst) > 0:
        chrom_type = chrom_type_lst[-1]
    else:
        chrom_type = 'UnkownChromatographyType'

    # set facet space and data filter pattern based on chromatography type
    if chrom_type == 'SEC':
        facet_space = {'x': [1], 'y': [1, 6]}
        data_filter_pattern = r'^UV [0-9]{1}_([1-9]{1}[0-9]{0,2})[^)M]+$|ystem flow$'
    else:
        facet_space = {'x': [1], 'y': [1, 1, 1, 4]}
        data_filter_pattern = r'^UV [0-9]{1}_([1-9]{1}[0-9]{0,2})[^)]+$|ystem flow$|^Cond$|Conc'
        
    # prepare data for plotting
    df = data.to_pandas()
    if 'data_type' not in df.columns:
        print('No data found in file {}'.format(path))
    df['data_type_unit'] = df['data_type'] + ' [' + df['unit'] + ']'

    # filter data
    df_filtered = df[df['data_name'].str.contains(data_filter_pattern)]
    df_filtered = df[df['data_name'].str.contains(data_filter_pattern)]

    # convert values to numeric
    df_filtered.loc[:, 'value'] = pd.to_numeric(df_filtered.loc[:, 'value'])

    # replace data names
    df_filtered['data_name'].replace(to_replace=r'^UV ([0-9]){1}_([1-9]{1}[0-9]{0,2})(.*)$', value=r'UV \1 (\2 nm) \3', regex=True, inplace=True)
    df_filtered['data_name'].replace(to_replace=r'^Conc B$', value=r'Concentration B', regex=True, inplace=True)
    df_filtered['data_name'].replace(to_replace=r'^Cond$', value=r'Conduction', regex=True, inplace=True)
    df_filtered['data_name'].replace(to_replace=r'^([^@]+).+SUB$', value=r'\1(blank corr.)', regex=True, inplace=True)

    # filter data for fractions
    data_fraction = df[df['data_name'].str.contains('Fraction')]
    data_fraction.drop(columns=['data_type_unit'], inplace=True) # drop column to avoid error when plotting
    data_fraction = data_fraction[~data_fraction['value'].str.contains('Frac')]
    data_fraction['data_type_unit'] = 'UV [mAU]' # plot in UV panel

    # generate plot
    plot = (p9.ggplot(df_filtered)
        + p9.geom_vline(p9.aes(xintercept='volume'), data=data_fraction, linetype='dashed', color='#777777')
        + p9.geom_point(p9.aes(x='volume', y='value', color='data_name'),size=0.1)
        + p9.facet_grid('data_type_unit ~ ', scales='free', space=facet_space)
        # + p9.scale_y_continuous(limits=make_limits)
        + p9.geom_text(p9.aes(x='volume',y=0, label='value'), data=data_fraction, size=5, colour='#777777', va='top', ha='left', angle=90)
        + p9.theme_light()
        + p9.labels.xlab('Volume [ml]')
        + p9.labels.ylab('')
        + p9.themes.theme(
            strip_text = p9.element_text(size = 4, margin = {'t':0, 'b':0, 'l':0, 'r':0}),
            legend_position = (0.5,-0.05),
            legend_title = p9.element_blank(),
            legend_key=p9.element_rect(color = "white"),
            title=p9.element_text(size=8, face="bold"),
        )
        + p9.ggtitle(filename)
    )

    # plot.draw(show=True)

    for size in [14]:
        # plot.save(filename+'_plot_{}.pdf'.format(size), width=size*16/9, height=size, units='cm', dpi=1000, limitsize=False)
        plot.save(filename+'_plot_{}.png'.format(size), width=size*16/9, height=size, units='cm', dpi=1000, limitsize=False)

    # clean up memory 
    del data
    del df
    del df_filtered
    del data_fraction
    del plot

# get list of files to plot
filelist = glob.glob('../../export/*.zip')

# if multiprocessing is not wanted, use this:
# for path in filelist:
#    plot_data(path)

# if multiprocessing is wanted, use this:
with multiprocessing.Pool(multiprocessing.cpu_count() if multiprocessing.cpu_count() < len(filelist) else len(filelist)) as pool:
    output = set(pool.map(plot_data, filelist))
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator, LogLocator)
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from matplotlib.offsetbox import AnchoredText


#general
def set_ticks(ax, xmajor, ymajor, xminor, yminor, labelsize):
    ax.xaxis.set_major_locator(MultipleLocator(xmajor))
    ax.yaxis.set_major_locator(MultipleLocator(ymajor))

    ax.xaxis.set_minor_locator(AutoMinorLocator(xminor))
    ax.yaxis.set_minor_locator(AutoMinorLocator(yminor))

    ax.tick_params(which = 'major', length = 10, direction = 'in', labelsize = labelsize)
    ax.tick_params(which = 'minor', length = 5, direction = 'in', labelsize = labelsize)

    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    return

#catastrophic failure
def cal_cf(zspec, zphoto, mask):
    
    fail_mask = np.abs((zspec[mask] - zphoto[mask]) / (1 + zspec[mask])) > 0.15

    return fail_mask

def plot_cf(ax, zspec, zphoto, bincount, zmin, zmax, mask, title = None, xlabel = None, ylabel = None):
    
    fail_mask = cal_cf(zspec, zphoto, mask)
    
    counts, xedges, yedges, im = ax.hist2d(zspec[mask], zphoto[mask], bins = [bincount, bincount], range = [[zmin, zmax], [zmin, zmax]], norm=mplcolors.LogNorm())

#     info = r'$F_{fail}$: ' + f'{sum(fail_mask) / len(fail_mask) * 100:.2f}%' + '\n'+ r'$n_{fail}/n$'+ f'= {sum(fail_mask)}'
    info = r'$F_{fail}$: ' + f'{sum(fail_mask) / len(fail_mask) * 100:.2f}%' + '\n'+ f'{sum(fail_mask)}/{len(fail_mask)}'

    anchorz = AnchoredText(info, prop = dict(size=13), frameon = True, loc = 'upper right')
    ax.add_artist(anchorz)

    if xlabel:
        ax.set_xlabel(xlabel, fontsize = 15) #'$z_{spec}$'
    if ylabel:
        ax.set_ylabel(ylabel, fontsize = 15) # '$z_{a}$'
    if title:
        ax.set_title(title, fontsize = 15)

    return im

#NMAD
def cal_nmad(zspec, zphoto, zrange, zstep, mask):
    
    med_zs = []
    for iz, zz in enumerate(zrange):
        mask_z = (zspec >= zz - zstep) & (zspec <= zz + zstep)
        mask_final = mask & mask_z
        med_z = 1.48 * np.abs((zspec[mask_final] - zphoto[mask_final] - np.median(zspec[mask_final] - zphoto[mask_final]))/ (1 + zspec[mask_final]) )
        med_zs.append(np.median(med_z))
    
    return med_zs

def plot_nmad(ax , zspec, zphoto, zrange, zstep, mask, color, label,edgecolor = 'none', marker = 'o', s = 5, title = None, xlabel = None, ylabel = None):
    
    med_zs = cal_nmad(zspec, zphoto, zrange, zstep, mask)
    
    ax.scatter(zrange, med_zs, s = s, marker = marker, color = color, edgecolor = edgecolor, label = label)
    
    if xlabel:
        ax.set_xlabel(xlabel, fontsize = 15) #'$z_{spec}$'
    if ylabel:
        ax.set_ylabel(ylabel, fontsize = 15) # '$z_{a}$'
    if title:
        ax.set_title(title, fontsize = 15)
        
#BIAS
def cal_bias(zspec, zphoto, zrange, zstep, mask):
    
    mean_zs = []
    for iz, zz in enumerate(zrange):
        
        mask_z = (zspec >= zz - zstep) & (zspec <= zz + zstep)
        mask_safe = np.abs((zspec - zphoto) / (1 + zspec) ) < 0.15
        mask_final = mask & mask_z & mask_safe
        
        med_z = ((zspec[mask_final] - zphoto[mask_final])/ (1 + zspec[mask_final]) )
        mean_zs.append(np.mean(med_z))

    return mean_zs


def plot_bias(ax , zspec, zphoto, zrange, zstep, mask, color, label, s = 5, marker = 'o',title = None, xlabel = None, ylabel = None):
 
    mean_zs = cal_bias(zspec, zphoto, zrange, zstep, mask)
    
    ax.scatter(zrange, mean_zs, s = s, color = color, label = label, marker = marker)
    
    if xlabel:
        ax.set_xlabel(xlabel, fontsize = 15) #'$z_{spec}$'
    if ylabel:
        ax.set_ylabel(ylabel, fontsize = 15) # '$z_{a}$'
    if title:
        ax.set_title(title, fontsize = 15)

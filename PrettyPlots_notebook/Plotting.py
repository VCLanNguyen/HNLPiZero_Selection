import numpy as np
import pandas as pd

import matplotlib 
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cbook as cbook
import matplotlib.patches as patches
 
from matplotlib.ticker import NullFormatter, MaxNLocator
from matplotlib.ticker import FormatStrFormatter, StrMethodFormatter, FuncFormatter

#
plt.rcParams['savefig.facecolor']='white'
plt.rcParams['axes.facecolor'] = 'white'

# set tick width
plt.rcParams['xtick.major.size'] = 8
plt.rcParams['xtick.major.width'] = 1
plt.rcParams['xtick.minor.size'] = 3
plt.rcParams['xtick.minor.width'] = 1

plt.rcParams['ytick.major.size'] = 8
plt.rcParams['ytick.major.width'] = 1
plt.rcParams['ytick.minor.size'] = 3
plt.rcParams['ytick.minor.width'] = 1

#---------------------------------------------------------------------------------------------------------------------------------#
def plot_title(
                ax,
                title, xtitle, ytitle,
                fontsize
                ):

    #PLOT TITLE, AXIS LABEL
    ax.set_title(title, fontsize = fontsize, loc='left')
    ax.set_xlabel(xtitle, fontsize = fontsize)
    ax.set_ylabel(ytitle, fontsize = fontsize)
    #ax.get_yaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
#---------------------------------------------------------------------------------------------------------------------------------#
def plot_tick(ax, fontsize):
   #ax.xaxis.set_major_locator(MaxNLocator(5))
   #ax.yaxis.set_major_locator(MaxNLocator(5))
   ax.tick_params(bottom = True, top = True, left = True, right = True)

   ax.tick_params(axis = 'x', labelsize = fontsize, direction = 'out')
   ax.tick_params(axis = 'y', labelsize = fontsize, direction = 'out')

#---------------------------------------------------------------------------------------------------------------------------------#
def plot_2dhist(
                dfx, dfy, 
                ax, fig, 
                xmin, xmax, ymin, ymax, xnbin, ynbin,
                xlimmin, xlimmax, ylimmin, ylimmax,
                title = "", xtitle = "" , ytitle = "",
                weights = None,
                iflog = False, cmin = None, cmax = None,
                color = 'plasma', ifcbar = False,
                iflegend = False, xtextloc = 0, ytextloc = 0,
                iftext = False,
                fontsize = 16
               ):
    
    if iflog:
        norm = colors.LogNorm(vmin = cmin, vmax = cmax)
    else:
        norm = None

   
    xedges = np.arange(xmin, xmax + 2*(xmax-xmin)/xnbin, (xmax-xmin)/xnbin)
    yedges = np.arange(ymin, ymax + 2*(ymax-ymin)/ynbin, (ymax-ymin)/ynbin)

    h, xbins, ybins, im = ax.hist2d(dfx, dfy, 
                                    range=[(xmin, xmax + (xmax - xmin)/xnbin), (ymin, ymax + (ymax - ymin)/ynbin)],
                                    bins = [xedges, yedges],
                                    weights = weights,
                                    norm=norm,
                                    cmap = color
                                    )
    #SET X/Y LIMT
    ax.set_xlim(xlimmin, xlimmax)
    ax.set_ylim(ylimmin, ylimmax)
        
    #COLOR BAR
    if ifcbar:
        cbar = fig.colorbar(im, ax=ax)
        cbar.ax.tick_params(labelsize=fontsize)
    
    plot_title(ax,
              title, xtitle, ytitle,
              fontsize)

    if iflegend:
        plot_2dhist_legend(dfx, dfy, ax,
                xmin, xmax, ymin, ymax,
                fontsize
                )

    if iftext:
        xbinw = (xmax-xmin)/xnbin / 2
        ybinw = (ymax-ymin)/ynbin / 2
        for i in range(len(ybins)-1):
            for j in range(len(xbins)-1):
                if h.T[i,j] > 0:
                    ax.text(xbins[j] + xbinw,
                            ybins[i] + ybinw, 
                            "{0:.0f}".format(h.T[i,j]), 
                            color="r", ha="center", va="center", fontsize = fontsize-4)
    
    plot_tick(ax, fontsize)

    return h, xbins, ybins, im

#---------------------------------------------------------------------------------------------------------------------------------#
def plot_1dhist(
                pltdf, 
                ax,
                xmin, xmax, xnbin,
                xlimmin, xlimmax,
                ifnorm = False, iforientation = False, 
                histtype = 'step', ifstacked = False,
                weights = None,
                color = 'indigo', 
                linecolor = 'indigo', linestyle = 'solid', linewidth = 0.5,
                title = "", xtitle = "", ytitle = "", fontsize = 16,
                ifysci = False, ifxsci = False,
                ifylim = False , ylimmin = None, ylimmax = None,
                ifstatbox = False, loc = 'best',
                iflabelbox = False, label = None
                ):
        
    orienation = 'vertical'    

    if iforientation:
    	orienation = 'horizontal'
      
    if iflabelbox:
        label = label

    n, bins, patches = ax.hist(
                            pltdf,
                            bins = np.arange(xmin, xmax + 2*(xmax-xmin)/xnbin, (xmax-xmin)/xnbin),
                            weights = weights,
                            density = ifnorm,
                            histtype=histtype,
                            stacked= ifstacked,
                            orientation = orienation,
                            color = color, 
                            edgecolor = linecolor,
                            linestyle = linestyle,
                            linewidth = linewidth,
                            label = label
                        )
    
    #SET X/Y LIMT
    ax.set_xlim(xlimmin, xlimmax)

    if ifylim:
    	ax.set_ylim(ylimmin, ylimmax)
    
    #SET Scientific format
    if ifysci:
        ax.set_yscale('log')

    #SET Scientific format
    if ifxsci:
        ax.set_xscale('log')

    #AXIS/PLOT TITLE
    plot_title(ax,
              title, xtitle, ytitle, fontsize)

    #LEGEND
    if ifstatbox:
        plot_1dhist_legendbox(pltdf, ax,
                       xmin, xmax, 
                       fontsize, loc,
                       ifstatbox)

    if iflabelbox:
        ax.legend(loc=loc, fontsize=fontsize-2, fancybox=True, ncol = 1)
    
    #TICK
    plot_tick(ax, fontsize)

    return n, bins, patches
#---------------------------------------------------------------------------------------------------------------------------------#
def plot_scatter(
                dfx, dfy, 
                ax,
                xlimmin, xlimmax, 
                color = 'indigo', 
                marker = '+',
                title = "", xtitle = "", ytitle = "", fontsize = 16,
                ifysci = False, ifxsci = False,
                ifylim = False , ylimmin = None, ylimmax = None,
                iflabelbox = False, label = None, loc = 'best'
                ):
        
    if iflabelbox:
        label = label

    ax.scatter(
              dfx, dfy,
              color = color,
              marker = marker,
                label = label
                        )
    
    #SET X/Y LIMT
    ax.set_xlim(xlimmin, xlimmax)

    if ifylim:
    	ax.set_ylim(ylimmin, ylimmax)
    
    #SET Scientific format
    if ifysci:
        ax.set_yscale('log')

    #SET Scientific format
    if ifxsci:
        ax.set_xscale('log')

    #AXIS/PLOT TITLE
    plot_title(ax,
              title, xtitle, ytitle, fontsize)

    #LEGEND
    if iflabelbox:
        ax.legend(loc=loc, fontsize=fontsize, fancybox=True)
    
    #TICK
    plot_tick(ax, fontsize)

#---------------------------------------------------------------------------------------------------------------------------------#
def plot_line(
                dfx, dfy, 
                ax,
                xlimmin, xlimmax, 
                color = 'indigo', 
                linestyle = '-', lw = 2,
                title = "", xtitle = "", ytitle = "", fontsize = 16,
                ifysci = False, ifxsci = False,
                ifylim = False , ylimmin = None, ylimmax = None,
                iflabelbox = False, label = None, loc = 'best'
                ):
        
    if iflabelbox:
        label = label

    ax.plot(
              dfx, dfy,
              color = color,
              linestyle = linestyle, linewidth = lw,
              label = label
                        )
    
    #SET X/Y LIMT
    ax.set_xlim(xlimmin, xlimmax)

    if ifylim:
    	ax.set_ylim(ylimmin, ylimmax)
    
    #SET Scientific format
    if ifysci:
        ax.set_yscale('log')

    #SET Scientific format
    if ifxsci:
        ax.set_xscale('log')

    #AXIS/PLOT TITLE
    plot_title(ax,
              title, xtitle, ytitle, fontsize)

    #LEGEND
    if iflabelbox:
        ax.legend(loc=loc, fontsize=fontsize, fancybox=True)
    
    #TICK
    plot_tick(ax, fontsize)

#---------------------------------------------------------------------------------------------------------------------------------#
def plot_step(
                dfx, dfy, 
                ax,
                xlimmin, xlimmax, 
                color = 'indigo', 
                title = "", xtitle = "", ytitle = "", fontsize = 16,
                ifysci = False, ifxsci = False,
                ifylim = False , ylimmin = None, ylimmax = None,
                iflabelbox = False, label = None, loc = 'best'
                ):
        
    if iflabelbox:
        label = label

    ax.step(
              dfx, dfy,
              color = color,
                label = label
                        )
    
    #SET X/Y LIMT
    ax.set_xlim(xlimmin, xlimmax)

    if ifylim:
    	ax.set_ylim(ylimmin, ylimmax)
    
    #SET Scientific format
    if ifysci:
        ax.set_yscale('log')

    #SET Scientific format
    if ifxsci:
        ax.set_xscale('log')

    #AXIS/PLOT TITLE
    plot_title(ax,
              title, xtitle, ytitle, fontsize)

    #LEGEND
    if iflabelbox:
        ax.legend(loc=loc, fontsize=fontsize, fancybox=True)
    
    #TICK
    plot_tick(ax, fontsize)

#---------------------------------------------------------------------------------------------------------------------------------#
def plot_bar(
             bins, height, 
             ax,
             width,
             xlimmin, xlimmax,
             bottom = None,
             color = 'indigo', edgecolor = 'white', linewidth = 0.05, 
             title = "", xtitle = "", ytitle = "", fontsize = 16,
             ifysci = False, ifxsci = False,
             ifylim = False , ylimmin = None, ylimmax = None,
             iflabelbox = False, label = None, loc = 'best'
             ):
        
    if iflabelbox:
        label = label

    ax.bar(
            x=bins, height=height,
            width = width,
            bottom = bottom,
            align = 'edge',
            color= color,
            edgecolor = edgecolor,
            linewidth = linewidth,
            label = label
                        )
    
    #SET X/Y LIMT
    ax.set_xlim(xlimmin, xlimmax)

    if ifylim:
    	ax.set_ylim(ylimmin, ylimmax)
    
    #SET Scientific format
    if ifysci:
        ax.set_yscale('log')

    #SET Scientific format
    if ifxsci:
        ax.set_xscale('log')

    #AXIS/PLOT TITLE
    plot_title(ax,
              title, xtitle, ytitle, fontsize)

    #LEGEND
    if iflabelbox:
        ax.legend(loc=loc, fontsize=fontsize, fancybox=True)
    
    #TICK
    plot_tick(ax, fontsize)
#---------------------------------------------------------------------------------------------------------------------------------#        
def plot_1dhist_legendbox(pltdf, ax,
                       xmin, xmax, 
                       fontsize, loc,
                       ifstatbox = False
                       ):
    # create a list with two empty handles (or more if needed)
    handles = [patches.Rectangle((0, 0), 1, 1, fc="white", ec="white", 
                                 lw=0, alpha=0)] * 3

    labels = []
    # create the corresponding number of labels (= the text you want to display)
    labels.append("Entries = {0:.4g}".format(len(pltdf)))
    labels.append("$\mu$ = {0:.3g}".format(pltdf.mean()))
    labels.append("$\sigma$ = {0:.3g}".format(pltdf.std()))

    # create the legend, supressing the blank space of the empty line symbol and the
    # padding between symbol and label by setting handlelenght and handletextpad
    ax.legend(handles, labels, loc=loc, fontsize=fontsize, 
          fancybox=True, framealpha=0.7, 
          handlelength=0, handletextpad=0)
        
#---------------------------------------------------------------------------------------------------------------------------------#        
def plot_1dhist_labelbox(pltdf, ax,
                       xmin, xmax, 
                       fontsize, loc):
    # create a list with two empty handles (or more if needed)
    handles = [patches.Rectangle((0, 0), 1, 1, fc="white", ec="white", 
                                 lw=0, alpha=0)] * 3

    # create the corresponding number of labels (= the text you want to display)
    labels = []
    labels.append("Entries = {0:.4g}".format(len(pltdf)))
    labels.append("$\mu$ = {0:.3g}".format(pltdf.mean()))
    labels.append("$\sigma$ = {0:.3g}".format(pltdf.std()))

    # create the legend, supressing the blank space of the empty line symbol and the
    # padding between symbol and label by setting handlelenght and handletextpad
    ax.legend(handles, labels, loc=loc, fontsize=fontsize, 
          fancybox=True, framealpha=0.7, 
          handlelength=0, handletextpad=0)
        
#---------------------------------------------------------------------------------------------------------------------------------#
def plot_2dhist_legend(dfx, dfy, ax,
                xmin, xmax, ymin, ymax,
                fontsize
                ):


    when = ((dfy > ymin) & (dfy < ymax) & (dfx > xmin) & (dfx < xmax))
    # create a list with two empty handles (or more if needed)
    handles = [patches.Rectangle((0, 0), 1, 1, fc="white", ec="white", 
                                 lw=0, alpha=0)] * 3

    # create the corresponding number of labels (= the text you want to display)
    labels = []
    labels.append("Entries = {0:.4g}".format(dfy.loc[when].count()))
    labels.append("$\mu$ = {0:.4g}".format(dfy.loc[when].mean()))
    labels.append("$\sigma$ = {0:.4g}".format(dfy.loc[when].std()))

    # create the legend, supressing the blank space of the empty line symbol and the
    # padding between symbol and label by setting handlelenght and handletextpad
    ax.legend(handles, labels, loc='best', fontsize=fontsize, 
          fancybox=True, framealpha=0.7, 
          handlelength=0, handletextpad=0)

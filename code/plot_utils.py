import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import numpy as np
import os
import utils as ut


def plot_ew(x, y, spectra, x_rang,
            line,
            lw=0.225, ts=1.9, ss=0.7,
            color_rang=[0,1],
            types=['O', 'B'],
            name='',
            savefold='../figures/ew/',
            filename='',
            plot_bars=True,
            what='line',
            fit=True):

    spec_filt = [
        'O4If', 'O4V', 'O7V', 'O7Ib(f)', 
        'B0V', 'B0Ib', 'B5Ia', 'B5V',
        'A0V', 'A0I', 'A5Ia', 'A5V', 
        'F0II', 'F0V', 'F5I', 'F5V', 
        'G0V', 'G0I', 'G5V', 'G4Ia',
        'K0V', 'K0Ib', 'K4p5Ib', 'K5V',
        'M0I', 'M0V', 'M4II', 'M4p5V'
        ]    
    
    spec_filt = [spec for spec in spec_filt if spec[0] in types]


    idx = np.where((x >= x_rang[0]) & (x <= x_rang[1]))[0]
    x, y = x[idx], y[idx]

    size, rat = 720, 1.15
    fig_w, fig_h = size*rat, size
    subplots = (2,1)
    fig, axs, fs, gs = initialize_figure(fig_w=fig_w, fig_h=fig_h, ratio=None, 
                                            theme=None, subplots=subplots, text_size=1,
                                            hspace=-1, hr=[1.5,1]
                                            )
    
    ax = axs[0][0]
    ax.plot(x, y, lw=lw*fs*1.5, color='k', label='abs', zorder=10)
    # x axis and labels on top
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    # change tick label size
    ax.tick_params(axis='both', which='major', labelsize=fs*ts)
    ax.set_xlabel('Wavelength (Å)', fontsize=fs*ts)
    ax.set_ylabel('Absorption Depth', fontsize=fs*ts)
    ax.text(0.1, 0.8, f'{line} {what}', fontsize=fs*ts*2, 
        transform=ax.transAxes, fontname='Arial')
    ax.set_xlim(x_rang)


    ew_dict = {}
    models = ['gaussian', 'gaussian', 'gaussian']
    k = 0

    dic_color = {}
    for i, key in enumerate(spectra.keys()):
        spec_type = key.split('_')[-1]
        color = ut.color_spectra(spec_type, cmap='turbo', range=color_rang)



        if any(substring in key for substring in spec_filt):
            if spec_type not in dic_color.keys():
                dic_color[spec_type[0]] = color

            x_, y_ = spectra[key]['x'], spectra[key]['abs']
            idx = np.where((x_ >= x_rang[0]) & (x_ <= x_rang[1]))[0]
            x_, y_ = x_[idx], y_[idx]
            ew, y_fit = ut.get_equivalent_width(x_, y_, x_rang=x_rang, return_fit=True, model='lorentz', fit=fit)
            ax.plot(x_, y_, lw=lw*fs, label=key, color=color)
            #ax.plot(x_, y_fit, color=color, ls='--', lw=lw*fs)
            ew_dict[key] = ew
            k+=1

    ax1 = axs[1][0]
    colors = [ut.color_spectra(key.split('_')[-1], cmap='turbo', range=color_rang) for key in ew_dict.keys()] 
    names = [key.split('_')[-1] for key in ew_dict.keys()] + [name]
    ew_dict['extra'] = ut.get_equivalent_width(x, y, x_rang=x_rang, return_fit=False, model='gaussian')
    colors = colors + ['k']


    sorted_items = sorted(enumerate(ew_dict.items()), key=lambda x: x[1][1])
    idxs, sorted_items = zip(*sorted_items)
    ew_dict = dict(sorted_items)
    colors = [colors[i] for i in idxs]
    names = [names[i] for i in idxs]


    if plot_bars==True:
        ax1.bar(ew_dict.keys(), ew_dict.values(), color=colors)
        ax1.set_ylabel('Equivalent Width (Å)', fontsize=fs*ts)
        ax1.set_xlabel('Spectral Type', fontsize=fs*ts, labelpad=0.25*fs*ts)
        ax1.set_xticklabels(names, rotation=45, ha='right')
        ax1.tick_params(axis='both', which='major', labelsize=fs*ts, pad=0.5*fs*ts)
        ax1.set_ylim(0, np.max(list(ew_dict.values()))*1.1)


        # plot on top of bars the value
        val_ref = ew_dict['extra']
        min_val = np.min(list(ew_dict.values()))
        max_val = np.max(list(ew_dict.values()))
        dif = max_val - min_val
        for i, (key, val) in enumerate(ew_dict.items()):
            ax1.text(i, val+dif*0.01, f'{val/val_ref:.2f}', fontsize=fs*ts*0.7, ha='center', va='bottom')

    flatten = lambda l: [item for sublist in l for item in sublist]
    for ax_ in flatten(axs):
        ax_.tick_params(axis='both', which='major', pad=0.5*fs*ts)

    colors_ = list(set(colors))+['k']
    colors_ = list(dic_color.values())[::-1] + ['k']
    lines = [mlines.Line2D([], [], color=c, marker='s', markersize=1*fs, label=t, linestyle='None', zorder=15) for t, c in zip(types+[name], colors_)]

    # Display the legend
    ax.legend(handles=lines, loc='upper right', 
              fontsize=fs*ts, framealpha=0.5, ncol=1)

    if filename != '':

        savefold = savefold
        if not os.path.exists(savefold):
            # create the folder if it does not exist
            os.makedirs(savefold)
        figname = f'{filename}.png'
        savepath = savefold + figname

        fig.savefig(savepath, dpi=300, bbox_inches='tight')


    return

def simp_figure(
    fig_size=100,
    fig_w=1080, fig_h=1080,
    text_size=1, grid=True, theme="dark",
    color='#000000',
    dpi=300,
    layout='constrained'
):
    ratio = fig_w / fig_h
    fig_width = fig_w / dpi
    fig_height = fig_h / dpi
    fig_size = fig_width * fig_height
    fs = np.sqrt(fig_size)
    fig = plt.figure(
        figsize=(fig_width, fig_height),
        dpi=dpi,  # Default dpi, will adjust later for saving
        layout=layout,
    )



    ax = fig.subplots()

    if theme == 'dark':
        fig.patch.set_facecolor(color)
        plt.rcParams.update({"text.color": "white"})
        ax.set_facecolor(color)
        ax.tick_params(colors="white")
        ax.spines["bottom"].set_color("white")
        ax.spines["top"].set_color("white")
        ax.spines["left"].set_color("white")
        ax.spines["right"].set_color("white")
        ax.xaxis.label.set_color("white")
        ax.yaxis.label.set_color("white")

    ax.xaxis.set_tick_params(which="minor", bottom=False)
    ax.tick_params(
        axis="both",
        which="major",
        labelsize=1.5 * text_size * fs,
        size=fs * 0.5,
        width=fs * 0.15,
    )
    if grid:
        ax.grid(
            which="major",
            linewidth=fs * 0.015,
            color="white" if theme == "dark" else "black",
        )
    for spine in ax.spines.values():
        spine.set_linewidth(fs * 0.15)

    # axes equal
    #ax.set_aspect("equal")
            
    return fig, ax, fs


def initialize_figure(
    fig_size=20, ratio=1,
    fig_w=512, fig_h=512,
    text_size=1, subplots=(1, 1), grid=True, theme="dark",
    color='#222222',
    dpi=300,
    wr=None, hr=None, hmerge=None, wmerge=None,
    layout='constrained',
    hspace=None, wspace=None,
    tick_direction='out',
    minor=False,
    top_bool=False
):
    """
    Initialize a Matplotlib figure with a specified size, aspect ratio, text size, and theme.

    Parameters:
    fig_size (float): The size of the figure.
    ratio (float): The aspect ratio of the figure.
    text_size (float): The base text size for the figure.
    subplots (tuple): The number of subplots, specified as a tuple (rows, cols).
    grid (bool): Whether to display a grid on the figure.
    theme (str): The theme for the figure ("dark" or any other string for a light theme).

    Returns:
    fig (matplotlib.figure.Figure): The initialized Matplotlib figure.
    ax (list): A 2D list of axes for the subplots.
    fs (float): The scaling factor for the figure size.
    """
    if ratio is not None:
        fs = np.sqrt(fig_size)
        fig = plt.figure(
            figsize=(np.sqrt(ratio * fig_size), np.sqrt(fig_size / ratio)),
            dpi=dpi,
            layout=layout,
        )
    else:
        dpi = dpi
        ratio = fig_w / fig_h
        fig_width = fig_w / dpi
        fig_height = fig_h / dpi
        fig_size = fig_width * fig_height
        fs = np.sqrt(fig_size)
        fig = plt.figure(
            figsize=(fig_width, fig_height),
            dpi=dpi,  # Default dpi, will adjust later for saving
            layout=layout,
        )

    if wr is None:
        wr_ = [1] * subplots[1]
    else:
        wr_ = wr
    if hr is None:
        hr_ = [1] * subplots[0]
    else:
        hr_ = hr
    

    gs = mpl.gridspec.GridSpec(subplots[0], subplots[1], figure=fig, width_ratios=wr_, height_ratios=hr_, hspace=hspace, wspace=wspace)


    ax = [[None] * subplots[1] for _ in range(subplots[0])]

    if theme == "dark":
        fig.patch.set_facecolor(color)
        plt.rcParams.update({"text.color": "white"})

    for i in range(subplots[0]):
        for j in range(subplots[1]):
            
            if hmerge is not None:
                if i in hmerge:
                    ax[i][j] = fig.add_subplot(gs[i, :])
                else:
                    ax[i][j] = fig.add_subplot(gs[i, j])
            elif wmerge is not None:
                if j in wmerge:
                    ax[i][j] = fig.add_subplot(gs[:, j])
                else:
                    ax[i][j] = fig.add_subplot(gs[i, j])
            else:
                ax[i][j] = fig.add_subplot(gs[i, j])

            if theme == "dark":
                ax[i][j].set_facecolor(color)
                ax[i][j].tick_params(colors="white")
                ax[i][j].spines["bottom"].set_color("white")
                ax[i][j].spines["top"].set_color("white")
                ax[i][j].spines["left"].set_color("white")
                ax[i][j].spines["right"].set_color("white")
                ax[i][j].xaxis.label.set_color("white")
                ax[i][j].yaxis.label.set_color("white")

            #ax[i][j].xaxis.set_tick_params(which="minor", bottom=False)

            if grid:
                ax[i][j].grid(
                    which="major",
                    linewidth=fs * 0.015,
                    color="white" if theme == "dark" else "black",
                )
            for spine in ax[i][j].spines.values():
                spine.set_linewidth(fs * 0.15)

            ax[i][j].tick_params(
                axis="both",
                which="major",
                labelsize=1.5 * text_size * fs,
                size=fs * 0.5,
                width=fs * 0.15,
                top=top_bool,
                direction=tick_direction
            )

            if minor:
                ax[i][j].minorticks_on()
                ax[i][j].tick_params(axis='both', which="minor", 
                direction=tick_direction,
                top=top_bool,
                size=fs * 0.25, width=fs * 0.15,)

    if hmerge is not None:
        for k in hmerge:
            for l in range(1, subplots[1]):
                fig.delaxes(ax[k][l])

    if wmerge is not None:
        for k in wmerge:
            for l in range(1, subplots[0]):
                fig.delaxes(ax[l][k])
            
    
    return fig, ax, fs, gs



def plot_lines(ax, spectral_lines, colors, xl, fs, lw=0.1, ts=1, ss=0.7,
               gap=-0.065, ym=[0.6, 1.3], mod=4):
    k = 0
    for i, line in enumerate(spectral_lines): 
        if any(substring in line for substring in colors.keys()):
            key = next(substring for substring in colors if substring in line)
            color = colors[key]
        else:
            color = 'b'  # default color if line is not in the dictionary
        vals = spectral_lines[line]
        if isinstance(vals, int) or isinstance(vals, float):
            vals = [vals]
        for val in vals:
            if (val < xl[1]) & (val > xl[0]):
                    m = np.mod(k,mod)
                    gap = gap
                    ymi, yma = ym[0], ym[1]
                    ax.plot([val, val], [ymi, yma], ls='--', lw=lw*fs*1.2, color=color, alpha=0.75, zorder=-1)
                    ax.text(val+2, yma+gap*m, line, fontsize=fs*ts*1.2,
                            ha='left', va='top', color=color, alpha=1,
                            fontname='Arial', zorder=40-i)
                    k+=1
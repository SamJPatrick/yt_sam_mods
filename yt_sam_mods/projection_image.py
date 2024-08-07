import copy
from matplotlib import pyplot, colors
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle
from matplotlib.ticker import FuncFormatter
import numpy as np
import os
import yt

#import locale
#locale.setlocale(locale.LC_ALL, 'en_US')

from grid_figure import GridFigure
from unyt import unyt_quantity

from yt_p2p.timeline import create_timeline

def powformat_OLD(number, cofmt="%.1f", end=False):
    exponent = 1
    if number > 1:
        exponent = int(np.ceil(np.log10(number)))
    else:
        exponent = int(np.floor(np.log10(number)))
    coef = cofmt % (number / 10**exponent)
    if float(coef) == 1:
        return "10$^{\\rm %d}$" % exponent
    if exponent == 0:
        return coef
    return "%s$\\times$10$^{\\rm %s}$" % (coef, exponent)

def powformat(number, end=False):
    cofmt = "%.1f"
    exponent = 1
    if number > 1:
        exponent = int(np.floor(np.log10(number)))
    else:
        exponent = int(np.floor(np.log10(number)))
    coef = cofmt % (number / 10**exponent)
    if not end:
        return "10$^{\\rm %d}$" % exponent
    return "%s$\\times$10$^{\\rm %s}$" % (coef, exponent)

def logformat(number, format="%d", end=False):
    return format % (np.log10(number))

def intcommaformat(number, end=False):
    return locale.format("%d", number, grouping=True)

def intcommaformat2(number, end=False):
    if number >= 1:
        return locale.format("%d", number, grouping=True)
    return locale.format("%.1f", number, grouping=True)
    
def intformat(number):
    return "%d" % number

def clean_image(data, threshold):

    size = 1
    my_filter = data > threshold
    x, y = np.where(my_filter)
    for i in range(x.size):
        neighbors = data[x[i] - size: x[i] + size + 1,
                         y[i] - size: y[i] + size + 1]
        cvalue = neighbors[neighbors < threshold].mean()
        data[x[i] - size: x[i] + size + 1,
             y[i] - size: y[i] + size + 1] = cvalue
        
    return data

def my_timeline(my_fig, config):

    height = config.get("height", 0.03)
    bottom = config.get("bottom_buffer", 0.08)
    h_buffer = config.get("horizontal_buffer", 0.03)

    left = h_buffer
    top = bottom + height
    width = 1.0 - 2 * h_buffer

    co = config["cosmology"]
    my_axes = my_fig.figure.add_axes((left, bottom, width, height),
                                     facecolor="black")
    create_timeline(my_axes, config["cosmology"],
                    config["current_time"], config["final_time"],
                    t_units="Myr",
                    t_major=co.arr(np.arange(0, 501, 50), "Myr"),
                    t_minor=co.arr(np.arange(0, 501, 10), "Myr"),
                    t_current=config["current_time"],
                    redshifts=np.array([100, 50, 40, 30, 25, 22,
                                        20, 18, 16, 15, 14, 14,
                                        13, 12, 11, 10]),
                    text_color="white")
    return my_axes

def single_image(panel, output_file, axes=None, fig=None,
                 figsize=(8, 8), dpi=200, fontsize=14,
                 x_range=None, x_label=None, x_tick_label_format='%d', 
                 n_x_ticks=5, x_label_position='bottom', x_label_clip=None, 
                 y_range=None, y_label=None, y_tick_label_format='%d', 
                 n_y_ticks=5, y_label_position='left', y_label_clip=None,
                 top_buffer=0.15, bottom_buffer=0.15, 
                 left_buffer=0.15, right_buffer=0.15,
                 show_cbar=True, show_cbar_label=True, 
                 cbar_length=0.9, cbar_width=0.05, cbar_buffer=0.01,
                 cbar_position='right', cbar_max_ticks=6,
                 text_color="black",
                 side="left"):
    
    if fig is None:
        fig = GridFigure(
            1, 1, figsize=figsize,
            top_buffer=top_buffer, bottom_buffer=bottom_buffer,
            left_buffer=left_buffer, right_buffer=right_buffer)

    cbar_tick_formatter = panel.get("cbar_tick_formatter")
    if cbar_tick_formatter is None:
        cbar_tick_formatter = logformat

    proj_data = None
    if "ds" in panel or ('filename' in panel and panel['filename'] is not None):
        if "ds" in panel:
            ds = panel["ds"]
        else:
            print ("Reading image data from %s." % panel['filename'])
            ds = yt.load(panel["filename"])
        proj_data = ds.data[panel['quantity']]
        if "length_bar" in panel and panel["length_bar"]:
            physical_width = ds.parameters["width"]
            current_redshift = ds.current_redshift
    elif 'data' in panel:
        proj_data = panel['data']
    else:
        proj_dim = [1.0, 1.0]

    if proj_data is not None:
        if 'scale_factor' in panel:
            proj_data *= panel['scale_factor']
        proj_dim = np.array(proj_data.shape).astype('float')
        # remove nans
        proj_data[proj_data != proj_data] = 0.0

    if panel.get("transpose"):
        proj_data = np.transpose(proj_data)
        
    if 'log_field' not in panel:
        log_field = True
    else:
        log_field = panel['log_field']

    if axes is None:
        axes = fig[0]
    panel['axes'] = axes

    if 'particle_overlay' in panel and \
            'cbar_position' in panel['particle_overlay']:
        if panel['particle_overlay']['cbar_position'] == y_label_position:
            y_label_position = None
        if panel['particle_overlay']['cbar_position'] == x_label_position:
            x_label_position = None

    if x_range is not None:
        def x_tick_format(x, pos):
            value = x * (float(x_range[1] - x_range[0]) / proj_dim[0]) + x_range[0]
            return x_tick_label_format % value
        panel['axes'].xaxis.set_major_formatter(FuncFormatter(x_tick_format))
        if n_x_ticks is not None:
            panel['axes'].xaxis.set_ticks(np.linspace(0, proj_dim[0], n_x_ticks))
        if x_label_position == 'top':
            panel['axes'].xaxis.tick_top()
        elif x_label_position == 'bottom':
             panel['axes'].xaxis.tick_bottom()
        panel['axes'].xaxis.set_ticks_position('both')
        panel["axes"].xaxis.set_tick_params(direction="in")
        tick_labels = panel['axes'].xaxis.get_ticklabels()
        for tick_label in tick_labels:
            tick_label.set_color(text_color)
            tick_label.set_size(fontsize)
            tick_label.set_color(text_color)
        if x_label_clip is not None:
            if x_label_clip == 'left' or x_label_clip == 'both':
                tick_labels[0].set_visible(False)
            if x_label_clip == 'right' or x_label_clip == 'both':
                tick_labels[-1].set_visible(False)
    if x_range is None or x_label_position is None:
        panel['axes'].xaxis.set_ticklabels([])

    if y_range is not None:
        def y_tick_format(y, pos):
            value = y * (float(y_range[1] - y_range[0]) / proj_dim[1]) + y_range[0]
            return y_tick_label_format % value
        panel['axes'].yaxis.set_major_formatter(FuncFormatter(y_tick_format))
        if n_y_ticks is not None:
            panel['axes'].yaxis.set_ticks(np.linspace(0, proj_dim[1], n_y_ticks))
        if y_label_position == 'right':
            panel['axes'].yaxis.tick_right()
        elif y_label_position == 'left':
            panel['axes'].yaxis.tick_left()
        panel['axes'].yaxis.set_ticks_position('both')
        panel["axes"].yaxis.set_tick_params(direction="in")
        tick_labels = panel['axes'].yaxis.get_ticklabels()
        for tick_label in tick_labels:
            tick_label.set_color(text_color)
            tick_label.set_size(fontsize)
            tick_label.set_color(text_color)
        if y_label_clip is not None:
            if y_label_clip == 'bottom' or y_label_clip == 'both':
                tick_labels[0].set_visible(False)
            if y_label_clip == 'top' or y_label_clip == 'both':
                tick_labels[-1].set_visible(False)
    if y_range is None or y_label_position is None:
        panel['axes'].yaxis.set_ticklabels([])

    if x_label is not None and x_label_position is not None: 
        panel['axes'].xaxis.set_label_position(x_label_position)
        panel['axes'].xaxis.set_label_text(x_label, fontsize=fontsize,
                                           color=text_color)
    if y_label is not None and y_label_position is not None:
        panel['axes'].yaxis.set_label_position(y_label_position)
        panel['axes'].yaxis.set_label_text(y_label, fontsize=fontsize,
                                           color=text_color)

    if x_range is None and y_range is None:
        panel['axes'].xaxis.set_visible(False)
        panel['axes'].yaxis.set_visible(False)

    if proj_data is not None:
        if panel['range'] is None:
            panel['range'] = ['min', 'max']
        if panel['range'][0] == 'min':
            if log_field:
                if (proj_data > 0).any():
                    image_min = proj_data[proj_data > 0].min()
                else:
                    image_min = None
            else:
                image_min = proj_data.min()
            panel['range'][0] = float(image_min)
        else:
            image_min = panel['range'][0]
            panel['range'][0] = image_min
        if panel['range'][1] == 'max':
            image_max = proj_data.max()
            panel['range'][1] = float(image_max)
        else:
            image_max = panel['range'][1]
        if panel["range"][0] > panel["range"][1]:
            panel["range"][1] = panel["range"][0] * 1.001
            image_max = panel["range"][1]

        if "floor" in panel:
            panel["range"][1] = max(panel["range"][1], panel["floor"])

        if "ceiling" in panel:
            panel["range"][1] = min(panel["ceiling"], panel["range"][1])
            image_max = panel["range"][1]

        panel["image_max"] = float(image_max)

        if log_field == 'double':
            if panel['negative_range'] is None:
                panel['negative_range'] = ['min', 'max']
            if panel['negative_range'][0] == 'min':
                if (proj_data < 0).any():
                    negative_image_min = proj_data[proj_data < 0].min()
                else:
                    negative_image_min = None
                panel['negative_range'][0] = negative_image_min
            else:
                negative_image_min = panel['negative_range'][0]
            if panel['negative_range'][1] == 'max':
                if (proj_data < 0).any():
                    negative_image_max = proj_data[proj_data < 0].max()
                else:
                    negative_image_max = None
                panel['negative_range'][1] = negative_image_max
            else:
                negative_image_max = panel['negative_range'][1]

        if image_min is None or \
            (log_field == 'double' and \
             negative_image_min is None):
            proj_data = None

    if proj_data is not None:
        if panel.get("clean"):
            print ("Cleaning:")
            proj_data = clean_image(proj_data, panel["clean_threshold"])
        
        if log_field:
            if (proj_data > 0).any():
                print ("%s: field range: %+e - %+e, image range: %+e - %+e." % \
                    (panel['quantity'],
                     proj_data[proj_data > 0].min(), proj_data[proj_data > 0].max(),
                     image_min, image_max))
                if log_field == 'double' and (proj_data < 0).any():
                    print ("%s(negative): field range: %+e - %+e, image range: %+e - %+e." % \
                        (' ' * (len(panel['quantity'])-10),
                         proj_data[proj_data < 0].min(), proj_data[proj_data < 0].max(),
                         negative_image_min, negative_image_max))
        else:
            print ("%s: field range: %e - %e, image range: %e - %e." % \
                (panel['quantity'], proj_data.min(), proj_data.max(),
                 image_min, image_max))

        if log_field == 'double':
            negative_proj_data = negative_image_min / proj_data

        if log_field == 'double':
            clip = False
        elif 'clip' in panel:
            clip = panel['clip']
        else:
            clip = True

        if clip: proj_data = proj_data.clip(image_min, image_max)

        if log_field:
            mynorm = colors.LogNorm(image_min, image_max)
        else:
            mynorm = colors.normalize(image_min, image_max)

        proj_data = np.flipud(proj_data)
        image_dims = proj_data.shape

        panel['image'] = panel['axes'].imshow(proj_data, cmap=panel['cmap'], 
                                              norm=mynorm, aspect='equal',
                                              interpolation='nearest', 
                                              extent=(0, proj_dim[0],
                                                      0, proj_dim[1]))#, origin='lower')

        if log_field == 'double':
            mynorm = colors.LogNorm(1.0, negative_image_min / negative_image_max)
            panel['negative_image'] = panel['axes'].imshow(negative_proj_data,
                                                           cmap=panel['negative_cmap'],
                                                           norm=mynorm, aspect='equal',
                                                           interpolation='nearest', 
                                                           extent=(0, proj_dim[0],
                                                                   0, proj_dim[1]))

        if show_cbar:
            if cbar_position == 'left' or cbar_position == 'right':
                cbar_orientation = 'vertical'
            elif cbar_position == 'top' or cbar_position == 'bottom':
                cbar_orientation = 'horizontal'
            else:
                print ("ERROR: cbar_position must be 'top', 'bottom', 'left', or 'right'.")
                return

            panel["cax"] = fig.add_cax(
                panel["axes"], cbar_position,
                buffer=cbar_buffer, length=cbar_length,
                width=cbar_width)

            if log_field:
                if log_field == 'double':
                    cbar_max_ticks = int(cbar_max_ticks / 2)
                tick_min = np.ceil(np.log10(panel['range'][0]))
                tick_max = np.floor(np.log10(panel['range'][1])) + 1
                tick_step = max(1, np.ceil((tick_max - tick_min) /
                                           np.float(cbar_max_ticks)))
                cbar_ticks = 10**np.arange(tick_min, tick_max, tick_step)
                if cbar_ticks.size == 0:
                    cbar_ticks = np.array(panel["range"])
                else:
                    if not np.allclose(cbar_ticks[0], panel["range"][0], atol=0):
                        cbar_ticks = np.concatenate([[panel["range"][0]], cbar_ticks])
                    if not np.allclose(cbar_ticks[-1], panel["range"][1], atol=0):
                        cbar_ticks = np.concatenate([cbar_ticks, [panel["range"][1]]])
                cbar_ticks = np.concatenate([cbar_ticks, [panel["range"][1]]])
            else:
                cbar_ticks = None
            panel['cbar'] = pyplot.colorbar(panel['image'], orientation=cbar_orientation, 
                                            cax=panel['cax'], ticks=cbar_ticks)
            panel['cbar'].solids.set_edgecolor("face")

            if log_field == 'double':
                tick_min = np.ceil(np.log10(-panel['negative_range'][1]))
                tick_max = np.floor(np.log10(-panel['negative_range'][0])) + 1
                tick_step = max(1, np.ceil((tick_max - tick_min) /
                                           np.float(cbar_max_ticks)))
                negative_cbar_ticks = 10**np.arange(tick_min, tick_max, tick_step)

                if (negative_image_min / negative_image_max) < 10.0:
                    negative_cbar_ticks = np.concatenate([[negative_image_min],
                                                          negative_cbar_ticks,
                                                          [negative_image_max]])
                negative_cbar_tick_values = -negative_image_min / negative_cbar_ticks
                panel['negative_cax'] = pyplot.axes([negative_cax_left, negative_cax_bottom, 
                                                     cbar_horizontal, cbar_vertical])
                panel['negative_cbar'] = pyplot.colorbar(panel['negative_image'],
                                                         orientation=cbar_orientation, 
                                                         cax=panel['negative_cax'],
                                                         ticks=negative_cbar_tick_values)

            if cbar_orientation == 'vertical':
                if cbar_position == 'left':
                    panel['cax'].yaxis.set_label_position('left')
                    panel['cax'].yaxis.tick_left()
                    if log_field == 'double':
                        panel['negative_cax'].yaxis.set_label_position('left')
                        panel['negative_cax'].yaxis.tick_left()
                else:
                    panel['cax'].yaxis.set_label_position('right')
                    panel['cax'].yaxis.tick_right()
                    if log_field == 'double':
                        panel['cax'].yaxis.set_label_position('right')
                        panel['cax'].yaxis.tick_right()
                panel["cax"].tick_params(axis='y', colors="black")
                ticklabels = panel['cbar'].ax.get_yticklabels()
                for ticklabel in ticklabels:
                    ticklabel.set_color(text_color)
                    ticklabel.set_size(0.75 * fontsize)
                    ticklabel.set_verticalalignment("center")
                    ticklabel.set_horizontalalignment(side)
                # ticklabels[0].set_size(0.75 * fontsize)
                # ticklabels[0].set_verticalalignment("top")
                # ticklabels[-2].set_size(0.75 * fontsize)
                # ticklabels[-2].set_verticalalignment("top")
                # ticklabels[-1].set_alpha(0.0)
                # ticklabels[-2].set_alpha(0.0)
                # ticklabels[0].set_alpha(0.0)
                
            else:
                if cbar_position == 'top':
                    panel['cax'].xaxis.set_label_position('top')
                    panel['cax'].xaxis.tick_top()
                    if log_field == 'double':
                        panel['negative_cax'].xaxis.set_label_position('top')
                        panel['negative_cax'].xaxis.tick_top()
                else:
                    panel['cax'].xaxis.set_label_position('bottom')
                    panel['cax'].xaxis.tick_bottom()
                    if log_field == 'double':
                         panel['negative_cax'].xaxis.set_label_position('bottom')
                         panel['negative_cax'].xaxis.tick_bottom()
                for ticklabel in panel['cbar'].ax.get_xticklabels(): 
                    ticklabel.set_color(text_color)
                    ticklabel.set_size(fontsize)
                ticks = panel["cbar"].ax.get_xticklabels()
                for tick in ticks:
                    tick.set_color(text_color)
                if log_field == 'double':
                    for ticklabel in panel['negative_cbar'].ax.get_xticklabels(): 
                        ticklabel.set_color(text_color)
                        ticklabel.set_size(fontsize)
                if cbar_ticks is not None:
                    panel['cbar'].ax.set_xticklabels(map(cbar_tick_formatter, 
                                                         cbar_ticks))
                    if log_field == 'double':
                        panel['negative_cbar'].ax.set_xticklabels(
                            map(cbar_tick_formatter, negative_cbar_ticks))

            if log_field == 'double':
                if 'positive_label' in panel:
                    panel['cbar'].set_label(panel['positive_label'],
                                            fontsize=0.75*fontsize)
                else:
                    panel['cbar'].set_label('(pos)', fontsize=0.75*fontsize)
                if 'negative_label' in panel:
                    panel['negative_cbar'].set_label(panel['negative_label'],
                                                     fontsize=0.75*fontsize)
                else:
                    panel['negative_cbar'].set_label('(neg)', fontsize=0.75*fontsize)
            elif show_cbar_label and 'label' in panel:
                if cbar_position in ["left", "right"]:
                    panel["cax"].yaxis.set_label_text(panel["label"],
                                                      fontsize=fontsize,
                                                      color=text_color)
                    if cbar_position == "right":
                        panel["cax"].yaxis.set_label_coords(4.2, 0.5)
                    elif cbar_position == "left":
                        panel["cax"].yaxis.set_label_coords(-3.2, 0.5)

                elif cbar_position in ["top", "bottom"]:
                    panel["cax"].xaxis.set_label_text(panel["label"],
                                                      fontsize=fontsize,
                                                      color=text_color)
                else:
                    raise ValueError(f"Can't do {cbar_position=}.")

    if "length_bar" in panel and panel["length_bar"]:
        if "length_bar_scale" in panel:
            print_length = panel["length_bar_scale"]
        else:
            print_length = 0.4 * physical_width

        if "length_bar_color" in panel:
            length_bar_color = panel["length_bar_color"]
        else:
            length_bar_color = "white"

        if "length_bar_units" in panel:
            if panel["length_bar_units"] == "auto":
                print_units = "kpc"
                if print_length.in_units(print_units) < 1.:
                    print_units = "pc"
                if print_length.in_units(print_units) < 1.:
                    print_units = "AU"
            else:
                print_units = panel["length_bar_units"]
        else:
            print_units = str(print_length.units)

        physical_width.convert_to_units(print_units)

        print_length_s = "%d %s" % \
          (print_length.in_units(print_units), print_units)
        print_units2 = "kpc"
        print_length2 = print_length / (1 + current_redshift)
        if print_length2.in_units(print_units2) < 1.:
            print_units2 = "pc"
        if print_length2.in_units(print_units2) < 1.:
            print_units2 = "AU"
        if print_units2 == "kpc":
            print_length_s2 = "%.1f %s'" % \
              (print_length2.in_units(print_units2), print_units2)
        else:
            print_length_s2 = "%d %s'" % \
              (print_length2.in_units(print_units2), print_units2)
            
        line_length = (print_length / physical_width).in_units("") * proj_data.shape[0]
        length_bar_vert = 725
        
        if "length_bar_left" in panel:
            bar_left = panel["length_bar_left"]
        elif "length_bar_right" in panel:
            bar_left = panel["length_bar_right"] - line_length
        else:
            bar_left = 50

        panel["axes"].plot([bar_left, bar_left + line_length],
                           [length_bar_vert, length_bar_vert], 
                           color=length_bar_color)
        panel["axes"].text(bar_left + 0.5 * line_length, length_bar_vert + 25, 
                           print_length_s, color=length_bar_color,
                           horizontalalignment="center", fontsize=12)
        panel["axes"].text(bar_left + 0.5 * line_length, length_bar_vert - 55, 
                           print_length_s2, color=length_bar_color,
                           horizontalalignment="center", fontsize=12)
        panel["axes"].set_xlim(0, proj_data.shape[0])
        panel["axes"].set_ylim(0, proj_data.shape[1])
                
    if 'particle_overlay' in panel and panel['particle_overlay'] is not None:
        particle_overlay(panel['particle_overlay'], panel['axes'], proj_dim, 
                         left_side, bottom_side, panel_width, fontsize=fontsize, 
                         cbar_length=cbar_length, cbar_width=cbar_width, 
                         cbar_tick_formatter=cbar_tick_formatter,
                         cbar_position=cbar_position)

    if 'text' in panel:
        for text in panel['text']:
            panel['axes'].text(text['x'], text['y'], text['s'], **text['kwargs'])

    if "ds" in panel:
        del ds
    del proj_data

    if output_file is not None:
        pyplot.savefig(output_file, bbox_inches='tight', dpi=dpi)

def multi_image(panels, output_file, n_columns=2, figsize=(8, None),
                dpi=200, fontsize=14,
                tick_range=None, tick_label=None,
                tick_label_format='%d', n_ticks=None,
                top_buffer=0.15, bottom_buffer=0.15,
                left_buffer=0.15, right_buffer=0.15,
                vertical_buffer=0, horizontal_buffer=0,
                cbar_orientation='vertical', fig_text=None,
                bg_color="white", text_color="black", timeline=None,
                **kwargs):

    n_rows = np.int(np.ceil(len(panels) / n_columns))

    my_fig = GridFigure(
        n_rows, n_columns, figsize=figsize, square=True,
        top_buffer=top_buffer, bottom_buffer=bottom_buffer,
        left_buffer=left_buffer, right_buffer=right_buffer,
        vertical_buffer=vertical_buffer,
        horizontal_buffer=horizontal_buffer)

    if timeline is not None:
        tl_axes = my_timeline(my_fig, timeline)

    x_label_clip = None
    y_label_clip = None
    if n_rows > 1: y_label_clip = 'both'
    if n_columns > 1: x_label_clip = 'both'

    for i, panel in enumerate(panels):
        if panel is None:
            continue

        my_row, my_column, _ = my_fig._get_index(i)
        my_axes = my_fig[i]

        show_cbar = True
        show_cbar_label = True
        cbar_position = None
        if panel.get("hide_cbar"):
            show_cbar = False
        else:
            if cbar_orientation == 'vertical':
                # show_cbar_label = True
                y_label = None
                y_label_position = None
                if my_column == 0:
                    cbar_position = 'left'
                elif my_column == n_columns - 1:
                    cbar_position = 'right'
                else:
                    show_cbar = False

                if my_row == n_rows - 1:
                    x_label = tick_label
                    x_label_position = 'bottom'
                elif my_row == 0:
                    x_label = tick_label
                    x_label_position = 'top'
                else:
                    x_label = None
                    x_label_position = None

            elif cbar_orientation == 'horizontal':
                x_label = None
                x_label_position = None
                if my_row == n_rows - 1:
                    cbar_position = 'bottom'
                elif my_row == 0:
                    cbar_position = 'top'
                else:
                    show_cbar = False

                if my_column == 0:
                    y_label = tick_label
                    y_label_position = 'left'
                elif my_column == n_columns - 1:
                    y_label = tick_label
                    y_label_position = 'right'
                else:
                    y_label = None
                    y_label_position = None
            else:
                show_cbar = False

        if my_column == 0:
            side = "right"
        else:
            side = "left"

        single_image(panel, None, fig=my_fig, axes=my_axes, fontsize=fontsize,
                     x_range=tick_range, x_label=x_label, 
                     x_label_position=x_label_position, 
                     n_x_ticks=n_ticks, 
                     x_tick_label_format=tick_label_format,
                     x_label_clip=x_label_clip,
                     y_range=tick_range, y_label=y_label, 
                     y_label_position=y_label_position, 
                     n_y_ticks=n_ticks, 
                     y_tick_label_format=tick_label_format,
                     y_label_clip=y_label_clip,
                     show_cbar=show_cbar, show_cbar_label=show_cbar_label,
                     cbar_position=cbar_position, text_color=text_color,
                     side=side,
                     **kwargs)

    if fig_text is not None:
        for f_text in fig_text:
            my_fig.text(f_text['x'], f_text['y'], f_text['s'], **f_text['kwargs'])

    if output_file is not None:
        my_fig.figure.savefig(output_file, dpi=dpi,
                              facecolor=bg_color, edgecolor='none')
        for panel in panels:
            if panel is not None:
                if 'image' in panel: del panel['image']
                if 'cbar' in panel: del panel['cbar']
                if 'cax' in panel: del panel['cax']
                if 'axes' in panel: del panel['axes']
    return my_fig

def particle_overlay(pdict, my_axes, proj_dim, 
                     left_side, bottom_side, panel_width, 
                     fontsize=14, show_cbar=True,
                     cbar_length=0.9, cbar_width=0.05, cbar_buffer=0.01,
                     cbar_position='left',
                     cbar_tick_formatter=logformat,
                     text_color="black"):

    cbar_tick_formatter = intformat
    p_color = None
    p_size = None
    color_filter = None
    if 'alpha' in pdict:
        alpha = pdict['alpha']
    else:
        alpha = 1.0
        
    x_pos = pdict['px']
    y_pos = pdict['py']
    p_color = pdict['pz']
    if "z_range" not in pdict:
        pdict["z_range"] = None
    if pdict['z_range'] is None:
        pdict['z_range'] = [p_color.min(), p_color.max()]
    print ("color field range: %e - %e, image range: %e - %e." % \
        (p_color.min(), p_color.max(),
         pdict['z_range'][0], pdict['z_range'][1]))
        
    p_size = pdict['pr']

    patches = []
    for i in range(x_pos.size):
        circle = pyplot.Circle((x_pos[i], y_pos[i]), p_size[i], 
                               alpha=alpha)
        patches.append(circle)
    p = PatchCollection(patches, cmap=pdict["cmap"], alpha=alpha)
    p.set_array(np.log10(p_color.d))
    if "facecolor" in pdict:
        p.set_facecolor(pdict["facecolor"])
    if "linewidth" in pdict:
        p.set_linewidth(pdict["linewidth"])
    my_axes.add_collection(p)
    my_axes.set_xlim(0, proj_dim[0])
    my_axes.set_ylim(0, proj_dim[1])

    if pdict['cbar_position'] is None:
        return
    
    if not 'cbar_position' in pdict:
        my_cbar_position = cbar_position
    else:
        my_cbar_position = pdict['cbar_position']
    if pdict['z_range'] is not None and show_cbar:
        if my_cbar_position == 'left' or my_cbar_position == 'right':
            cbar_orientation = 'vertical'
        elif my_cbar_position == 'top' or my_cbar_position == 'bottom':
            cbar_orientation = 'horizontal'
        else:
            print ("ERROR: cbar_position must be 'top', 'bottom', 'left', or 'right'.")
            return

        if cbar_orientation == 'vertical':
            if my_cbar_position == 'left':
                cax_left = left_side - (cbar_width * panel_width) - cbar_buffer
            else:
                cax_left = left_side + panel_width + cbar_buffer
            cax_bottom = bottom_side + (0.5 * panel_width * (1.0 - cbar_length))
            cbar_horizontal = cbar_width * panel_width
            cbar_vertical = cbar_length * panel_width
        else:
            if my_cbar_position == 'top':
                cax_bottom = bottom_side + panel_width + cbar_buffer
            else:
                cax_bottom = bottom_side - (panel_width * cbar_width + cbar_buffer)
            cax_left = left_side + (0.5 * panel_width * (1.0 - cbar_length))
            cbar_horizontal = cbar_length * panel_width
            cbar_vertical = cbar_width * panel_width

        my_cax = pyplot.axes([cax_left, cax_bottom, 
                              cbar_horizontal, cbar_vertical])
        cbar_ticks = np.arange(np.ceil(np.log10(pdict['z_range'][0])), 
                               (np.floor(np.log10(pdict['z_range'][1])) + 1), 1)
        my_cbar = pyplot.colorbar(p, orientation=cbar_orientation, 
                                        cax=my_cax, ticks=cbar_ticks)
        my_cbar.solids.set_edgecolor("face")

        if cbar_orientation == 'vertical':
            if my_cbar_position == 'left':
                my_cax.yaxis.set_label_position('left')
                my_cax.yaxis.tick_left()
            else:
                my_cax.yaxis.set_label_position('right')
                my_cax.yaxis.tick_right()
            for ticklabel in my_cbar.ax.get_yticklabels(): 
                ticklabel.set_color(text_color)
                ticklabel.set_size(fontsize)
            my_cbar.ax.set_yticklabels(map(cbar_tick_formatter, cbar_ticks))
        else:
            if my_cbar_position == 'top':
                my_cax.xaxis.set_label_position('top')
                my_cax.xaxis.tick_top()
            else:
                my_cax.xaxis.set_label_position('bottom')
                my_cax.xaxis.tick_bottom()
            for ticklabel in my_cbar.ax.get_xticklabels(): 
                ticklabel.set_color(text_color)
                ticklabel.set_size(fontsize)
            my_cbar.ax.set_xticklabels(map(cbar_tick_formatter, cbar_ticks))

        if 'cbar_label' in pdict:
            my_cbar.set_label(pdict['cbar_label'], fontsize=fontsize)

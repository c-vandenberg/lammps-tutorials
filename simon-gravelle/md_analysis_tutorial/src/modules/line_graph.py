import numpy
import matplotlib.pyplot as pyplot
from numpy import ndarray
from typing import List, Union, Tuple


class LineGraph:
    @staticmethod
    def line_graph_subplots(data_arrays: List[ndarray], subplot_titles: List, x_labels: List, y_labels: List,
                            y_lims: Union[Tuple, List[Tuple]], x_lim: Union[Tuple, List[Tuple]], figure_title: str):
        subplot_number: int = len(data_arrays)

        line_graphs, subplot_axes = pyplot.subplots(subplot_number, 1, figsize=(10, 3 * subplot_number))

        # If only one subplot, convert `axes` to list
        if subplot_number == 1:
            subplot_axes: List = [subplot_axes]

        line_graphs.patch.set_facecolor('black')

        # Iterate through all subplots, stylise and plot data
        for key, axes in enumerate(subplot_axes):
            # Stylise subplot
            axes.set_facecolor('black')
            axes.tick_params(colors='white', which='both')
            for spine in axes.spines.values():
                spine.set_edgecolor('white')

            # Plot data
            axes.plot(data_arrays[key][:, 0], data_arrays[key][:, 1], color='cyan')
            axes.set_title(subplot_titles[key], color='white')
            axes.set_xlabel(x_labels[key], color='white')
            axes.set_ylabel(y_labels[key], color='white')
            axes.set_ylim(y_lims[key])
            axes.set_xlim(x_lim)

        pyplot.subplots_adjust(hspace=0.5)

        line_graphs.text(0.5, 0.0005, figure_title, ha='center', va='center', color='white', fontsize=12)

        pyplot.tight_layout()
        pyplot.show()

    @staticmethod
    def single_line_graph(data_arrays: List[ndarray], line_labels: List, line_colours: List, x_label: str,
                          y_label: str, y_lim: Tuple, x_lim: Tuple, graph_title: str, figure_text: str):
        # Set up the plot
        fig, ax = pyplot.subplots(figsize=(10, 6))
        fig.patch.set_facecolor('black')
        ax.set_facecolor('black')

        # Customize ticks and spines
        ax.tick_params(colors='white', which='both')
        for spine in ax.spines.values():
            spine.set_edgecolor('white')

        # Extract data and plot
        for key, data in enumerate(data_arrays):
            y_axis_data = data[0]
            x_axis_data = data[1]
            ax.plot(y_axis_data, x_axis_data, color=line_colours[key], label=line_labels[key])

        # Set axes labels, limits and graph title
        ax.set_xlabel(x_label, color='white')
        ax.set_ylabel(y_label, color='white')
        ax.set_ylim(y_lim)
        ax.set_xlim(x_lim)

        # Add legend
        legend = ax.legend(loc='upper right', frameon=False, fontsize=12)
        for text in legend.get_texts():
            text.set_color('white')

        # Add figure title
        fig.text(0.5, 0.0005,
                 figure_text,
                 ha='center', va='center', color='white', fontsize=12)

        # Show plot
        pyplot.tight_layout()
        pyplot.show()

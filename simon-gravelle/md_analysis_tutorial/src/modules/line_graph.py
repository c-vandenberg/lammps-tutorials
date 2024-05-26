import numpy
import matplotlib.pyplot as pyplot
from numpy import ndarray
from typing import List, Union, Tuple


class LineGraph:
    @staticmethod
    def line_graph_subplots(data_arrays: List[ndarray], subplot_titles: List, x_labels: List, y_labels: List,
                            y_lims: Union[Tuple, List[Tuple]], x_lim: Union[Tuple, List[Tuple]], figure_title: str):
        subplot_number = len(data_arrays)

        fig, subplot_axes = pyplot.subplots(subplot_number, 1, figsize=(10, 3 * subplot_number))

        # If there's only one subplot, `axes` won't be a list, so we need to make it a list.
        if subplot_number == 1:
            subplot_axes = [subplot_axes]

        fig.patch.set_facecolor('black')

        for i, axes in enumerate(subplot_axes):
            axes.set_facecolor('black')
            axes.tick_params(colors='white', which='both')
            for spine in axes.spines.values():
                spine.set_edgecolor('white')

            # Plotting the data
            axes.plot(data_arrays[i][:, 0], data_arrays[i][:, 1], color='cyan')
            axes.set_title(subplot_titles[i], color='white')
            axes.set_xlabel(x_labels[i], color='white')
            axes.set_ylabel(y_labels[i], color='white')
            axes.set_ylim(y_lims[i])
            axes.set_xlim(x_lim)

        pyplot.subplots_adjust(hspace=0.5)

        fig.text(0.5, 0.0005, figure_title, ha='center', va='center', color='white', fontsize=12)

        pyplot.tight_layout()
        pyplot.show()

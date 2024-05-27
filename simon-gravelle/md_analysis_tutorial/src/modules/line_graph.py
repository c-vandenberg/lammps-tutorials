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
    def plot_bond_length_distributions(starting_bond_length_distribution,
                                       maximum_deformation_bond_length_distribution):

        # Extract data for plotting
        starting_bin_centers = starting_bond_length_distribution[0]
        starting_histogram = starting_bond_length_distribution[1]

        deformation_bin_centers = maximum_deformation_bond_length_distribution[0]
        deformation_histogram = maximum_deformation_bond_length_distribution[1]

        # Set up the plot
        fig, ax = pyplot.subplots(figsize=(10, 6))
        fig.patch.set_facecolor('black')
        ax.set_facecolor('black')

        # Customize ticks and spines
        ax.tick_params(colors='white', which='both')
        for spine in ax.spines.values():
            spine.set_edgecolor('white')

        # Plot the data
        ax.plot(starting_bin_centers, starting_histogram, color='cyan', label='At Start (Frames 1 - 20)')
        ax.plot(deformation_bin_centers, deformation_histogram, color='orange',
                label='During Maximum Deformation (Frames 200 - 220)')

        # Set labels and title
        ax.set_xlabel('Bond Length (Å)', color='white')
        ax.set_ylabel('Probability', color='white')
        ax.set_ylim([0.00, 0.13])
        ax.set_xlim([1.30, 1.65])

        # Add legend
        legend = ax.legend(loc='upper right', frameon=False, fontsize=12)
        for text in legend.get_texts():
            text.set_color('white')

        # Add figure title
        fig.text(0.5, 0.0005,
                 r'$\bf{Fig\ 2}$ Bond length distribution carbon nanotube (CNT) at start of simulation & at maximum deformation.',
                 ha='center', va='center', color='white', fontsize=12)

        # Show plot
        pyplot.tight_layout()
        pyplot.show()

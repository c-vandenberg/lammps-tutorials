import matplotlib.pyplot as pyplot
from matplotlib.legend import Legend
from numpy import ndarray
from typing import List, Union, Tuple


class LineGraph:
    @staticmethod
    def line_graph_subplots(data_arrays: List[ndarray], subplot_titles: List, x_labels: List, y_labels: List,
                            x_lim: Union[Tuple, List[Tuple]], y_lims: Union[Tuple, List[Tuple]], graph_title: str,
                            figure_text: str):
        """
        Create multiple line graph subplots.

        Parameters
        ----------
        data_arrays : List[ndarray]
            List of NumPy data arrays where each array contains x and y data for a subplot.
        subplot_titles : List[str]
            List of titles for each subplot.
        x_labels : List[str]
            List of x-axis labels for each subplot.
        y_labels : List[str]
            List of y-axis labels for each subplot.
        x_lim : Union[Tuple[float, float], List[Tuple[float, float]]]
            Limits for the x-axis. Can be a single tuple for a single subplot or a list of tuples for multiple subplots.
        y_lims : Union[Tuple[float, float], List[Tuple[float, float]]]
            Limits for the y-axis. Can be a single tuple for a single subplot or a list of tuples for multiple subplots.
        graph_title : str
            Title for the graph.
        figure_text : str
            Additional text to display at the bottom of the figure.

        Returns
        -------
        None
        """
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

        line_graphs.suptitle(graph_title, color='white')
        line_graphs.text(0.5, 0.0005, figure_text, ha='center', va='center', color='white', fontsize=12)

        pyplot.tight_layout()
        pyplot.show()

    @staticmethod
    def single_line_graph(data_arrays: List[ndarray], figure_size: Tuple[int, int],
                          line_colours: List, x_label: str, y_label: str, y_lim: Tuple, x_lim: Tuple, graph_title: str,
                          figure_text: str, figure_text_font_size: int, font_size: int, label_size: int, line_width: int,
                          line_labels: Union[List, None] = None):
        """
        Create a single line graph with multiple lines.

        Parameters
        ----------
        data_arrays : List[ndarray]
            List of NumPy data arrays where each array contains x and y data for a line plot.
        figure_size : Tuple[int, int]
            Size of the figure (width, height)
        line_colours : List[str]
            List of colors for each line.
        x_label : str
            Label for the x-axis.
        y_label : str
            Label for the y-axis.
        y_lim : Tuple[float, float]
            Limits for the y-axis.
        x_lim : Tuple[float, float]
            Limits for the x-axis.
        graph_title : str
            Title for the graph.
        figure_text : str
            Additional text to display at the bottom of the figure.
        figure_text_font_size: int
            Font size for the text displayed at the bottom of the figure.
        font_size : int
            Font size for the axis labels, tick labels and title.
        label_size : int
            Font size for the tick labels.
        line_width : int
            Width of the lines in the plot.
        line_labels : Union[List, None]
            List of labels for each line (optional).

        Returns
        -------
        None
        """
        # Set up the line graph
        line_graph, line_graph_axes = pyplot.subplots(figsize=figure_size)
        line_graph.patch.set_facecolor('black')
        line_graph_axes.set_facecolor('black')

        # Stylise ticks and spines
        line_graph_axes.tick_params(colors='white', which='both', labelsize=label_size)
        for spine in line_graph_axes.spines.values():
            spine.set_edgecolor('white')

        # Extract data and plot
        for key, data in enumerate(data_arrays):
            x_axis_data = data[0]
            y_axis_data = data[1]
            if line_labels and len(line_labels) > 0:
                line_graph_axes.plot(x_axis_data, y_axis_data, color=line_colours[key], linewidth=line_width,
                                     label=line_labels[key])
            else:
                line_graph_axes.plot(x_axis_data, y_axis_data, color=line_colours[key], linewidth=line_width)

        # Set axes labels, limits and graph title
        line_graph_axes.set_xlabel(x_label, fontsize=font_size, color='white')
        line_graph_axes.set_ylabel(y_label, fontsize=font_size, color='white')
        line_graph_axes.set_ylim(y_lim)
        line_graph_axes.set_xlim(x_lim)
        line_graph_axes.set_title(graph_title, fontsize=font_size, color='white')

        # Add legend
        legend: Legend = line_graph_axes.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False,
                                                fontsize=font_size)
        for text in legend.get_texts():
            text.set_color('white')

        # Add figure title
        line_graph.text(0.5, 0.0005,
                        figure_text,
                        ha='center', va='center', color='white', fontsize=figure_text_font_size)

        # Show plot
        pyplot.tight_layout()
        pyplot.show()

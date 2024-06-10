import matplotlib.pyplot as pyplot
from matplotlib.legend import Legend
from numpy import ndarray
from typing import List, Union, Tuple


class LineGraph:
    @staticmethod
    def single_line_graph(data_arrays: List[ndarray], figure_size: Tuple[int, int],
                          line_colours: List, x_label: str, y_label: str, x_lim: Tuple, y_lim: Tuple, graph_title: str,
                          figure_text: str, figure_text_font_size: Union[int, float],
                          figure_text_x_coord: Union[int, float], figure_text_y_coord: Union[int, float],
                          font_size: Union[int, float], tick_label_size: Union[int, float],
                          line_width: Union[int, float], save_path: Union[str, None] = None,
                          line_labels: Union[List, None] = None,
                          line_labels_position: str = None,
                          dashed_lines: List[Tuple[str, Union[int, float]]] = None):
        """
        Create a single line graph with single or multiple lines.

        Parameters
        ----------
        data_arrays : List[ndarray]
            List of NumPy data arrays where each array contains x and y data.
        figure_size : Tuple[int, int]
            Size of the figure (width, height).
        line_colours : List[str]
            List of colors for each line.
        x_label : str
            Label for the x-axis.
        y_label : str
            Label for the y-axis.
        x_lim : Tuple[float, float]
            Limits for the x-axis.
        y_lim : Tuple[float, float]
            Limits for the y-axis.
        graph_title : str
            Title for the graph.
        figure_text : str
            Additional text to display at the bottom of the figure.
        figure_text_font_size: int
            Font size for the text displayed at the bottom of the figure.
        figure_text_x_coord : Union[int, float]
            Figure text x-coordinate.
        figure_text_y_coord : Union[int, float]
            Figure text y-coordinate.
        font_size : Union[int, float]
            Font size for the axis labels, tick labels and title.
        tick_label_size : Union[int, float]
            Font size for the tick labels.
        line_width : Union[int, float]
            Width of the lines in the plot.
        save_path : Union[str, None]
            Directory path where the plot will be saved (optional).
        line_labels : Union[List, None]
            List of labels for each line (optional).
        line_labels_position : str
            Position for the line labels (optional).
        dashed_lines : List[Tuple[str, Union[int, float]]]
            Dashed lines for the x-axis and/or y-axis (optional)

        Returns
        -------
        None
        """
        # Set up the line graph
        line_graph, line_graph_axes = pyplot.subplots(figsize=figure_size)
        line_graph.patch.set_facecolor('black')
        line_graph_axes.set_facecolor('black')

        # Stylise ticks and spines
        line_graph_axes.tick_params(colors='white', which='both', labelsize=tick_label_size)
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

        if dashed_lines:
            for axis, coord in dashed_lines:
                if axis.lower() not in ['x', 'y']:
                    raise ValueError('Dashed line axis must be "x" or "y"')

                if axis.lower() == 'x':
                    line_graph_axes.axvline(coord, color='white', linestyle='--', linewidth=1.5)
                else:
                    line_graph_axes.axhline(coord, color='white', linestyle='--', linewidth=1.5)

        if line_labels:
            label_position = line_labels_position if line_labels_position is not None else 'upper right'

            # Add legend
            legend: Legend = line_graph_axes.legend(loc=label_position, frameon=False,
                                                    fontsize=font_size)
            for text in legend.get_texts():
                text.set_color('white')

        # Add figure text
        line_graph.text(figure_text_x_coord, figure_text_y_coord, figure_text, ha='center', va='center',
                        color='white', fontsize=figure_text_font_size)

        if save_path:
            pyplot.savefig(save_path, bbox_inches='tight')

        pyplot.show()

import matplotlib.pyplot as pyplot
import mpl_toolkits.mplot3d as mplot3d
import mpl_toolkits.mplot3d.art3d as art3d


class ScatterPlot:
    @staticmethod
    def two_dimensional_scatter_plot(x_axis_values: list, y_axis_values: list, x_axis_label: str, y_axis_label: str,
                                     title: str):
        # Create graph
        two_dimension_figure, two_dimension_axes = pyplot.subplots(figsize=(10, 6))

        # Stylise & plot the x and y values
        two_dimension_figure.patch.set_facecolor('black')
        two_dimension_axes.set_facecolor('black')
        two_dimension_axes.scatter(x_axis_values, y_axis_values, c='cyan', marker='o', alpha=0.6)
        two_dimension_axes.set_ylim(12, two_dimension_axes.get_ylim()[1])

        pyplot.title(title, color='white')
        pyplot.xlabel(x_axis_label, color='white')
        pyplot.ylabel(y_axis_label, color='white')
        two_dimension_axes.tick_params(colors='white', which='both')

        for spine in two_dimension_axes.spines.values():
            spine.set_edgecolor('white')

        pyplot.show()

    @staticmethod
    def three_dimensional_scatter_plot(x_axis_values: list, y_axis_values: list, z_axis_values: list,
                                       marker_colours: list, x_axis_label: str, y_axis_label: str, z_axis_label: str,
                                       marker_colours_label: str, title: str):
        three_dimension_figure: pyplot.Figure = pyplot.figure(figsize=(10, 6))
        three_dimension_axes: mplot3d.Axes3D = three_dimension_figure.add_subplot(111, projection='3d')

        three_dimension_scatter: art3d.Path3DCollection = three_dimension_axes.scatter(
            x_axis_values,
            y_axis_values,
            z_axis_values,
            c=marker_colours,
            cmap='viridis_r',
            marker='o'
        )

        # Add color bar to show frame mapping
        colour_bar: pyplot.colorbar = pyplot.colorbar(three_dimension_scatter, ax=three_dimension_axes, shrink=0.5)
        colour_bar.set_label(marker_colours_label)

        # Set titles and labels
        three_dimension_axes.set_title(title)
        three_dimension_axes.set_xlabel(x_axis_label)
        three_dimension_axes.set_ylabel(y_axis_label)
        three_dimension_axes.set_zlabel(z_axis_label)

        # Show the plot
        pyplot.show()

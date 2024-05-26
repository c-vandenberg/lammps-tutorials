import matplotlib.pyplot as pyplot


class ScatterPlot:
    @staticmethod
    def two_dimensional_scatter_plot(
            x_axis_values: list,
            y_axis_values: list,
            x_axis_label: str,
            y_axis_label: str,
            title: str
    ):
        # Create graph
        two_d_figure, two_d_axes = pyplot.subplots(figsize=(10, 6))

        # Stylise & plot the x and y values
        two_d_figure.patch.set_facecolor('black')
        two_d_axes.set_facecolor('black')
        two_d_axes.scatter(x_axis_values, y_axis_values, c='cyan', marker='o', alpha=0.6)
        two_d_axes.set_ylim(12, two_d_axes.get_ylim()[1])

        pyplot.title(title, color='white')
        pyplot.xlabel(x_axis_label, color='white')
        pyplot.ylabel(y_axis_label, color='white')
        two_d_axes.tick_params(colors='white', which='both')

        for spine in two_d_axes.spines.values():
            spine.set_edgecolor('white')

        pyplot.show()

    @staticmethod
    def three_dimensional_scatter_plot(
            x_axis_values: list,
            y_axis_values: list,
            z_axis_values: list,
            marker_colours: list,
            x_axis_label: str,
            y_axis_label: str,
            z_axis_label: str,
            marker_colours_label: str,
            title: str
    ):
        three_d_figure: pyplot.Figure = pyplot.figure(figsize=(10, 6))
        three_d_axes: pyplot.Axes = three_d_figure.add_subplot(111, projection='3d')

        three_dimension_scatter = three_d_axes.scatter(
            x_axis_values,
            y_axis_values,
            z_axis_values,
            c=marker_colours,
            cmap='viridis_r',
            marker='o'
        )

        # Add color bar to show frame mapping
        cb = pyplot.colorbar(three_dimension_scatter, ax=three_d_axes, shrink=0.5)
        cb.set_label(marker_colours_label)

        # Set titles and labels
        three_d_axes.set_title(title)
        three_d_axes.set_xlabel(x_axis_label)
        three_d_axes.set_ylabel(y_axis_label)
        three_d_axes.set_zlabel(z_axis_label)

        # Show the plot
        pyplot.show()

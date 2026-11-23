""" Functions for applying plot options (annotate, log, etc.) to plotly figures """


def format_plot_options(fig, plot_options):
    """
    Applies the plot options to a plotly figure, as format_plot_options (plot_formatting.py)
    does to matplotlib axes, calling the function of each option.

    Parameters
    ----------
    fig : plotly.graph_objects.Figure
        Figure of the plot.
    plot_options : list
        Options to configure plot.
    """

    # log axes; values <= 0 are not drawn, as in matplotlib, which applies them always
    # (log_validity drops those values before checking)
    for log_ax in ["logx", "logy"]:
        if log_ax in plot_options:
            log_axes(fig, log_ax)


def log_axes(fig, log_ax, undo=False):
    """
    Log or un-log an axis of a plotly figure, as log_axes (plot_options_matplotlib.py)
    does to a matplotlib axis.

    Parameters
    ----------
    fig : plotly.graph_objects.Figure
        Figure of the plot.
    log_ax : str
        Which axis to log. Options are 'logx' or 'logy'.
    undo : bool, optional
        If True, sets the axis scale back to linear. Default is False.
    """

    axis_type = "linear" if undo else "log"
    if log_ax == "logx":
        fig.update_xaxes(type=axis_type)
    if log_ax == "logy":
        fig.update_yaxes(type=axis_type)

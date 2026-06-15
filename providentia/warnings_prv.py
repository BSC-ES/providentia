""" Warning function """


def show_message(read_instance, msg, from_conf=None, deactivate=False, print=False):
    """
    Display or log warning messages and direct to the correct output.

    Parameters
    ----------
    read_instance : object
        Stores the instance of the current mode being used, such as
        'dashboard', 'download', 'report' or 'interpolation'.
    msg : str
        Warning message to be displayed or logged.
    from_conf : bool, optional
        Indicates whether the message originates from a configuration file.
        If True or None, the message may be shown in the dashboard.
    deactivate : bool, optional
        If True, suppresses user-facing warnings.
    print : bool, optional
        If True, forces the message to be logged regardless of execution mode.
    """

    # variable used to control when the warnings don't need to be shown
    if deactivate:
        return

    if (
        read_instance.mode in ["report", "library", "interpolation", "download"]
        or print
    ):
        read_instance.logger.warning("Warning: " + msg)

    else:
        # there are some warnings that will only be shown if we launch the dashboard
        # using a configuration file (those in filter.py, read.py and configuration.py)
        if (from_conf is None) or (from_conf is True):
            if not read_instance.delay:
                from .dashboard_elements import MessageBox

                MessageBox(msg)
            else:
                read_instance.delayed_warnings.append(msg)

""" Warning functions """


def show_message(read_instance, msg, msg_offline=None, from_conf=None, deactivate=False):
    # variable used to control when the warnings don't need to be shown
    if deactivate:
        return

    if (read_instance.offline) or (read_instance.interactive) or (read_instance.interpolation) or (read_instance.download):
        if msg_offline is not None:
            print('Warning: ' + msg_offline)
        else:
            print('Warning: ' + msg)
    
    else:
        # there are some warnings that will only be shown if we launch the dashboard
        # using a configuration file (those in filter.py, read.py and configuration.py)
        if (from_conf is None) or (from_conf):
            from dashboard_elements import MessageBox
            MessageBox(msg)

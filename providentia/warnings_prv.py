""" Warning functions """
def show_message(read_instance, msg, from_conf=None, deactivate=False, print=False):
    # variable used to control when the warnings don't need to be shown
    if deactivate:
        return

    if (read_instance.report) or (read_instance.library) or (read_instance.interpolation) or (read_instance.download) or print:
       read_instance.logger.warning('Warning: ' + msg)
    
    else:
        # there are some warnings that will only be shown if we launch the dashboard
        # using a configuration file (those in filter.py, read.py and configuration.py)
        if (from_conf is None) or (from_conf is True):
            if not read_instance.delay:
                from dashboard_elements import MessageBox
                MessageBox(msg)
            else:
                read_instance.delayed_warnings.append(msg)
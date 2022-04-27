import datetime
import pandas as pd

start_date = '20120101'
end_date = '20120220'
resolution = 'monthly'

if (resolution == 'hourly') or (resolution == 'hourly_instantaneous'):
    active_frequency_code = 'H'
elif (resolution == '3hourly') or (resolution == '3hourly_instantaneous'):
    active_frequency_code = '3H'
elif (resolution == '6hourly') or (resolution == '6hourly_instantaneous'):
    active_frequency_code = '6H'
elif resolution == 'daily':
    active_frequency_code = 'D'
elif resolution == 'monthly':
    active_frequency_code = 'MS'

str_active_start_date = str(start_date)
str_active_end_date = str(end_date)
time_array = pd.date_range(start=datetime.datetime(int(str_active_start_date[:4]),
                                                   int(str_active_start_date[4:6]),
                                                   int(str_active_start_date[6:8])),
                           end=datetime.datetime(int(str_active_end_date[:4]),
                                                 int(str_active_end_date[4:6]),
                                                 int(str_active_end_date[6:8])),
                           freq=active_frequency_code)[:-1]

print(time_array)
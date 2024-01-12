# load Providentia interactive class
from providentia import Interactive

# read and filter data from .conf
inst = Interactive(conf='test_interactive.conf')

# save data to netCDF
inst.save(format='nc')

# return data in memory (in netCDF format)
data = inst.get_data(format='nc')
print(data)

# make a Providentia plot 
ax = inst.make_plot('timeseries')
print(ax.get_xydata())
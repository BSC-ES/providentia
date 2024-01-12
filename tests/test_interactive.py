# load Providentia interactive class
from providentia import Interactive

# read and filter data from .conf
inst = Interactive(conf='test_local.conf')

# save data to netCDF
inst.save(format='nc')

# return data in memory (in netCDF format)
data = inst.get_data(format='nc')

# make a Providentia plot 
ax = inst.make_plot('timeseries')
print('Test is finished')
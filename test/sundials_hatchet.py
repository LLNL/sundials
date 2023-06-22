import os
import sys
import platform
import datetime as dt
from IPython.display import HTML, display

import hatchet as ht

# find way to load in a folder

# Setup desired cali query.
grouping_attribute = "function"
default_metric = "sum(sum#time.duration),inclusive_sum(sum#time.duration)"
query = "select function,%s group by %s format json-split" % (
    default_metric,
    grouping_attribute,
)

# cali_file_path = "/usr/workspace/pan13/sundials/adiak-testing/Benchmarking/arkode_diffusion_2D_mpi_d2d_arkode_serial/arkode_diffusion_2D_mpi_d2d_arkode_serial.06122023_151704.cali"

cali_file_path = "/usr/workspace/pan13/sundials/adiak-testing/Benchmarking/advection_reaction_3D_arkdirk_newton/advection_reaction_3D_arkdirk_newton.06092023_143842.cali"

gf = ht.GraphFrame.from_caliper(cali_file_path, query)

print(gf.dataframe)

# use the cali-query to aggregate files into a new file 

# array-to-store-graphframes
# for file in caliper folder (can you do this?)
#     create graph frames with cali files

# sort the graphframes by date
# aggregate all but the last one (most recent one)

# compare results
#     - if past a certain threshold, send an email to some list
#         - find a package for this  
#         - how acquire emails, can one pull from a list, or send to a mailing list? (like archwhisky)

# do some other analysis cause why not

# to do
# set up .cali file loads -> likely only do for the latest job (which you can get by checking the environment ID. How would you do locally? Well this wasn't really meant for personal use. one can probably make a second version of this script)
# set up function to aggregate example performance data - how? probably by checking the ci job value in the cali files
# this should be run whenever the benchmark script is run (like after the whole benchmarking and example profiling)
# compare to cached data -> if anomaly in the worst way, send out an email

# send email for general performance

# set up a way to store cached performance data. Most likely a CSV or json of some sort with CI job id, date, and relevant data. Or could just calculate straight up
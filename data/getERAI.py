#!/usr/bin/env python

# A simple script to download ERA-interim reanalysis data
# Before running you must:
#   have both ecmwfapi and cdo bindings installed
#       sudo port install py-pip (and follow instructions about port --select)
#       sudo pip install cdo
#       sudo pip install ecmwf-api-client
#   have generated .ecmwfapirc according to https://software.ecmwf.int/wiki/display/WEBAPI/Access+ECMWF+Public+Datasets#AccessECMWFPublicDatasets-key (this page as a lot more info too)
#
# The script assumes we're doing this year by year, so you can select firstYear and lastYear
#
# NB! ECMWF's public interface is unstable and the retrieval of files will
# periodically hang. I have added time-outs and retries to this script to
# address this, but you may still have to restart the script if that mechanism
# fails. Before you restart the script log on to ECMWF's web site, go to
# http://apps.ecmwf.int/webmars/joblist/ and delete all the jobs in your job
# list. This prevents new jobs from being put on hold because your job list is
# full.

# Imports
from ecmwfapi import ECMWFDataServer
from cdo import *
import datetime, calendar
import os, signal, tempfile, atexit

# Loop parameters
firstYear = 2010
lastYear  = 2015

# Request parameters
# The analysis
param_an = "151.128/165.128/166.128/167.128/168.128"
step_an  = "0"
time_an  = "00:00:00/06:00:00/12:00:00/18:00:00"

# The forecast
param_fc = "144.128/169.128/175.128/228.128"
step_fc  = "6/12"
time_fc  = "00:00:00/12:00:00"

# Constants - this is always the same for both forecast and analysis (as long as its ERAI)
exp_class = "ei"
dataset = "interim"
expver = "1"
grid = "0.5/0.5"
area = "90/-180/40/179.5"
levtype = "sfc"
stream = "oper"
fformat = "netcdf"

# Timeout duration (s)
timeout = 15*60
# Number of retries
retries = 3

# Module instances
cdo = Cdo()
server = ECMWFDataServer()

# Handler for the inevitable time out
def handler(signum, frame):
    raise Exception(__file__ + ": Connection timed out.")

signal.signal(signal.SIGALRM, handler)

# Temporary output files
tmp_fc = tempfile.mkstemp(suffix='.nc')
os.close(tmp_fc[0])
tmp_an = tempfile.mkstemp(suffix='.nc')
os.close(tmp_an[0])

# Handler to clean up - also in the case of keyboard interrupt
def cleaner(file1, file2):
    os.remove(file1)
    os.remove(file2)
    print "Gracefull exit with temporary files removed."

atexit.register(cleaner, tmp_fc[1], tmp_an[1])

########################################################################
# Loop over years and months
########################################################################

for year in range(firstYear,lastYear+1):
    for month in range(1,13):

        date1 = datetime.date(year,month,1)
        eom   = calendar.monthrange(year, month)[1]
        date2 = datetime.date(year,month,eom)

        # Check if the output file exists and skip if it does
        # Usefull if we had to restart the script for some reason
        filename = "erai.6h."  + date1.strftime('%Y%m') + ".nc"
        if os.path.isfile(filename):
            print __file__ + ": " + filename + " exists. Skipping."
            continue

        # We first pick up the analysis
        datestr = date1.strftime('%Y-%m-%d') + '/to/' + date2.strftime('%Y-%m-%d')
        for attempt in range(0,retries):
            signal.alarm(timeout)
            try:
                server.retrieve({
                    "class": exp_class,
                    "dataset": dataset,
                    "date": datestr,
                    "expver": expver,
                    "grid": grid,
                    "area": area,
                    "levtype": levtype,
                    "param": param_an,
                    "step": step_an,
                    "stream": stream,
                    "time": time_an,
                    "type": "an",
                    "format": fformat,
                    "target": tmp_an[1],
                })
                break
            except Exception, exc:
                print exc
                if attempt+1 < retries:
                    print __file__ + ": Retrying. Attempt " + str(attempt+1)
        else:
            continue

        # And then the forecast. For this we need to go one day back to get midnight on the first day of the month
        date1_fc = date1 - datetime.timedelta(1)
        datestr  = date1_fc.strftime('%Y-%m-%d') + '/to/' + date2.strftime('%Y-%m-%d')
        for attempt in range(0,retries):
            signal.alarm(timeout)
            try:
                server.retrieve({
                    "class": exp_class,
                    "dataset": dataset,
                    "date": datestr,
                    "expver": expver,
                    "grid": grid,
                    "area": area,
                    "levtype": levtype,
                    "param": param_fc,
                    "step": step_fc,
                    "stream": stream,
                    "time": time_fc,
                    "type": "fc",
                    "format": fformat,
                    "target": tmp_fc[1],
                })
                break
            except Exception, exc:
                print exc
                if attempt+1 < retries:
                    print __file__ + ": Retrying. Attempt " + str(attempt+1)
        else:
            continue

        # Use cdo to select the right date range and merge
        tmp_fc2 = cdo.seldate(date1.strftime('%Y-%m-%d'), date2.strftime('%Y-%m-%d'), input = tmp_fc[1])
        cdo.merge(input = tmp_fc2 + " " + tmp_an[1], output = filename)


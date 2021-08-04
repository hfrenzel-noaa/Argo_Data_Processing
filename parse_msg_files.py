#!/usr/bin/env python
#
# This script parses at least one msg file from an Argo float (core or BGC) for its
# engineering data and creates one netcdf output file per float
# (e.g., for use with ERDDAP).
#
# H. Frenzel, CICOES, University of Washington // NOAA-PMEL
#
# Current version: August 3, 2021
#
# First version: April 27, 2021

import argparse
import datetime
import netCDF4 as nc
import numpy as np
import os
import re
from datetime import date

import pdb

def parse_input_args():
    '''Parse the command line arguments and return them as an object.'''
    parser = argparse.ArgumentParser()

    # required argument:
    parser.add_argument('filename_in', nargs='+', help='name of the input file')
    # options:
    parser.add_argument('-d', '--directory', default='.', type=str,
                        help='working directory (default: cwd)')
    parser.add_argument('-o', '--output_directory', default='.', type=str,
                        help='output directory (default: cwd)')
    parser.add_argument('-v', '--verbose', default=False, action="store_true",
                        help='if set, display more progress updates')

    args = parser.parse_args()
    # note that input and output files are always list objects,
    # even if there is just one file
    return args


def change_cwd(dir):
    '''Change to the specified directory, unless it is None.'''
    if dir:
        os.chdir(dir)


def get_filename_out(filename_in, output_dir):
    '''Determine and return the corresponding name of the netcdf output file for
    the given name of the input file. Issue a warning message and return None
    if the input filename does not conform to the expected naming standard.'''
    regex = re.compile(r'(\d+)\.\d+.msg')
    match_obj = regex.search(filename_in)
    if match_obj:
        return '{0:s}/{1:s}.nc'.format(output_dir, match_obj.group(1))
    else:
        print('File "{0:s"} has an unexpected name'.format(filename_in))
        return None


def parse_msg_file(filename):
    '''Parse the given msg file and return a dictionary with 
    engineering variables as the keys and its values and units as
    the values. Also return a dictionary with time, lon, lat values.
    Return "None, None" if the file doesn't exist or if its size is 0.'''
    if not os.path.exists(filename) or not os.path.getsize(filename):
        return None, None, None
    manu, program = detect_float_type(filename)
    if not manu:
        return None, None, None
    # the 'GPS fix obtained' set of 3 lines appears in different places
    # for different kinds of floats, so parse the whole file
    coords = parse_msg_gps_fix(filename)
    vars = dict()
    # DEBUG print('parsing {0:s}'.format(filename)) # DEBUG
    fp = open(filename)
    line = fp.readline()
    if line.startswith('$'):
        # DEBUG print('has $ header')
        parse_msg_header(fp, vars)
     # DEBUGelse:
         # DEBUGprint('no $ header') # DEBUG
    parse_msg_footer(fp, vars)
    fp.close()
    return (vars, coords, program)


def detect_float_type(filename):
    '''Determine if the given msg file belongs to an APEX or Navis float.
    Return manufacturer ('apex' or 'navis') and program ('core' or 'bgc').
    Assumption: A line with "Apf" or "Npf" is found in the first few lines
    (usually in the first or fourth line of an msg file).'''
    # set defaults first
    program = 'core' # only changed if IsusInit or DuraInit is found in file
    manu = None
    with open(filename) as fp:        
        for line in fp:
            if 'Apf' in line and 'FwRev' in line:
                manu = 'apex'
            elif 'Npf' in line and 'FwRev' in line:
                manu = 'navis'
            elif 'IsusInit' in line or 'DuraInit' in line:
                program = 'bgc'
                if manu:
                    break
            elif 'DebugBits' in line and manu:
                break
    if not manu:
        print('Float type not detected in "{0:s}"'.format(filename))
    return manu, program
        

def scan_file_pattern(fp, pattern, pattern2=None):
    '''Look for the given pattern in the open file with the given file
    pointer fp. Return True and the line if found or False, None if not.'''
    for line in fp:
        if pattern in line and (not pattern2 or pattern2 in line):
            return True, line
    return False, None


def parse_msg_gps_fix(filename):
    '''Open and parse the msg file with the given name for a set of lines 
    starting with "# GPS fix obtained". 
    Extract lon/lat/time information and add it to the coords dictionary, which
    is returned.'''
    coords = dict()
    with open(filename) as fp:
        if scan_file_pattern(fp, 'GPS fix obtained')[0]:
            line = fp.readline() # this line has column headers only
            line = fp.readline().strip() # GPS fix data are 2 lines later
            regex = re.compile(r'^Fix:\s+([\d\.\-]+)\s+([\d\.\-]+)\s+(\d+/\d+/\d+\s+\d+)')
            match_obj = regex.search(line)
            if match_obj:
                coords['lon'] = match_obj.group(1)
                coords['lat'] = match_obj.group(2)
                coords['time'] = get_time_since_1970(match_obj.group(3))
                return coords
            
    # try an alternate pattern for time only
    with open(filename) as fp:
        success, line = scan_file_pattern(fp, 'Profile', 'terminated')
        if success:
            print(line)
            regex = re.compile(r'Profile.*terminated:\s*(.+)')
            match_obj = regex.search(line)
            if match_obj:
                coords['time'] = get_time_since_1970(match_obj.group(1))
                print('lon/lat not found')
            else:
                print(line)
                print('Error: time not found')
                pdb.set_trace()
        else:
            print('Error: GPS fix/Profile terminated with lon/lat/time not found')
            
    for crd in ['lon', 'lat', 'time']:
        if crd not in coords:
            coords[crd] = np.nan
    return coords


def parse_msg_header(fp, vars):
    '''Parse the header of the msg file with file pointer fp and fill
    dictionary vars with all successfully parsed variable names and their
    values and units.'''
    # the first line does not contain engineering variables
    fp.readline()

    regex1 = re.compile(r'\$\s+([\w]+)\((.+)\)\s+\[([\w]+)\]') # with units
    regex2 = re.compile(r'\$\s+([\w]+)\((.+)\)') # no units
    line = fp.readline().strip()
    # A line with only a '$' at its beginning marks the end of the header
    while line.strip() != '$':
        match_obj = regex1.search(line)
        if match_obj:
            vars[match_obj.group(1)] = (match_obj.group(2),match_obj.group(3))
        else:
            match_obj = regex2.search(line)
            if match_obj:
                vars[match_obj.group(1)] = (match_obj.group(2),'')
            else: #FIXME
                print(line)
                print('NO MATCH') # DEBUG
                pdb.set_trace()
        line = fp.readline().strip()


def parse_msg_footer(fp, vars):
    '''Parse the footer of the msg file with file pointer fp and fill
    dictionary vars with all successfully parsed variable names and their
    values and units. Note that units follow the values immediately.'''
    # skip a few more lines; more engineering variables are listed
    # after a line that starts with "Apf*FwRev" or "Npf*FwRev"
    line = fp.readline().strip()
    last_pos = fp.tell()
    while not (('Apf' in line or 'Npf' in line) and 'FwRev' in line):
        line = fp.readline().strip()
        if fp.tell() > last_pos:
            last_pos = fp.tell()
        else: # end of file was reached, no more engineering info
            return
    # the remaining lines contain more engineering information
    line = fp.readline().strip()
    last_pos = fp.tell()
    regex1 = re.compile(r'([\w]+)=([\d\.\-]+)(.*)')
    regex2 = re.compile(r'([\w]+)\[([\d]+)\]=(.+)') # array variables
    regex3 = re.compile(r'([\w]+)=(.+)') # anything else
    array_vars = dict() # for variables that occur in multiple lines
    # variable names are on the lhs, rhs will always be a tuple
    while line and '<EOT>' not in line:
        #pdb.set_trace()
        if line.startswith('# GPS fix obtained'):
            #DEBUG print('GPS fix found in footer, skipping the rest') # FIXME
            break
        
        match_obj = regex1.search(line)
        if match_obj:
            vars[match_obj.group(1)] = (match_obj.group(2), match_obj.group(3))
        else:
            # special treatment for several lines like this:
            # ParkDescentP[0]=6
            match_obj = regex2.search(line)
            if match_obj:
                if match_obj.group(1) not in array_vars:
                    array_vars[match_obj.group(1)] = (list(),)
                # assume that they are always listed in ascending order of index
                array_vars[match_obj.group(1)][0].append(match_obj.group(3))
            else:
                match_obj = regex3.search(line)
                if match_obj:
                    vars[match_obj.group(1)] = (match_obj.group(2),)
                else: #FIXME
                    print(line)
                    print('NO MATCH') # DEBUG
                    #pdb.set_trace()
        line = fp.readline().strip()
        if fp.tell() > last_pos:
            last_pos = fp.tell()
        else: # end of file was reached, no more engineering info
            break

    for key in array_vars.keys():
        vars[key] = array_vars[key]


def get_time_since_1970(datetime_string):
    '''Convert a string in the form "Thu Nov 12 00:08:02 2020" to a time value
    in seconds since Jan 1, 1970 midnight and return it as a float.
    Also try "Nov 12 2020 00:08:02" as an alternate format.'''
    try:
        # e.g., "Thu Nov 12 00:08:02 2020"
        date_time_obj = datetime.datetime.strptime(datetime_string.strip(),'%c')
    except ValueError:
        try:
            # e.g., "Nov 12 2020 00:08:02"
            date_time_obj = datetime.datetime.strptime(datetime_string.strip(),
                                                       '%b %d %Y %H:%M:%S')
        except ValueError:
            # e.g., "07/20/2021 104130"
            date_time_obj = datetime.datetime.strptime(datetime_string.strip(),
                                                       '%m/%d/%Y %H%M%S')

    date_time_1970 = datetime.datetime(1970, 1, 1, 0, 0, 0)
    dt = date_time_obj - date_time_1970
    return dt.total_seconds()


def create_nc_file(filename_out, filename_in, vars, program, verbose):
    '''Create netcdf file "filename_out" if it does not exist yet 
    and write simple numbers from the "vars" dictionary into it. 
    FIXME: Lists and strings are skipped for now.'''
    if os.path.exists(filename_out):
        if verbose:
            print('"{0:s}" exists already!'.format(filename_out))
        return
    ncfile = nc.Dataset(filename_out, 'w', format='NETCDF4')
    # global attributes
    ncfile.history = 'Created by parse_msg_file.py'
    ncfile.source = filename_in
    ncfile.date = date.today().strftime('%Y-%m-%d %H:%M:%S')
    # floatid is a single variable, it doesn't need a dimension
    floatid = get_floatid(filename_in)
    floatid_var = ncfile.createVariable('floatid', np.int32, ())
    # same for the program type
    prog_var = ncfile.createVariable('Program', str, ())
    # time dimension and (xyt) grid variables
    time_dim = ncfile.createDimension('time', None)
    time_var = ncfile.createVariable('time', np.float64, ('time',))
    time_var.units = 'seconds since 01-JAN-1970 00:00:00'
    time_var.time_origin = '01-JAN-1970 00:00:00'
    time_var.calendar = 'gregorian'
    time_var.FillValue_ = np.nan
    lon_var = ncfile.createVariable('longitude', np.float32, ('time',))
    lon_var.units = 'degrees_east'
    lon_var.FillValue_ = np.nan
    lat_var = ncfile.createVariable('latitude', np.float32, ('time',))
    lat_var.units = 'degrees_north'
    lat_var.FillValue_ = np.nan

    floatid_var[:] = floatid
    ncfile['Program'][0] = program

    # all time-dependent variables
    for var in vars:
        if 'DialCmd' in var:
            continue
        # a few values are lists instead of strings
        if isinstance(vars[var][0], str):
            nc_var = ncfile.createVariable(var, np.float32, ('time',))
            nc_var.FillValue_ = np.nan           
            if var == 'TimeStartTelemetry':
                # the values is in seconds, the second entry is a full
                # date, e.g. "Nov 22 2020 02:07:56"
                # ERDDAP would choke on that if used as a unit
                nc_var.units = 'seconds since 01-JAN-1970 00:00:00'
            else:
                try:
                    nc_var.units = vars[var][1]
                except: # ParkObs and SurfaceObs have a different format
                    pass # FIXME 
    ncfile.close()


def write_nc_file(filename_out, vars, coords, verbose):
    if not os.path.exists(filename_out):
        print('"{0:s} does not exist!'.format(filename_out))
        return
    ncfile = nc.Dataset(filename_out, 'a')
    times = ncfile['time'][:]
    # check if new time is not yet present in the file
    if len(times) == 0 or coords['time'] > times[-1]:
        nt = len(times)
        ncfile['time'][nt] = coords['time']
        ncfile['longitude'][nt] = coords['lon']
        ncfile['latitude'][nt] = coords['lat']
        # all time-dependent variables
        for var in vars:
            #DEBUG print(var)
            try:
                nc_var = ncfile[var]
            except IndexError:
                if 'DialCmd' in var:
                    continue
                #DEBUG print('{0:s} not yet present in file'.format(var))
                #pdb.set_trace()
                nc_var = ncfile.createVariable(var, np.float32, ('time',))
                nc_var.FillValue_ = np.nan            
                try:
                    nc_var.units = vars[var][1]
                except:
                    pass # FIXME
            # a few values are lists instead of strings
            if isinstance(vars[var][0], str):
                try:
                    this_num = float(vars[var][0])
                    ncfile[var][nt] = this_num
                except ValueError:
                    if verbose:
                        print('cannot convert: {0:s}'.format(vars[var][0]))
    ncfile.close()
    

def get_floatid(filename_in):
    '''Determine the floatid from the given filename and return it as an int.'''
    regex = re.compile(r'([\d]+)\.[\d]+\.msg')
    match_obj = regex.search(filename_in)
    if match_obj:
        return int(match_obj.group(1))
    else:
        raise ValueError('could not determine floatid of file' +
                         '{0:s}'.format(filename_in))


if __name__ == '__main__':
    args = parse_input_args()
    change_cwd(args.directory)
    for file in args.filename_in:
        filename_out = get_filename_out(file, args.output_directory)
        if args.verbose:
            print('Processing "{0:s}", writing to "{1:s}"'.format(file,
                                                                  filename_out))
        (vars, coords, program) = parse_msg_file(file)
        if vars:
            create_nc_file(filename_out, file, vars, program, args.verbose)
            #if not coords['time']:
                # this may happen for the first profile
                # the actual value of TimeStartTelemetry is in its units field
                #coords['time'] = get_time_since_1970(vars['TimeStartTelemetry'][1])
            write_nc_file(filename_out, vars, coords, args.verbose)
        elif args.verbose:
            print('{0:s} not found, skipping...'.format(file))

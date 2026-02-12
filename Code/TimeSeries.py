"""
Author:             Kyle Enns
Created:            2/1/2018
License:            Creative Commons Attribution 4.0 International (CC BY 4.0)
                    http://creativecommons.org/licenses/by/4.0/
Python version:     Tested on Python 3.6x (x64)


PURPOSE
------------------------------------------------------------------------------
The Magnetotelluric(MT) ascii to edi Python software was developed by the U.S.
Geological Survey Data at Risk project (DAR)for the purpose of converting Ascii
files into the EDI format industry standard for MT interpretation software.


U.S. GEOLOGICAL SURVEY DISCLAIMER
------------------------------------------------------------------------------
This software has been approved for release by the U.S. Geological Survey
(USGS). Although the software has been subjected to rigorous review, the
USGS reserves the right to update the software as needed pursuant to further
analysis and review. No warranty, expressed or implied, is made by the USGS
or the U.S. Government as to the functionality of the software and related
material nor shall the fact of release constitute any such warranty.
Furthermore, the software is released on condition that neither the USGS nor
the U.S. Government shall be held liable for any damages resulting from its
authorized or unauthorized use.

Any use of trade, product or firm names is for descriptive purposes only and
does not imply endorsement by the U.S. Geological Survey.
------------------------------------------------------------------------------


NOTES
------------------------------------------------------------------------------
"""

# Import standard python modules.
import os
import sys
import datetime
import math
import csv
from datetime import date, timedelta

# Import pyproj modules.
from pyproj import Proj

# -----------------------Beginning of File Locations---------------------------

# INPUT LOCATIONS:
# Directory for Time series site directories.
#   ex: r'...Validation_Resources\Magnetotelluric_Data_Files\Magnetotelluric_Time-Series_ASCII'
ts_ = None

# Directory for Cross power site directories.
#   ex: r'...Validation_Resources\Magnetotelluric_Data_Files\Magnetotelluric_Cross-Power_ASCII'
cp_ = None

# Directory for all QC CSV files.
#   ex: r'...Validation_Resources\QC_Site-Sheets'
qcDir_ = None

# Geospatial CSV File.
#   ex: r'...Validation_Resources\Magnetotelluric_Site-Location-Information.csv'
locations_ = None

# OUTPUTS LOCATION:
# Points at an output folder for the Edi files.
edi_ = None


# ----------------------------Beginning of Functions---------------------------

def create_qc_dic(qc_dir):
    """
    Converts a csv into a dictionary that defines the quality control
      for each site.
    :param qcDir:
    :return:
    Dictionary with quality control information for each site
    """
    qcFiles = os.listdir(qc_dir)
    qcDict = {}
    for qcf in qcFiles:

        # qcDict['fileName'][0] = Rating, [1] = GoodFrom, [2] = GoodTo,
        #   [3] = Comments
        # Ex: qcDict['PARD02_2'] = ('4', '0.2924', '2.5126',
        #   'Stronger signal along y-azimuth.')
        with open(qc_dir + '/' + qcf, 'r') as csvfile:
            reader = csv.reader(csvfile)
            for row in reader:
                if row[0].startswith('PAR') is True:
                    qcDict[row[1].strip()] = row[2], row[3], row[4], row[5]
        csvfile.close()
    return qcDict


def create_coord_dic(locations):
    """
    This is a dictionary that contains the coordinates for the different
      sites where the USGS MT truck was set up to
      measure the Electric and Magnetic fields.
    :param locations:
    :return:
    """
    lat, long, elev = {'Latitude': {}}, {'Longitude': {}}, {'Elevation': {}}
    with open(locations, 'r') as locFile:
        reader = csv.reader(locFile)
        count = 0
        for row in reader:
            count += 1
            if count > 1:
                if len(row[1]) < 2:
                    site = 'site0' + str(row[1])
                elif len(row[1]) == 2:
                    site = 'site' + str(row[1])
                long['Longitude'].update({site: float(row[2])})
                lat['Latitude'].update({site: float(row[3])})
                elev['Elevation'].update({site: float(row[5])})
    locFile.close()
    lat_long = [long, lat, elev]
    return lat_long


def info1(site_, types_, band_, run_, science_base, file_name):
    """
    EDI file Description.

    :param site_: site of time series file
    :param types_: time series or Cross power
    :param band_: high, medium, low, or lower
    :param run_: run number
    :param science_base: sciencebase web page url
    :param file_name: the name of the input file
    :return:
    Block string
    """

    info = ''' Description=Magnetotelluric (MT) Time Series (TS) Metadata

 ProductId=USA-California-San_Andreas_Fault-Parkfield-1990-MT''' \
           + site_ + types_ + '''_''' + band_ + '''_''' + run_ + '''.edi

 Tags=timeseries

 ExternalUrl Description=U.S. Geological Survey (USGS) MetaData
 ExternalUrl Url=''' + science_base + '''

 PrimaryData GroupKey=0
 PrimaryData OrderKey=0

 Attachment Filename=''' + file_name + '''
 Attachment Description=The original USGS MT Truck System time series file

 =====================================================================\n\n'''
    return info


def info2():
    """
    Survey purpose.

    :return:
    Block string
    """

    info = ''' Survey Purpose Description:
 This dataset includes the magnetotelluric (MT) sounding data collected in 1990
 across the San Andreas fault, Parkfield, California. The San Andreas fault at
 Parkfield, California has been and still is a region of global importance for
 earthquake prediction studies, since it provides the shortest cyclical
 earthquake occurrences in the world. The U.S. Geological Survey conducted the
 MT survey to improve understanding of the San Andreas fault near Parkfield,
 California.

 Data Description:
 The purpose of this data release is to provide MT station locations and data
 to the public. The primary goal of the MT survey is to map changes in
 electrical resistivity with depth that are related to differences in rock
 type. The structural geologic setting across the fault can be characterized
 by these variations in rock types.

 =====================================================================\n\n'''
    return info


def info3(todays_date, conversion_file):
    """
    Time this script runs and creates the Time series EDI files.

    :param todays_date: YYYY-MM-DD HH:mm:ss UTC-0700
    :param conversion_file: Python File Conversion Utility x.x
    :return:
    Block string
    """
    info = ''' CreateTime=''' + todays_date + '''

 CreatingApplication=''' + conversion_file + '''
 Source Code Author=Kyle Enns, USGS, Fort Collins, Colorado, USA

 =====================================================================\n\n'''
    return info


def info4():
    """
    File Creator for Primary(timeseries.ts) and secondary(edi) data.

    :return:
    Block string
    """

    info = ''' FILE CREATOR:
 Creator Name=Kyle Enns
 Creator Email=kenns@usgs.gov
 Creator Org=U.S. Geological Survey
 Creator OrgUrl=http://www.fort.usgs.gov/

 FILE SUBMITTER:
 Submitter Name=Brian D. Rodriguez
 Submitter Email=brod@usgs.gov
 Submitter Org=U.S. Geological Survey
 Submitter OrgUrl=http://crustal.usgs.gov/

 ReleaseStatus=Unrestricted Release

 =====================================================================\n\n'''
    return info


def info5():
    """
    USGS Conditions of Use.

    :return:
    Block string
    """
    info = ''' ConditionsOfUse=All data and metadata for this survey are
 available free of charge and may be copied freely, duplicated and further
 distributed provided this data set is cited as the reference. While the
 author(s) strive to provide data and metadata of best possible quality,
 neither the author(s) of this data set, nor U.S. Geological Survey make any
 claims, promises, or guarantees about the accuracy, completeness, or adequacy
 of this information, and expressly disclaim liability for errors and omissions
 in the contents of this file. Guidelines about the quality or limitations of
 the data and metadata, as obtained from the author(s), are included for
 informational purposes only.

 Users of this dataset must agree to the following two conditions:

 (1) he/she is required to express acknowledgements to the authors and the U.S.
 Geological Survey in his/her publication (including proceedings) and

 (2) he/she is required to not give these data to anyone. All users should get
 the data from the U.S. Geological Survey, so that everyone has the same data -
 it is easy to "edit" data and forget that you have done so, and then pass it
 on to someone else as the original dataset.

 =====================================================================\n\n'''
    return info


def info6(site_id, site_lat, site_long, elevation_, site_x_north, site_y_east,
          start_, end_, rating_, good_from, good_to, comments_):
    """
    Project and site information used for the Spud format.

    :param site_id: site location number
    :param site_lat: site latitude
    :param site_long: site longitude
    :param elevation_: site elevation
    :param site_x_north: X north
    :param site_y_east: Y east
    :param start_: when collection of data for the file began
    :param end_: when collection of data for the file ended
    :param rating_: QC for time series data collection
    :param good_from: QC for time series data collection
    :param good_to: QC for time series data collection
    :param comments_: comments
    :return:
    Block string
    """

    info = ''' Project=USGS California Coast Range Geoelectrical Studies
 Survey=San Andreas fault, CA
 YearCollected=1990
 Country=USA
 Site Id=''' + site_id + '''
 Area Name=Parkfield
 Ellipsoid=Clarke 1866
 Location datum=NAD27 CONUS
 SITE LATITUDE=''' + str(site_lat) + '''
 SITE LONGITUDE=''' + str(site_long) + '''
 Elevation units="meters"=''' + str(elevation_) + '''
 SITE XNORTH=''' + str(site_x_north) + '''
 SITE YEAST=''' + str(site_y_east) + '''
 UTM zone=Zone 10 North
 Declination epoch="1990.0"=15.0
 AcquiredBy=USGS
 Start=''' + str(start_) + '''
 End=''' + str(end_) + '''
 RunList=''' + str(site_id) + '''
 DataQualityNotes Rating=''' + str(rating_) + '''
 DataQualityNotes GoodFromPeriod=''' + str(good_from) + '''
 DataQualityNotes GoodToPeriod=''' + str(good_to) + '''
 DataQualityNotes Comments author=Brian D. Rodriguez, ''' + str(comments_) + '''
 DataQualityWarnings Flag=0

 =====================================================================\n\n'''
    return info


def info7(magnetometer_, azm_1, azm_2, azm_5, azm_6, manufacturer_1,
          manufacturer_2, dipole_length_ex, dipole_length_ey, hx_1, hx_2,
          hy_1, hy_2):
    """
    Instrument Manufacture information and parameters used for the
    spud format.

    :param magnetometer_: name of the magnetometer
    :param azm_1: sensor geographic azimuth units for Hx1 magnetometer
    :param azm_2: sensor geographic azimuth units for Hx2 magnetometer
    :param azm_5: sensor geographic azimuth units for dipole Ex
    :param azm_6: sensor geographic azimuth units for dipole Ey
    :param manufacturer_1: pre-amp manufacturer id for Ex-1
    :param manufacturer_2: pre-amp manufacturer id for Ey-1
    :param dipole_length_ex: dipole length in meters for Ex
    :param dipole_length_ey: dipole length in meters for Ey
    :param hx_1: magnetometer id for Hx-1
    :param hx_2: magnetometer id for Hx-2
    :param hy_1: magnetometer id for Hy-1
    :param hy_2: magnetometer id for Hy-2
    :return:
    Block string
    """

    info = ''' Instrument Manufacturer=USGS, Denver, CO, USA
 Instrument Name=INFOTEK A/D
 Instrument Id=1985

 Magnetometer type="induction"
 Magnetometer Manufacturer=USGS, Denver, CO, USA
 Magnetometer Name=''' + str(magnetometer_) + '''
 Magnetometer Hx-1  Id=''' + str(hx_1) + '''
 Magnetometer Hx-1  Sensor Geographic Azimuth units="degrees"=''' + str(azm_1) + '''
 Magnetometer Hx-2  Id=''' + str(hx_2) + '''
 Magnetometer Hx-2  Sensor Geographic Azimuth units="degrees"=''' + str(azm_1) + '''
 Magnetometer Hy-1  Id=''' + str(hy_1) + '''
 Magnetometer Hy-1  Sensor Geographic Azimuth units="degrees"=''' + str(azm_2) + '''
 Magnetometer Hy-2  Id=''' + str(hy_2) + '''
 Magnetometer Hy-2  Sensor Geographic Azimuth units="degrees"=''' + str(azm_2) + '''

 Dipole name="Ex" type="wire"
 Dipole Manufacturer=Tinker & Rasor, San Gabriel, CA, USA
 Dipole Length units="meters"= ''' + str(dipole_length_ex) + '''
 Dipole Sensor Geographic Azimuth units="degrees"=''' + str(azm_5) + '''
 Electrode location="N" number="unknown"=CuCu-SO4 Fat Boy Porous Pot
 Electrode location="S" number="not used"=not applicable
 Ex-1  Pre-amp Manufacturer Id=''' + str(manufacturer_1) + '''

 Dipole name="Ey" type="wire"
 Dipole Manufacturer=Tinker & Rasor, San Gabriel, CA, USA
 Dipole Length units="meters"= ''' + str(dipole_length_ey) + '''
 Dipole Sensor Geographic Azimuth units="degrees"=''' + str(azm_6) + '''
 Electrode location="E" number="unknown"=CuCu-SO4 Fat Boy Porous Pot
 Electrode location="W" number="not used"=not applicable
 Ey-1  Pre-amp Manufacturer Id=''' + str(manufacturer_2) + '''

 =====================================================================\n\n'''
    return info


def info8():
    """
    Authors comments and/or errors

    :return:
    Block string
    """

    info = ''' Comments author=none.

 Errors=none.

 =====================================================================\n\n'''
    return info


def info9(freq_, start_, end_, lp_band, hp_band, ts_length, nch_,
          gain_, sensitivity_):
    """
    Information about the individual ascii TS file (Number of channels,
      gain, and sensitivities for those channels) as well as acquisition
      software information.

    :param freq_: sampling rate in hz
    :param start_: when collection of data for the file began
    :param end_: when collection of data for the file ended
    :param lp_band: LP filter
    :param hp_band: HP filter
    :param ts_length: FFT length
    :param nch_: number of channels
    :param gain_: the gain values for each channel
    :param sensitivity_: the sensitivity values for each channel
    :return:
    Block string
    """

    info = ''' SamplingRate units="Hz":''' + str(freq_) + '''
 Start=''' + str(start_) + '''
 End=''' + str(end_) + '''

 Number of Time Series per Block= 1
 Number of Blocks= 1
 LP Filter= ''' + str(lp_band) + '''
 HP Filter= ''' + str(hp_band) + '''
 FFT Length= ''' + str(ts_length) + '''

 Number of Channels= ''' + str(nch_) + '''
 Channel 1=Hx-1
 Channel 2=Hx-2
 Channel 3=Ey-1
 Channel 4=Hy-1
 Channel 5=Hy-2
 Channel 6=Ex-1

''' + gain_ + '''

''' + sensitivity_ + '''

 AcquisitionSoftware Name=MTR02B.HPB
 AcquisitionSoftware LastMod=1986-06-06
 AcquisitionSoftware Author=Vic F. Labson, USGS, Denver, CO, USA

 =====================================================================\n\n'''
    return info


# Acquisition hardware and software
def info10():
    """
    USGS MT Truck System Time Series Format Description

    :return:
    Block string
    """

    info = ''' USGS MT Truck System Time Series Format Description:

 The time series file begins with the time series file name in HP BASIC format,
 followed by header information on some of the data acquisition settings
 (center frequency, number of channels used, number of time series data
 points sampled, gains for 10 channels, sensitivities for 10 channels,
 and ends with the time series data values. The corresponding cross-power
 file includes the balance of the data acquisition settings (calibration
 file names for sensor corrections and sensor direction) needed to convert
 the time series data into the frequency domain.

 The time series data is listed as one line per channel for each time series
 data point sampled. In other words,

 time series data value 1 = channel 1 of sample 1
 time series data value 2 = channel 2 of sample 1
 time series data value 3 = channel 3 of sample 1
                        .
                        .
                        .
 time series data value Nch = channel Nch of sample 1
 time series data value Nch+1 = channel 1 of sample 2
 time series data value Nch+2 = channel 2 of sample 2
 time series data value Nch+3 = channel 3 of sample 2
                        .
                        .
                        .
 time series data value Nch x N = channel Nch of sample N

 (where Nch is the total number of channels and N is the total number
 of samples)

 The total number of time series data records is the total number of
 channels recorded multiplied by the total number of time series data
 points sampled.

 The USGS MT time series data is written out as voltage (in units of
 millivolts or mV) with no corrections applied for sensors, dipole
 lengths, direction, or gain. The following steps are required to convert
 the voltage data into MT units - (nT for magnetic field  or mV/km for
 electric field).

 This voltage value (mV) corresponds to the voltage on each channel
 sample/hold amplifier. To convert to sensor units, each channel must
 first be corrected for the channel gain (G) applied by the system
 acquisition board (1, 2, 5, 10, 20, 50, 100, 200, 500 or 1000).

 The input voltage (mVinp) is thus:

    mVinp= mV/G

 In addition, the units of the electric field are normally presented in mV/km.
 The voltages must be divided by the dipole length (L) in meters. The
 sensitivity value of each channel stored in the header is the product of
 the dipole length in units of megameters (Mm). The electric field (E) in
 mV/km is thus:

    E = mVinp/(1000 x sensitivity)

 =====================================================================\n\n'''
    return info


def dd_to_dms(dd):
    """
    Changes geospatial values from decimal degrees to degree minutes seconds.

    :param dd: input of either latitude or longitude
    :return:
    Degree minutes seconds.
    """

    is_positive = dd >= 0
    dd = abs(dd)
    minutes, seconds = divmod(dd * 3600, 60)
    degrees, minutes = divmod(minutes, 60)
    degrees = degrees if is_positive else -degrees
    return degrees, minutes, seconds


# ---------------------------End of Functions---------------------------------
# ------------------Beginning of variables and workflow-----------------------

# Sets variables that don't change based on the different TS files inserted.
def main(ts, cp, edi, qc_dir, locations):
    """
    Uses the 5 dir/files to run all the appropriate functions to output
      all the necessary information for the Time series EDI files.

    :return:
    1 Time series EDI file per input ASCII file
    """

    indent = '  '
    nl = '\n'
    maxLines = '1000'
    equal = ('================================================================='
             + '========')
    conversionFile = 'Python File Conversion Utility 1.0'
    locName = '\\USA-California-San_Andreas_Fault-Parkfield-1990-MT'
    sbUrl = 'http://doi.org/10.5066/F71Z43P9'
    stdVers = '1.0'
    progVers = '1.0'

    qcDict = create_qc_dic(qc_dir)
    LatLong = create_coord_dic(locations)

    if 'TS_edi' not in os.listdir(edi):
        os.mkdir(edi + '\\TS_edi')

    sites = os.listdir(ts)

    # Loops through all of the site folders
    for s in sites:

        if s not in os.listdir(edi + '\\TS_edi'):
            os.mkdir(edi + '\\TS_edi\\' + s)

        # Define variables for time series(TS) and cross power(CP) sites.
        fileNames = []                               # List of all file names (TS).
        location = 0                                 # Keeps track of file count.
        run = 0                                      # The file's run number.
        timeSeries = ts + '\\' + s                   # TS site folder.
        crossPower = cp + '\\' + s                   # CP site folder.
        tsFiles = os.listdir(timeSeries)             # List of all the ts files.
        cpFiles = os.listdir(crossPower)             # List of all the cp files.
        countL, countD, countX, countK = 0, 0, 0, 0  # Band counts.
        addTimeD = 0                                 # Tracks count of ts files.
        # Prints the site where this TS file is located.
        print(s)

        # Below re-orders the files to ascending, which allows the script to add
        #   appropriate start and stop times for the their edi counterpart file.
        #   .TS and T.S were supposed to be collected as .TS but the single digit
        #   Medium band (D(which stands for Dead band)) had the 1-9 runs.  Truck
        #   Software/Hardware only captured one digit which lead to the files
        #   ending in '.TS' and 'T.S'.

        # Creating a blank slate for TS files.
        for d in range(0, len(tsFiles)):
            fileNames.append(' ')

        for f in tsFiles:
            if f[-3:] == '.TS':
                run = f[6:8]
            elif f[-3:] == 'T.S':
                run = f[6]

            if f[3] == 'D':
                fileNames[int(run) - 1] = f
                location += 1
            elif f[3] == 'K':
                fileNames[location] = f
                location += 1
            elif f[3] == 'L':
                fileNames[location] = f
                location += 1
            elif f[3] == 'X':
                fileNames[location] = f
                location += 1

        if run == 0:
            print('This doesn\'nt end in .TS or T.S')
        else:
            print(fileNames)
        # Now the script runs through each .TS file and creates a list within
        #   python that will store the file line by line into a list we can
        #   pull values from.
        for n in fileNames:
            # Define variables used to build TS specific output edi file name.
            types = '-TS'                                # TS name.
            tsList = []                                  # TS data list.
            cpList = []                                  # CP data list.
            count = 0                                    # Rows in TS ascii file.
            tsFile = open((timeSeries + '\\' + n), 'r')  # TS ascii file.
            tsRunWrite = 0                               # TS run number.

            for row in tsFile:
                tsList.append(row.strip())
                count += 1
            tsFile.close()

            # Checks to see which Time Series file it is (ending with .TS or T.S)
            #   to store values from the name (ex: PARD0121.TS) that tell us
            #   the Run #, Band, and Site information for the TS file.

            if n[-3:] == '.TS':
                # Extracting TS file name into name, band, run, and site number
                #   from the time series file name.
                name, tsBand, tsRun, tsSite = \
                      n.split('.')[0], n[3], n[6:8], n[4:6]
                tsRunWrite = str(tsRun)

                # Labeling tsBand as Mid,Low,or High based on TS file name.
                if tsBand == 'D':
                    band = 'Mid'
                elif tsBand == 'L' or 'X':
                    band = 'Low'
                elif tsBand == 'K':
                    band = 'High'

                # Creating output TS edi file name.
                ediFile = open(edi + '\\TS_edi' + '\\' + s
                               + locName + str(tsSite) + types + '_'
                               + band + '_' + tsRun + '.edi', 'w')

            elif n[-3:] == 'T.S':
                name, tsBand, tsRun, tsSite = \
                      n.split('.')[0][0:-1], n[3], n[6], n[4:6]
                tsRunWrite = '0' + str(tsRun)

                if tsBand == 'D':
                    band = 'Mid'
                elif tsBand == 'L' or 'X':
                    band = 'Low'
                elif tsBand == 'K':
                    band = 'High'

                ediFile = open(edi + '\\TS_edi' + '\\' + s
                               + locName + str(tsSite) + types + '_'
                               + band + '_' + tsRunWrite + '.edi', 'w')

            # Adds the Quality Control values associated with the file
            #   band (from the qcDict created from the Qc files) and
            #   run number.
            if tsBand == 'D':
                countD += 1
                rating, goodFrom, goodTo, comments = \
                    'unknown', 'n/a', 'n/a', 'n/a'

            elif tsBand == 'K':
                countK += 1
                rating, goodFrom, goodTo, comments = \
                    'unknown', 'n/a', 'n/a', 'n/a'

            elif tsBand == 'L':
                countL += 1
                rating, goodFrom, goodTo, comments = \
                    qcDict[name + 'ts'][0], qcDict[name + 'ts'][1], \
                    qcDict[name + 'ts'][2], qcDict[name + 'ts'][3]

            elif tsBand == 'X':
                countX += 1
                rating, goodFrom, goodTo, comments = \
                    qcDict[name + 'ts'][0], qcDict[name + 'ts'][1], \
                    qcDict[name + 'ts'][2], qcDict[name + 'ts'][3]
            else:
                rating, goodFrom, goodTo, comments = \
                    'unknown', 'n/a', 'n/a', 'n/a'

            # Loops through the associated CP files for the Site and matches the
            #   TS files to their corresponding CP file.  Now we can grab essential
            #   information from the CP files that don't live in the TS file.
            medBand = False
            for a in cpFiles:
                cpBand, cpRun, cpSite = a[3], int(a[7]), a[4:6]
                if cpBand == 'D' and tsBand == 'D':
                    if ((cpRun - 1) * 8) + 1 <= int(tsRun) <= (int(cpRun) * 8):
                        asciiFileCP = open(crossPower + '\\' + a, 'r')
                        for row in asciiFileCP:
                            cpList.append(row.strip())
                        lpFilter, hpFilter, medBand = '2.0', '0.1', True
                        asciiFileCP.close()

                elif cpBand == 'K' and tsBand == 'K':
                    if int(tsRun) - countD == cpRun:
                        asciiFileCP = open(crossPower + '\\' + a, 'r')
                        for row in asciiFileCP:
                            cpList.append(row.strip())
                        lpFilter, hpFilter = '202.0', '1.0'
                        asciiFileCP.close()

                elif cpBand == 'L' and tsBand == 'L':
                    if int(tsRun) - (countD + countK) == cpRun:
                        asciiFileCP = open(crossPower + '\\' + a, 'r')
                        for row in asciiFileCP:
                            cpList.append(row.strip())
                        lpFilter, hpFilter = '0.202', '0.003'
                        asciiFileCP.close()

                elif cpBand == 'X' and tsBand == 'X':
                    if int(tsRun) - (countD + countK + countL) == cpRun:
                        asciiFileCP = open(crossPower + '\\' + a, 'r')
                        for row in asciiFileCP:
                            cpList.append(row.strip())
                        lpFilter, hpFilter = '0.202', '0.003'
                        asciiFileCP.close()

            # Setting most of the variables for the Head/Info Sections that are
            #   found in the TS/CP files. These variables come from either tsList
            #   (which is the time series file) or cpList (which is the Cross
            #   Power file associated with that Time Series File).

            # Labelling magnetometer id's for each site.
            if str(tsSite) in ['01', '02', '03', '04', '05', '11']:
                hx1, hx2, hy1, hy2 = 1, 5, 7, 9
            elif str(tsSite) in ['06']:
                hx1, hx2, hy1, hy2 = 1, 9, 5, 7
            elif str(tsSite) in ['07']:
                hx1, hx2, hy1, hy2 = 1, 5, 7, 4
            elif str(tsSite) in ['08', '09']:
                hx1, hx2, hy1, hy2 = 1, 4, 5, 7
            else:
                hx1, hx2, hy1, hy2 = 1, 5, 7, 9

            # Sets todays date and start and end times for the
            #    TS files (Time variables).
            todaysDate = datetime.datetime.now()  # Set correct format.
            d = datetime.datetime.strptime(cpList[0], '%H:%M:%S %d %b %Y')
            if medBand is True:
                startTime = d + timedelta(seconds=int((int(tsList[3])
                                          / (4 * float(tsList[1])) * addTimeD)))
                endTime = startTime + timedelta(seconds=int((int(tsList[3]))
                                                / (4 * (float(tsList[1])))))
                addTimeD += 1
            else:
                startTime = d
                endTime = startTime + timedelta(seconds=int((int(tsList[3]))
                                                / (4 * (float(tsList[1])))))

            # Sets data ID, contractor, who created the CP file, and
            #   the date when this file is created.
            dataID = tsList[0]
            contractor, fileby = 'USGS', '"Jay Sampson"'
            acquireDate = d.strftime('%m/%d/%Y')

            # Reformat the start and end date.
            start = str(startTime.strftime('%Y-%m-%d')) + 'T' \
                + str(startTime.strftime('%H:%M:%S')) + ' UTC/GMT'
            end = str(endTime.strftime('%Y-%m-%d')) + 'T' \
                + str(endTime.strftime('%H:%M:%S')) + ' UTC/GMT'

            # Location variables
            country = 'USA'
            state = 'California'
            if tsSite == '09':
                county = 'Fresno'
            else:
                county = 'Monterey'

            # Conversion of Spatial Data variables (lat long) from
            #   decimal degrees(DD) to degrees minutes seconds (DMS)qcDict
            units = 'M'
            elevation = LatLong[2]['Elevation']['site' + str(tsSite)]
            latDD = LatLong[1]['Latitude']['site' + str(tsSite)]
            latDMS = dd_to_dms(latDD)
            latDMSf = str(int(latDMS[0])) + ':' + str(int(latDMS[1])) \
                + ':' + str(int(latDMS[2]))
            longDD = LatLong[0]['Longitude']['site' + str(tsSite)]
            longDMS = dd_to_dms(longDD)
            longDMSf = str(int(longDMS[0])) + ':' + str(int(longDMS[1])) \
                + ':' + str(int(longDMS[2]))

            # Number of channels, number of stacked TS files, and SR (a TS
            #   value used in TSERIESSECT section).
            nch = int(tsList[2])
            tslength = tsList[3]
            SR = float(tsList[1])

            # Using pyproj module project lat/long values into UTM zone 10s.
            myProj = Proj("+proj=utm +zone=10S, +north +ellps=clrk66 "
                          + "+datum=NAD27 +units=m +no_defs")
            lon, lat = myProj(longDD, latDD)

            # Calculates length of dipole along EX and EY.
            dipoleLengthEX = int(1E6 * float(cpList[5]))
            dipoleLengthEY = int(1E6 * float(cpList[8]))

            # Sets magnetometer, manufacturers, and azm cardinal directions
            #   for Ex,Ey and Rx, Ry.
            magnetometer = cpList[150]
            manufacturer1, manufacturer2 = cpList[155], cpList[152]
            azm1, azm2, azmEx, azmEy = \
                cpList[163], cpList[164], int(cpList[161]), int(cpList[162])

            # Using math module to calculate Ex and Ey distances based on
            #   length of dipole and cardinal direction (0-360 degrees).
            exX2 = int(math.sin(math.radians(azmEx)) * dipoleLengthEX)
            exY2 = int(math.cos(math.radians(azmEx)) * dipoleLengthEX)
            eyX2 = int(math.sin(math.radians(azmEy)) * dipoleLengthEY)
            eyY2 = int(math.cos(math.radians(azmEy)) * dipoleLengthEY)

            # Creates a gains and sensitivity string for part of the INFO section.
            #   num1/2 are to count up and keep track of channels.
            num1, num2, gain, sens, gainList = 0, 0, '', '', []
            for row in tsList[4:(4 + nch)]:
                num2 += 1
                gain += (' Channel ' + str(num2) + ' GAIN=' + str(row) + nl)
                gainList.append(row)
            for row in tsList[14:(14 + nch)]:
                num1 += 1
                sens += (' Channel ' + str(num1) + ' Sensitivity=' + str(row) + nl)

            # Header for the file.
            #   >HEAD EDI Section.
            ediFile.write('>HEAD\n\n')
            ediFile.write(indent + 'DATAID=' + dataID + nl)
            ediFile.write(indent + 'ACQBY=' + contractor + nl)
            ediFile.write(indent + 'FILEBY=' + fileby + nl)
            ediFile.write(indent + 'ACQDATE=' + acquireDate[0:6]
                          + acquireDate[-2:] + nl)
            ediFile.write(indent + 'FILEDATE=' + acquireDate[0:6]
                          + acquireDate[-2:] + nl)
            ediFile.write(indent + 'ENDDATE=' + acquireDate[0:6]
                          + acquireDate[-2:] + nl)
            ediFile.write(indent + 'COUNTRY=' + country + nl)
            ediFile.write(indent + 'STATE=' + state + nl)
            ediFile.write(indent + 'COUNTY=' + county + nl)
            ediFile.write(indent + 'LAT=' + latDMSf + nl)
            ediFile.write(indent + 'LONG=' + longDMSf + nl)
            ediFile.write(indent + 'UNITS=' + units + nl)
            ediFile.write(indent + 'STDVERS=' + stdVers + nl)
            ediFile.write(indent + 'PROGVERS=' + progVers + nl)
            ediFile.write(indent + 'PROGDATE='
                          + str(todaysDate.strftime('%m/%d/%Y')) + nl)
            ediFile.write(nl)

            # Information about the file.
            #   >INFO EDI Section.
            ediFile.write('>INFO   MAXLINES=' + maxLines + nl + nl)

            # Writes to the output edi file all of the info sections from the
            #   info functions above.
            ediFile.write(info1(tsSite, types, band, tsRunWrite, sbUrl, n))
            ediFile.write(info2() + info3(todaysDate.strftime('%Y/%m/%d %H:%M:%S ')
                                          + 'UTC-0700', conversionFile))
            ediFile.write(info4() + info5())
            ediFile.write(info6(name, latDD, longDD, elevation, round(lat, 2),
                                round(lon, 2), start, end, rating, goodFrom,
                                goodTo, comments))
            ediFile.write(info7(magnetometer, azm1, azm2, azmEx, azmEy,
                                manufacturer1, manufacturer2, dipoleLengthEX,
                                dipoleLengthEY, hx1, hx2, hy1, hy2))
            ediFile.write(info8())
            ediFile.write(info9(tsList[1], start, end, lpFilter, hpFilter,
                                tsList[3], nch, gain, sens))
            ediFile.write(info10())

            # Defining measurements.
            #   >=DEFINEMEAS EDI Section.
            ediFile.write('>=DEFINEMEAS' + nl + nl)
            ediFile.write(indent + 'MAXCHAN=7' + nl)
            ediFile.write(indent + 'MAXRUN=999' + nl)
            ediFile.write(indent + 'MAXMEAS=99999' + nl)
            ediFile.write(indent + 'UNITS=M' + nl)
            ediFile.write(indent + 'REFLAT=' + latDMSf + nl)
            ediFile.write(indent + 'REFLONG=' + longDMSf + nl + nl)

            if len(tsSite) == 2 and tsSite[0] == '0':
                tsSite = tsSite[1]

            # Writes out the >HMEAS lines for each channel to output edi file.
            ediFile.write('>HMEAS ID=  ' + str(tsSite) + str(1) + '.'
                          + '001' + ' ' + 'CHTYPE=HX  ' + 'X=  0 Y=  0 AZM=  '
                          + azm1 + ' ' + 'GAIN=' + gainList[0] + nl)
            ediFile.write('>HMEAS ID=  ' + str(tsSite) + str(2) + '.'
                          + '001' + ' ' + 'CHTYPE=HX  ' + 'X=  0 Y=  0 AZM=  '
                          + azm1 + ' ' + 'GAIN=' + gainList[1] + nl)
            ediFile.write('>EMEAS ID=  ' + str(tsSite) + str(3) + '.'
                          + '001' + ' ' + 'CHTYPE=EY  ' + 'X=  0 Y=  0 X2=  '
                          + str(eyX2) + ' Y2= ' + str(eyY2) + ' GAIN='
                          + gainList[2] + nl)
            ediFile.write('>HMEAS ID=  ' + str(tsSite) + str(4) + '.'
                          + '001' + ' ' + 'CHTYPE=HY  ' + 'X=  0 Y=  0 AZM=  '
                          + azm2 + ' ' + 'GAIN=' + gainList[3] + nl)
            ediFile.write('>HMEAS ID=  ' + str(tsSite) + str(5) + '.' + '001'
                          + ' ' + 'CHTYPE=HY  ' + 'X=  0 Y=  0 AZM=  ' + azm2
                          + ' ' + 'GAIN=' + gainList[4] + nl)
            ediFile.write('>EMEAS ID=  ' + str(tsSite) + str(6) + '.' + '001'
                          + ' ' + 'CHTYPE=EX  ' + 'X=  0 Y=  0 X2=  ' + str(exX2)
                          + ' Y2= ' + str(exY2) + ' GAIN=' + gainList[5] + nl + nl)

            # Creating data values in arrays.
            #   >=TSERIESSECT EDI Section.
            ediFile.write('>=TSERIESSECT' + nl)
            ediFile.write(indent + 'SECTID=' + '-' + str(tsSite) + nl)
            ediFile.write(indent + 'NCHAN=' + str(nch) + nl)
            ediFile.write(indent + 'MAXBLKS=' + str(100) + nl)
            ediFile.write('  //' + str(nch) + nl)
            ediFile.write('  ')

            # Builds the arrays and their parameters to meet EDI standard.
            for i in range(0, int(nch)):
                ediFile.write(str(tsSite) + str(i + 1) + '.' + '001' + '  ')
            ediFile.write(nl + nl)

            ediFile.write('>TSERIES SECTID= -' + str(tsSite) + ' NPTS= '
                          + str(tslength) + ' ' + 'SR= ' + str(SR) + ' '
                          + 'MPX= TIME' + ' //'
                          + str(int(tslength) * int(nch)) + nl)

            add, end = 1, 0
            for x in tsList[24:]:
                if type(x) is str:
                    end += 1
                    if end == 1:
                        ediFile.write('      ')
                        col = 9 - len(x)
                        ediFile.write(x + (' ' * col))
                    elif end < 6:
                        col = 9 - len(x)
                        ediFile.write(x + (' ' * col))
                    elif end == 6:
                        col = 9 - len(x)
                        ediFile.write(x + (' ' * col) + nl)
                        end = 0
                        add += 1

            ediFile.write(nl + nl)
            ediFile.write('>END')

            # print(rating, goodFrom, goodTo, comments)
            # print('NewList ', cpList[1])
            # print(startTime, nl, endTime)
            # print(tsList[0])
            # print('Lines in TS file: ', count, nl, equal)

        ediFile.close()

if __name__ == '__main__':
    for item in [ts_, cp_, edi_, qcDir_, locations_]:
        if item is None:
            print ('Error: One of the inputs is None.')
            break
    main(ts_, cp_, edi_, qcDir_, locations_)

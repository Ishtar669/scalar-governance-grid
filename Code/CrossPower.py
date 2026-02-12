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
from decimal import Decimal
from datetime import date, timedelta

# Import numpy module.
import numpy as np

# Import pyproj module.
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

    info = ''' Description=Magnetotelluric (MT) Transfer Functions (TF)

 ProductId=USA-California-San_Andreas_Fault-Parkfield-1990-MT''' \
           + site_ + types_ + '''_''' + band_ + '''_''' + run_ + '''.edi

 SubType=MT_TF

 Tags=impedance,tipper,apparent_resistivity,impedance_phase,impedance_strike,
 impedance_skew,coherency,signal,noise

 ExternalUrl Description=U.S. Geological Survey (USGS) MetaData
 ExternalUrl Url=''' + science_base + '''

 PrimaryData Filename=''' + file_name + '''.jpg

 PrimaryData GroupKey=0
 PrimaryData OrderKey=0

 Attachment Filename=''' + file_name + '''
 Attachment Description=The original USGS MT Truck system cross-power file

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
    :param rating_: QC for cross power data collection
    :param good_from: QC for cross power data collection
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


def info9(start_, end_, freq_, cp_file_name, ts_file_name):
    """
    Information about the CP file and it's associated TS files
      that it was derived from processing information.

    :param start_: when collection of data for the file began
    :param end_: when collection of data for the file ended
    :param freq_: sampling rate in hz
    :param cp_file_name: cross power file name
    :param ts_file_name: all time series that apply to to this cp file
    :return:
    Block string
    """

    info = ''' SamplingRate units="Hz":''' + str(freq_) + '''
 Start=''' + str(start_) + '''
 End=''' + str(end_) + '''

 ProcessingInfo:
 SignConvention=exp(+iwt)
 RemoteRef type=none

 ProcessedBy=Jay A. Sampson
 ProcessingSoftware Name=MTR02B.HPB
 ProcessingSoftware LastMod=1986-06-06
 ProcessingSoftware Author=Victor F. Labson, USGS, Denver, CO, USA
 ProcessingTag=''' + str(cp_file_name) + '''

 ProcessingTimeSeriesUsed:
''' + str(ts_file_name) + '''

 =====================================================================\n\n'''
    return info


def info10():
    """
    USGS MT Truck System Cross Power Format Description.

    :return:
    Block string
    """

    info = ''' USGS MT Truck System Cross Power Format Description:

 Exactly (NCHAN x NCHAN) data values are included in the data set. The
 data set represents an estimate of the auto and cross power spectra
 for a set of NCHAN measurements over a particular frequency range.
 Crosspower matrices are Hermitian and can be packed as shown in the
 following example. The complex auto and cross power spectra for 4
 channels A, B, C, and D:

    A   B   C   D
 A <AA*> <AB*> <AC*> <AD*>
 B <BA*> <BB*> <BC*> <BD*>
 C <CA*> <CB*> <CC*> <CD*>
 D <DA*> <DB*> <DC*> <DD*>

 can be stored as follows:

    A   B   C   D
 A <AA*> Imag<AB*> Imag<AC*> Imag<AD*>
 B Real<AB*> <BB*> Imag<BC*> Imag<BD*>
 C Real<AC*> Real<BC*> <CC*> Imag<CD*>
 D Real<AD*> Real<BD*> Real<CD*> <DD*>

 with the real part in the lower left triangle and the imaginary part
 in the upper right triangle. Note that the auto spectra are real.  The
 only ambiguity is the sign of the imaginary values. This is defined
 such that the sign associated with the imaginary values is in the upper
 right triangle of the original matrix is preserved in the compressed
 matrix. The data set is this compressed matrix, read by row, e.g.
 <AA*>, Imag<AB*>,  Imag<AC*>, ...  Real<CD*>, <DD*>.

 The data appear in the >=SPECTRASECT block below. The value of the FREQ
 option gives the center frequency of the spectra estimates in Hz. The
 BW option specifies the bandwidth in Hz between the (half-power) cut-off
 frequencies.

 The value of the AVGT option is the number of independent estimates in
 time which were averaged to make these spectra estimates. Cascade
 decimation averages constant percentage bandwidth spectra estimates
 which are independent in time. Weighted averaging can lead to a number
 of independent samples which is not an integer. The AVGF option is the
 number of independent estimates in frequency (number of harmonics)
 which were averaged to make these spectra estimates. Spectra estimates
 generated from an FFT are averaged in frequency to produce constant
 percentage bandwidth spectra. Note that if a number of FFT runs are
 averaged, both AVGT and AVGF can be greater than 1.

 ===================================================================== \n\n'''
    return info


def dd_to_dms(dd_):
    """
    Changes geospatial values from decimal degrees to degree minute seconds.

    :param dd_: input of either latitude or longitude
    :return:
    """

    is_positive = dd_ >= 0
    dd = abs(dd_)
    minutes, seconds = divmod(dd * 3600, 60)
    degrees, minutes = divmod(minutes, 60)
    degrees = degrees if is_positive else -degrees
    return degrees, minutes, seconds


def s_notation(num_, dec_):
    """
    Takes a float number and changes it into scientific notation based on the
      number of decimal places required.  dec_ = 5 or 6 have specific
      requirements as part of a consistent formatted output

    :param num_: number to be put into scientific notation
    :param dec_: number of decimal places desired
    :return:
    """
    ns = ''
    if float(num_) < 0:
        num = ('%.' + str(dec_) + 'E') % Decimal(num_)
        if dec_ == 6:
            ns = '' + str(num)
        elif dec_ == 5:
            ns = ' ' + str(num)
        return ns

    else:
        num = ('%.' + str(dec_) + 'E') % Decimal(num_)

        if dec_ == 6:
            ns = ' ' + str(num)
        elif dec_ == 5:
            ns = '  ' + str(num)
        return ns


def cmat(cp_values, nf_, nch_):
    """
    Creates an array from the ascii file using imaginary and real values.
      This array is specific to the 1990 software that was run on the USGS
      truck mounted system. This is how to unpack this ascii file into an
      array that is formatted appropriately for the EDI format.

    :param cp_values: a list of values from the cross power ascii file
    :param nf_:
    :param nch_:
    :return:
    """
    null, z = '100000000000000000000000000000000', 0
    nf2 = nf_ * 2
    cp_array_list = []
    for k in range(1, nf_ + 1):
        a = np.empty([2, 7, 7], dtype=object)
        p = k * 2 - 2
        for i in range(nch_):
            for j in range(i + 1):
                try:
                    a[0, j, i] = cp_values[p]
                    a[0, i, j] = cp_values[p]

                    if i == j:
                        p = p + nf2

                    else:
                        z = cp_values[p + 1]
                        z = float(z)
                        z = z * -1.0
                        neg_z = str(z)
                        a[1, j, i] = cp_values[p + 1]
                        a[1, i, j] = neg_z
                        p = p + nf2

                except ValueError:
                    print('Something went wrong at the creation of the Array.')

        # Important to know:  The XP spectra mats.xls gives us the appropriate
        #   array to include in the CP edi files. The values in that array are
        #   reversed here below and the values in python count from 0-5 where
        #   the spectra mats file counts from 1-6.  ...append(a[0,0,0]) first
        #   value is real/imaginary, second is column, and the third is the
        #   row.

        # Row 1
        cp_array_list.append(a[0, 0, 0])  # Column 1
        cp_array_list.append(a[1, 3, 0])  # Column 2
        cp_array_list.append(null)        # Column 3
        cp_array_list.append(a[1, 5, 0])  # Column 4
        cp_array_list.append(a[1, 2, 0])  # Column 5
        cp_array_list.append(a[1, 1, 0])  # Column 6
        cp_array_list.append(a[1, 4, 0])  # Column 7

        # Row 2
        cp_array_list.append(a[0, 3, 0])
        cp_array_list.append(a[0, 3, 3])
        cp_array_list.append(null)
        cp_array_list.append(a[1, 5, 3])
        cp_array_list.append(a[1, 2, 3])
        cp_array_list.append(a[1, 1, 3])
        cp_array_list.append(a[1, 4, 3])

        # Row 3
        cp_array_list.append(null)
        cp_array_list.append(null)
        cp_array_list.append(null)
        cp_array_list.append(null)
        cp_array_list.append(null)
        cp_array_list.append(null)
        cp_array_list.append(null)

        # Row 4
        cp_array_list.append(a[0, 5, 0])
        cp_array_list.append(a[0, 5, 3])
        cp_array_list.append(null)
        cp_array_list.append(a[0, 5, 5])
        cp_array_list.append(a[1, 2, 5])
        cp_array_list.append(a[1, 1, 5])
        cp_array_list.append(a[1, 4, 5])

        # Row 5
        cp_array_list.append(a[0, 2, 0])
        cp_array_list.append(a[0, 3, 2])
        cp_array_list.append(null)
        cp_array_list.append(a[0, 5, 2])
        cp_array_list.append(a[0, 2, 2])
        cp_array_list.append(a[1, 1, 2])
        cp_array_list.append(a[1, 4, 2])

        # Row 6
        cp_array_list.append(a[0, 1, 0])
        cp_array_list.append(a[0, 3, 1])
        cp_array_list.append(null)
        cp_array_list.append(a[0, 5, 1])
        cp_array_list.append(a[0, 2, 1])
        cp_array_list.append(a[0, 1, 1])
        cp_array_list.append(a[1, 4, 1])

        # Row 7
        cp_array_list.append(a[0, 4, 0])
        cp_array_list.append(a[0, 4, 3])
        cp_array_list.append(null)
        cp_array_list.append(a[0, 5, 4])
        cp_array_list.append(a[0, 4, 2])
        cp_array_list.append(a[0, 4, 1])
        cp_array_list.append(a[0, 4, 4])
    return cp_array_list


def get_ts_list(ts_file):
    """
    Takes the time series file as the argument and grabs each line as
      a part of a new list that is used to create the cross power EDI.
      (grabs one of the corresponding TS file for essential data)

    :param ts_file: time series file
    :return:
    """

    ts_list = []
    for ascii_value in ts_file:
        ts_list.append(ascii_value.strip())
    return ts_list


# ---------------------------End of Functions---------------------------------
# ---------------Beginning of variables and calling functions-----------------

# Sets variables that don't change based on the different CP files inserted.

def main(ts, cp, edi, qc_dir, locations):
    """
    Uses the 5 dir/files to run all the appropriate functions to output
      all the necessary information for the Cross power EDI files.

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

    if 'CP_edi' not in os.listdir(edi):
        os.mkdir(edi + '\\CP_edi')

    sites = os.listdir(cp)

    # Loops through all of the site folders
    for s in sites:
        if s not in os.listdir(edi + '\\CP_edi'):
            os.mkdir(edi + '\\CP_edi\\' + s)

        if s[-2:] == '11':
            site = '10'

        # Define variables for time series(TS) and cross power(CP) sites.
        location = 0                                 # Keeps track of file count.
        crossPower = cp + '\\' + s                   # CP site folder.
        timeSeries = ts + '\\' + s                   # TS site folder.
        cpFileNames = os.listdir(crossPower)         # List of all the CP files.
        countL, countD, countX, countK = 0, 0, 0, 0  # Band counts.
        addTimeD = 0                                 # Tracks count of TS files.

        print(cpFileNames)

        # Now the script runs through each CP file and creates a list within
        #   python that will store the file line by line into a list we can
        #   pull values from.
        for n in cpFileNames:
            # Define variables used to build CP specific output edi file name.
            types = '-CP'                                # CP name.
            count = 0                                    # Rows in CP ascii file.
            cpFile = open((crossPower + '\\' + n), 'r')  # CP ascii file.
            cpList = []                                  # Rows in CP file.
            tsList = []                                  # Rows in TS file.
            name = n

            for row in cpFile:
                cpList.append(row.strip())
                count += 1
            cpFile.close()

            # Checks band for the CP file and opens associated TS file for it (D,
            #   K, L, and X).  Adds the Quality Control values associated with
            #   the file band (from the qcDict created from the Qc files) and
            #   run number.

            # Extracts site, band, and run from the cross power file name.
            cpSite, cpBand, cpRun = n[4:6], n[3], n[-1]
            # Checks for the Medium band 'D'
            if cpBand == 'D':
                countD += 1
                band, tsFileNames, tsStacked, plus = 'Mid', '', int(cpList[2]), ''
                rating, goodFrom, goodTo, comments, tsCount = \
                    qcDict[n][0], qcDict[n][1], qcDict[n][2], qcDict[n][3], \
                    countD * tsStacked
                lpFilter, hpFilter = '2.0', '0.1'

                # Adjusts for T.S and .TS endings to open the correct TS file.
                if countD == 1:
                    tsFile, ext = open(ts + '\\' + s + '\\' + 'PARD' + cpSite
                                       + str(int(cpRun) * tsStacked) + 'T.S',
                                       'r'), 'T.S'
                else:
                    tsFile, ext = open(ts + '\\' + s + '\\' + 'PARD' + cpSite
                                       + str(int(cpRun) * tsStacked) + '.TS',
                                       'r'), '.TS'

                # Loops through the number of time series stacked within the
                #   CP file and creates that list of TS files associated
                #   with this current CP file (ex: for this MT data, TS stacked
                #   is 8 for every CP file).
                for x in range(1, tsStacked + 1):
                    if (x + ((countD - 1) * tsStacked)) < 10:
                        ext = 'T.S'
                    else:
                        ext = '.TS'
                    plus = ' PARD' + cpSite + str(x + ((countD - 1) * tsStacked)) \
                           + ext + nl
                    tsFileNames = tsFileNames + plus
                tsList = get_ts_list(tsFile)
                tsFile.close()

            # Checks for the High band 'K'
            elif cpBand == 'K':
                countK += 1
                band = 'High'
                rating, goodFrom, goodTo, comments = \
                    qcDict[n][0], qcDict[n][1], qcDict[n][2], qcDict[n][3]
                lpFilter, hpFilter = '202.0', '1.0'
                tsFileNames = ' N/A'

            # Checks for the Low band 'L'
            elif cpBand == 'L':
                countL += 1
                band, tsCount = 'Low', ((countD * tsStacked) + countL)
                rating, goodFrom, goodTo, comments = \
                    qcDict[n][0], qcDict[n][1], qcDict[n][2], qcDict[n][3]
                lpFilter, hpFilter = '0.202', '0.003'
                tsFile = open(ts + '\\' + s + '\\' + 'PARL' + cpSite
                              + str(tsCount) + '.TS', 'r')
                tsFileNames = ' PARL' + cpSite + str(tsCount) + '.TS'
                tsList = get_ts_list(tsFile)
                tsFile.close()

            # Checks for the Lower band 'X'
            elif cpBand == 'X':
                countX += 1
                band, tsCount = 'xLow', ((countD * tsStacked) + countL + countX)
                rating, goodFrom, goodTo, comments = \
                    qcDict[n][0], qcDict[n][1], qcDict[n][2], qcDict[n][3]
                lpFilter, hpFilter = '0.202', '0.003'
                tsFile = open(ts + '\\' + s + '\\' + 'PARX' + cpSite
                              + str(tsCount) + '.TS', 'r')
                tsFileNames = ' PARX' + cpSite + str(tsCount) + '.TS'
                tsList = get_ts_list(tsFile)
                tsFile.close()

            # Adds a 0 before the run number if it's only 1 digit.
            if len(cpRun) == 1:
                cpRunWrite = '0' + cpRun
            elif len(cpRun) == 2:
                cpRunWrite = cpRun

            # Opens the CP output file with the correct format name for the
            #   edi standard.
            ediFile = open(edi + '\\CP_edi' + '\\' + s
                           + locName + str(cpSite) + types + '_' + band
                           + '_' + cpRunWrite + '.edi', 'w')

            # Sets todays date and a start time from the CP file (cpList[0])
            todaysDate = datetime.datetime.now()
            d = datetime.datetime.strptime(cpList[0], '%H:%M:%S %d %b %Y')
            startTime = d

            # Creates an end time based on frequency and amount of data
            #   collected (cpList[2]).
            if band == 'High':
                freq = '250'
                endTime = startTime \
                    + timedelta(seconds=1024 / (4 * .01 * int(cpList[2])))
            else:
                freq = tsList[1]
                endTime = startTime \
                    + timedelta(seconds=int(int(tsList[3])) / (
                                (4 * float(tsList[1])) * int(cpList[2])))

            # Sets data ID, contractor, who created the CP file, and
            #   the date when this file is created.
            dataID = '"' + cpList[1] + '"'
            contractor, fileby = 'USGS', '"Jay Sampson"'
            acquireDate = d.strftime('%m/%d/%Y')

            # Reformat the start and end date.
            start = str(startTime.strftime('%Y-%m-%d')) + 'T' + str(
                startTime.strftime('%H:%M:%S')) + ' UTC/GMT'
            end = str(endTime.strftime('%Y-%m-%d')) + 'T' + str(
                endTime.strftime('%H:%M:%S')) + ' UTC/GMT'

            # Location variables.
            country = 'USA'
            state = 'California'
            if cpSite == '09':
                county = 'Fresno'
            else:
                county = 'Monterey'

            # Conversion of Spatial Data variables (lat long) from
            #   decimal degrees(DD) to degrees minutes seconds (DMS).
            units = 'M'
            elevation = LatLong[2]['Elevation']['site' + str(cpSite)]
            latDD = LatLong[1]['Latitude']['site' + str(cpSite)]
            latDMS = dd_to_dms(latDD)
            latDMSf = str(int(latDMS[0])) + ':' + str(int(latDMS[1])) \
                + ':' + str(int(latDMS[2]))
            longDD = LatLong[0]['Longitude']['site' + str(cpSite)]
            longDMS = dd_to_dms(longDD)
            longDMSf = str(int(longDMS[0])) + ':' + str(int(longDMS[1])) \
                + ':' + str(int(longDMS[2]))

            rotSpec = '0'
            avgT = cpList[2]

            # Number of channels and frequencies.
            nch = int(cpList[33])
            nf = cpList[44]

            # Calculates length of dipole along EX and EY.
            dipoleLengthEX = int(1E6 * float(cpList[5]))
            dipoleLengthEY = int(1E6 * float(cpList[8]))

            # Using pyproj module project lat/long values into UTM zone 10s.
            myProj = Proj("+proj=utm +zone=10S, +north +ellps=clrk66 "
                          + "+datum=NAD27 +units=m +no_defs")
            lon, lat = myProj(longDD, latDD)

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

            # Creates the bw list and avg frequency list that is added in
            #  before each array section in the SPECTRASECT block.
            bwList, avgfList = [], []
            for x in cpList[80:(80 + int(nf))]:
                bwList.append(x)
            for x in cpList[115:(115 + int(nf))]:
                avgfList.append(x)

            # Labelling magnetometer id's for each site.
            if str(cpSite) in ['01', '02', '03', '04', '05', '11']:
                hx1, hx2, hy1, hy2 = 1, 5, 7, 9
            elif str(cpSite) in ['06']:
                hx1, hx2, hy1, hy2 = 1, 9, 5, 7
            elif str(cpSite) in ['07']:
                hx1, hx2, hy1, hy2 = 1, 5, 7, 4
            elif str(cpSite) in ['08', '09']:
                hx1, hx2, hy1, hy2 = 1, 4, 5, 7
            else:
                hx1, hx2, hy1, hy2 = 1, 5, 7, 9

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
            ediFile.write(info1(cpSite, types, band, cpRunWrite, sbUrl, n))
            ediFile.write(info2())
            ediFile.write(info3(todaysDate.strftime('%Y/%m/%d %H:%M:%S ')
                                + 'UTC-0700', conversionFile))
            ediFile.write(info4() + info5())
            ediFile.write(info6(name, latDD, longDD, elevation, round(lat, 2),
                                round(lon, 2), start, end, rating, goodFrom,
                                goodTo, comments))
            ediFile.write(info7(magnetometer, azm1, azm2, azmEx, azmEy,
                                manufacturer1, manufacturer2, dipoleLengthEX,
                                dipoleLengthEY, hx1, hx2, hy1, hy2))
            ediFile.write(info8())
            ediFile.write(info9(start, end, freq, name, tsFileNames))
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

            if band == 'High':
                gainList = ['N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', ]

            if len(cpSite) == 2 and cpSite[0] == '0':
                cpSite = cpSite[1]

            # Writes out the >HMEAS lines for each channel to output edi file.
            ediFile.write('>HMEAS ID=  ' + str(cpSite) + str(1) + '.' + '001'
                          + ' ' + 'CHTYPE=HX  ' + 'X=  0 Y=  0 AZM=  ' + azm1
                          + ' ' + nl)
            ediFile.write('>HMEAS ID=  ' + str(cpSite) + str(2) + '.' + '001'
                          + ' ' + 'CHTYPE=HY  ' + 'X=  0 Y=  0 AZM=  ' + azm2
                          + ' ' + nl)
            ediFile.write('>HMEAS ID=  ' + str(cpSite) + str(3) + '.' + '001'
                          + ' ' + 'CHTYPE=HZ  ' + 'X=  0 Y=  0 AZM=  0' + nl)
            ediFile.write('>EMEAS ID=  ' + str(cpSite) + str(4) + '.' + '001'
                          + ' ' + 'CHTYPE=EX  ' + 'X=  0 Y=  0 X2=  ' + str(exX2)
                          + ' Y2= ' + str(exY2) + nl)
            ediFile.write('>EMEAS ID=  ' + str(cpSite) + str(5) + '.' + '001'
                          + ' ' + 'CHTYPE=EY  ' + 'X=  0 Y=  0 X2=  ' + str(eyX2)
                          + ' Y2= ' + str(eyY2) + nl)
            ediFile.write('>HMEAS ID=  ' + str(cpSite) + str(6) + '.' + '001'
                          + ' ' + 'CHTYPE=HX  ' + 'X=  0 Y=  0 AZM=  ' + azm1
                          + ' ' + nl)
            ediFile.write('>HMEAS ID=  ' + str(cpSite) + str(7) + '.' + '001'
                          + ' ' + 'CHTYPE=HY  ' + 'X=  0 Y=  0 AZM=  ' + azm2
                          + ' ' + nl + nl)

            # Creating data values in arrays.
            #   >=SPECTRASECT EDI section.
            ediFile.write('>=SPECTRASECT' + nl)
            ediFile.write(indent + 'SECTID= ' + cpSite + nl)
            ediFile.write(indent + 'NCHAN= ' + str(nch + 1) + nl)
            ediFile.write(indent + 'NFREQ= ' + nf + nl)
            ediFile.write(indent + 'MAXBLKS=100' + nl + nl)
            ediFile.write(indent + '//' + '7' + nl)
            ediFile.write(indent + str(cpSite) + '1.001  ' + str(cpSite)
                          + '2.001  ' + str(cpSite) + '3.001  ' + str(cpSite)
                          + '4.001  ' + str(cpSite) + '5.001  ' + str(cpSite)
                          + '6.001  ' + str(cpSite) + '7.001  ' + nl)

            # Looping through the frequencies in the Cross Powers Ascii file
            #   and putting them in a format for the EDI and adding them to a
            #   list to pull from in the spectra section.
            add, end, freqList = 1, 0, []
            for x in cpList[45:(45 + int(nf))]:
                freqList.append(s_notation(x, 6))

            CpArrayList = cmat(cpList[165:-1], int(nf), int(nch))
            n = 0
            for freq in reversed(freqList):
                n += 1
                ediFile.write(nl + '>SPECTRA FREQ=' + s_notation(freq, 5)[1:]
                              + ' ROTSPEC= ' + rotSpec + ' BW='
                              + s_notation(bwList[freqList.index(freq)], 5)[1:]
                              + ' AVGT= ' + avgT + ' AVGF= '
                              + avgfList[freqList.index(freq)]
                              + ' //' + str((int(nch) + 1) * (int(nch) + 1)) + nl)
                add, end = 0, 0

                # Outputs the arrays and their parameters to meet EDI standard.
                for x in CpArrayList[(int(nf) - n) * 49: 49 * ((int(nf) + 1) - n)]:
                    if x is None:
                        x = '100000000000000000000000000000000'
                    x = s_notation(str(x), 5)
                    x = str(x)
                    end += 1
                    if end == 1:
                        ediFile.write('')
                        col = 9 - len(x)
                        ediFile.write(str(x) + (' ' * col))
                    elif end < 7:
                        col = 9 - len(x)
                        ediFile.write(str(x) + (' ' * col))
                    elif end == 7:
                        col = 9 - len(x)
                        ediFile.write(str(x) + (' ' * col) + nl)
                        end = 0
                        add += 1

            ediFile.write(nl + nl + '>END')
            ediFile.close()

if __name__ == '__main__':
    for item in [ts_, cp_, edi_, qcDir_, locations_]:
        if item is None:
            print ('Error: One of the inputs is None.')
            break
    main(ts_, cp_, edi_, qcDir_, locations_)
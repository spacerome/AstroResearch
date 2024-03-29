#!/usr/bin/env python

# Simulation Script for DL ASTRO

import astropy as ap
import numpy as np
from astropy.io import fits
import math as m
from astropy.stats import sigma_clip
from astropy.stats import sigma_clipped_stats
from types import *
import copy as copypkg
import os, shutil, glob, sys, string, re, math, operator, random
import glob, argparse, random

#####
#Gather List of Files
#####

###############################################################################

#####
#Keywords
#####

# def keywords():
#     return None
# observe = 'obs'
# truesky = 'tru'
# cray = 'crm'
# lacray = 'lac'
# crfrackey = 'CRFRAC'
# crvext = 'crv'
# maxpix = 65565
# blocksize = 64
# obsext = '_obs'
# truext = '_tru'
# crmext = '_crm'
# lacext = '_lac'
# crvext = '_crv'

###############################################################################
#
# Dictionary of Keywords
#
###############################################################################

keywords = {

    'observe'   :   'obs',
    'crmask'    :   'crm',
    'lacosmic'  :   'lac',
    'crayval'   :   'crv',
    'maxpix'    :   65565,
    'blocksize' :   64,
    'obsext'    :   '_obs',
    'truext'    :   '_tru',
    'crmext'    :   '_crm',
    'lacext'    :   '_lac',
    'crvext'    :   '_crv'

}

###############################################################################
#
# Modifying This May Combine Later
#
###############################################################################

def check_exist(filename, status):

    """
    check_exist(filename, status)
    checks to see if filename exists
    if status==r, must exist, otherwise prints error + returns False
    if status==w, if exists and clobber=no then prints error + returns False
    else deletes + returns True
    """

    if (status == "r"):
        # check to see if it exists for reading
        # (i.e. must be present)
        if (not (os.path.exists(filename))):
            print(f"Couldn't open input file: {filename}.")
            return False
    else:
        # check to see if it exists for writing
        # (i.e. must not exist or clobber=yes)
        if (os.path.exists(filename)):
            if (status == "w"):
                return True
            else:
                return False

    return True

###############################################################################

def issorted(ims):

    fit = sorted(ims, key=lambda item: (int(item.partition(' ')[0])
                                  if item[0].isdigit()
                                  else float('inf'), item))

    return fit

###############################################################################

def openfits(fitsIm):

    # Open Given Fits File Assuming the calling function is calling the name
    # from list

    # try:
    #     check_exist(fitsIm)
    # except:
    #     pass # Will Change Later

    fitsIm = fits.open(fitsName)

    return fitsIm

###############################################################################

def closefits(fitsIm):

    #Closes Fits File Image Name
    fitsIm.close()

    return None

###############################################################################

def fileList(img):

    img = issorted(img) # Sort Image List
    for i in range(len(img)):

        break # Will remove later

    return None #Will change later

###############################################################################

def randomImage(noutput):

    choice = random.randrange(noutput)

    return choice

###############################################################################

###############################################################################
#
# crm to crv
#
###############################################################################

def crm2crv(fileName):
    fileIndex = fileName.strip('_obs.fits')
    crmIndex = fileIndex + '_crm.fits'
    crvIndex = fileIndex + '_crv.fits'
    truIndex = fileIndex + '_tru.fits'
    lacIndex = fileIndex + '_lac.fits'
    obsFile = fits.open(fileName)
    crmFile = fits.open(crmIndex)
    truFile = fits.open(truIndex)
    obsData = obsFile.data
    crmData = crmFile.data
    truData = truFile.data
    crayVal = obsData - truData
    crvData = numpy.where(crmData ==1)
    # RAW - SKY (Gets Cray Values)
    # Np Where cdata1 = 1 (Cray)
    # if np where = 1 then we input the cray value
    # extra condition: IF value is -1 then we multiply it by -1

    return None

###############################################################################

def crmv():

    # Master loop to check number of outputs reached and then
    # Pick random seed image

    return None

###############################################################################

def usage():
    (xdir,xname) = os.path.split(sys.argv[0])
    print(f"Usage: {xname} [-h] [-b blocksize] [-n noutput] [-f crfrac]\n")

###############################################################################
#
# Updated Version of Function from iqutils
#
###############################################################################

def shiftImage(input, output, shift, border = 0, crmkey = "CRM", crmnew = "",
               skysec = "SKYSEC", clobber = globclob, verbose = globver):

    '''

        shifts input image.  WCS is preserved.

        input   name of input image
        output  name of output image
        shift   the shift:  [dx,dy] in pixels
        border  value to set in border regions [0]
        bpmkey  bad pixel mask header keyword [BPM]
        bpmnew  new bad pixel mask for shifted image [none]

        verbose print messages about actions [yes]

    '''

    # Input checking

    if not os.path.exists(input) and not os.path.exists(input + ".fits"):
        print(f"Unable to open input image {input}.")
        return None
    elif os.path.exists(input + ".fits"):
        input += ".fits"
    if input == output:
        print("Unable to shift to same filename.")
        return None
    check_exist(output, 'w')

    # Open the input image as objects

    fimg        =   fits.open(input)
    hdr         =   fimg[0].header
    d           =   fimg[0].data
    (iny,inx)   =   D.shape

    # SHIFT

    dx,dy       =   int(shift[0]), int(shift[1])
    szx,szy     =   inx - abs(dx), iny - abs(dy)

    # Coordinate ranges for zero/positive shifts

    x0,y0       =   0, 0
    newx0,newy0 =   dx, dy

    # Change for negative shifts

    if dx < 0:
        x0 =- dx
        newx0 = 0
    if dy < 0:
        y0 =- dy
        newy0 = 0

    # Make the output data as Numarray

    DX = 0 * D + border

    # Insert Data from input image

    DX[newy0:newy0 + szy, newx0:newx0 + szx] = D[y0:y0 + szy, x0:x0 + szx]

    fimg[0].data = DX

    hdr.update("SHIFTED", f"Shifted from {input}")

    # Correct WCS (CRPIX) settings


    return None

###############################################################################
#
# WILL CHANGE LATER AND USE ASTOPY
#
###############################################################################

def rotate():

    return None

###############################################################################
#
# Main script
#
###############################################################################

def main():

    # Parse Command line

    try:
        opts, args = getopt.getopt(sys.argv[1:],
                     "hdb:n:f:",
                    ["help", "debug", "block = ", "nout = ", "crfrac = "])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    # Defaults

    debug   =   0
    noutput =   1
    block   =   64
    crfrac  =   0.10

    # Process details (sometimes useful)

    (xdir,xname)    = os.path.split(sys.argv[0])
    pid             = os.getpid()

    # Options parsing

    for opt, val in opts:
        # Usage info only
        if opt in ("-h", "--help"):
            usage()
            sys.exit(1)
        # Debug mode
        elif opt in ("-d", "--debug"):
            debug = 1
        # Blocksize (in pixels)
        elif opt in ("-b","--block"):
            block = val
        # Number of output files
        elif opt in ("-n", "--nout"):
            noutput = val
        # Target CR fraction
        elif opt in ("-f", "--crfrac"):
            crfrac = val
        else:
            # sys.stderr.write("Unmatched option %s\n" % opt)
            sys.stderr.write(f"Unmatched option {opt}\n")

    # Arguments are same-size input files to use in construction

    if len(args) < 1:
        usage()
        sys.exit(2)

    infiles =   args
    nfiles  =   len(infiles)

    # read initial image input files and will gather crms

    # store initial crvs within a list

    '''

    Using CRAY class to store data? Not sure if this is a good Idea as of yet

    '''

    '''

    crmv=crmfiles
    crv=crvfiles
    MIGHT MAKE THIS EASIER TO PREVENT FURTHER CONFUSION
    Use resultant name say
    fitsfile(listoffitsfiles)

    '''
    fitsFile = glob.glob('*.fits')

    obs = glob.glob('*_obs.fits')
    tru = glob.glob('*_tru.fits')
    crm = glob.glob('*_crm.fits')
    lac = glob.glob('*_lac.fits')

    # Sort Images

    obs = issorted(obs)
    tru = issorted(tru)
    crm = issorted(crm)
    lac = issorted(lac)

    # Makes crvfiles
    for i in range(len(obs)):
        crm2crv(obs)
    ###########################################################################

    for i in range(noutput):

        # pick base image

        n = randomImage(nfiles) # calls function to generate random value
        imageName = obs[n] # Chooses name of fits Image

        # Indexing
        fIndex = imageName.strip('_obs.fits')
        crmName = fileIndex + '_crm.fits'
        crvName = fileIndex + '_crv.fits'
        truName = fileIndex + '_tru.fits'
        lacName = fileIndex + '_lac.fits'
        # read in data, CRM

        crmFile = fits.open(crmName)
        crmData = crmFile.data


        # calculate initial XFRAC



        while (xfrac < crfrac):

            # pick add-on image



            # choose reflection / rotation



            # choose shift



            # add cosmic rays according to selected reflection, rotation, and
            # shift



            # augment CR mask appropriately



            # calculate new XFRAC



            # (include some logic to avoid infinite loop?)



        # define name for the output file

        # obsnew = ''

        # write the output image file

        '''

        use following class to output fitsfile to reduce excess data
        loss using FitsFileImages() class for this process

        '''

        # write the output CRM, TRU, ...

        '''

        use following class to output fitsfile to reduce excess data
        loss using FitsFileImages() class for this process

        '''

    #################################

    # FIND RAW CRM CRV ETC
    # Make CRV
    # Output file Clobber (Flag to overwrite)
    # Clobber can be used to overwrite
    # Write file


###############################################################################
#
# This will be modified or removed
#
###############################################################################

# def usage():
#
#     (xdir,xname)=os.path.split(sys.argv[0])
#
#     print(" Usage: %s [-n name] [-e epoch] <ra> <dec> [size]" % xname)
#     print("    <ra> and <dec> are sexagesimal hours/deg or decimal deg/deg")
#     print("    <size> is length of a side in arcmin (%s)" % def_radius)
#     print("    -n name gives the source name for output files (%s)" % \
#           def_source)
#     print("    -e epoch gives epoch, J2000 or B1950 (%s)" % def_epoch)

if __name__=="__main__":
    main()

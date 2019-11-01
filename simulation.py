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
import glob
import argparse
import random

#####
#Gather List of Files
#####

###############################################################################



###############################################################################

#  
# fitsImage = glob.glob('NGC6652_0???.fits') #Change to other file names
# fitsImage = sorted(fitsImage, key=lambda item: (int(item.partition(' ')[0])
#                                     if item[0].isdigit()
#                                     else float('inf'), item))
# 

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

# sorts resultant list

# Initial Global Variables

globclob = yes
globver = yes

###############################################################################

# Planned Class for ALL Fits Types WIP

###############################################################################

class FitsFileImages():

    def __init__(self, fname = None, fext = None, ftype = None, mpix= 65565,
                 bsize = 64, cfrac = 0.10):

        # Initial Null Properties

        self.fname = None
        self.fext =  None
        self.ftype = None

        # Non-null Properties

        self.mpix =     int(mpix)
        self.bsize =    int(bsize)
        self.cfrac =    float(cfrac)

        # Reset Properties

        if fname:   self.fname =    str(fname)
        if fext:    self.fext =     str(fext)
        if ftype:   self.ftype =    str(ftype)


###############################################################################

# Class for cosmic ray images WIP

###############################################################################

class CosmicRayImages():

    def __init__(self, fname = None, fext = None, ftype = None, mpix= 65565,
                 bsize = 64):

        # Initial Null Properties

        self.fname = None
        self.fext =  None
        self.ftype = None

        # Non-null Properties

        self.mpix =     int(mpix)
        self.bsize =    int(bsize)

        # Reset Properties

        if fname:   self.fname =    str(fname)
        if fext:    self.fext =     str(fext)
        if ftype:   self.ftype =    str(ftype)


###############################################################################
# Modifying This 

def fitsfile(file):
    outfile=""

    if file.endswith('.fits') or file.endswith('.imh'):
        if os.path.exists(file):
            outfile=file
        else:
            print("Can't find requested file '%s'" % file)
    else:
        if os.path.exists(file+'.fits'):
            outfile=file+'.fits'
        elif os.path.exists(file):
            outfile=file
        elif os.path.exists(file+'.imh'):
            outfile=file+'.imh'
        else:
            print("Can't find requested file '%s' or variants" % file)
        
    return outfile

###############################################################################
# Modifying This
def check_exist(filename, status, clobber=globclob):

    """ 

    check_exist(filename, status, clobber=yes)
    checks to see if filename exists
    if status==r, must exist, otherwise prints error + returns False
    if status==w, if exists and clobber=no then prints error + returns False
    else deletes + returns True

    """     

    if (status == "r"):
        # check to see if it exists for reading
        # (i.e. must be present)
        if (not (os.path.exists(filename))):
            print("Couldn't open input file: %s" % filename)
            return False
    else:
        # check to see if it exists for writing
        # (i.e. must not exist or clobber=yes)
        if (os.path.exists(filename)):
            if (clobber):
                os.remove(filename)
            else:
                print("File %s already exists and clobber=no" % filename)
                return False

    return True

###############################################################################

def issorted(ims):

    fit = sorted(ims, key=lambda item: (int(item.partition(' ')[0])
                                  if item[0].isdigit()
                                  else float('inf'), item))

    return fit

###############################################################################
# Will be modified Just a Place Holder
def iraffiles(files,nfiles=0):

    if type(files) is not StringType:
        print "Input filelist is not a string"
        exit(1)

    fout=[]
    fmult=files.split(",")
    for fcand in fmult:
        re1=re.search("^(.+//)?@(.+)(//.+)?$",fcand)
        re2=re.search("[\*\?]",fcand)
        if re1:
            # Using the IRAF "@file.lis" convention
            flist=re1.group(2)
            if os.path.exists(flist):
                fflist=getlines(flist)
                for fmem in fflist:
                    if re1.group(1):
                        fmem=re1.group(1)[:-2]+fmem
                    if re1.group(3):
                        fmem=fmem+re1.group(3)[2:]
                    if (fitsfile(fmem)!=""):
                        fout.append(fitsfile(fmem))
        elif re2:
            # Using UNIX wildcards
            flist=glob.glob(fcand)
            for fmem in flist:
                if (fitsfile(fmem)!=""):
                    fout.append(fitsfile(fmem))
        else:
            # Just plain filenames (?)
            if fitsfile(fcand)!="":
                fout.append(fitsfile(fcand))
            
    return fout

###############################################################################

def openfits(fitsIm):

    # Open Given Fits File Assuming the calling function is calling the name 
    # from list
    try: 
        check_exist(fitsIm)
    except:

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


    return None #Will change later

###############################################################################

def rando(fitsIm,noutput):

    return None

###############################################################################

def crm2crv():
    
####
# crm to crv
####

    # RAW - SKY (Gets Cray Values)
    # Np Where cdata1 = 1 (Cray)
    # if np where = 1 then we input the cray value
    # extra condition: IF value is -1 then we multiply it by -1 
    return None

###############################################################################

def crmv():
    #####
    # Master loop to check number of outputs reached and then
    # Pick random seed image
    ####
    return None

###############################################################################

def usage():
    (xdir,xname)=os.path.split(sys.argv[0])
    print("Usage: %s [-h] [-b blocksize] [-n noutput] [-f crfrac]\n" % \
          xname)

###############################################################################

def main():

    # Parse Command line
    try:
        opts, args = getopt.getopt(sys.argv[1:], 
                     "hdb:n:f:",
                    ["help","debug","block=","nout=","crfrac="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    # Defaults
    debug = 0
    noutput = 1
    block = 64
    crfrac = 0.10

    # Process details (sometimes useful)
    (xdir,xname)=os.path.split(sys.argv[0])
    pid=os.getpid()

    # Options parsing
    for opt, val in opts:
        # Usage info only
        if opt in ("-h", "--help"):
            usage()
            sys.exit(1)
        # Debug mode
        elif opt in ("-d", "--debug"):
            debug=1
        # Blocksize (in pixels)
        elif opt in ("-b","--block"):
            block=val
        # Number of output files
        elif opt in ("-n", "--nout"):
            noutput=val
        # Target CR fraction
        elif opt in ("-f", "--crfrac"):
            crfrac=val
        else:
            sys.stderr.write("Unmatched option %s\n" % opt)

    # Arguments are same-size input files to use in construction

    if len(args) < 1:
        usage()
        sys.exit(2)

    infiles=args
    nfiles=len(infiles)

    ###########################################################################

    for i in range(noutput):

        # pick base image

        ''' 
        rando(image,noutput) # Chooses a random image within the given number 
                             # of outputs



        '''

        # read in data, CRM



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



        # write the output image file



        # write the output CRM, TRU, ...



    #################################

    #return None

###############################################################################
# This will be modified or removed
def usage():

    (xdir,xname)=os.path.split(sys.argv[0])

    print(" Usage: %s [-n name] [-e epoch] <ra> <dec> [size]" % xname)
    print("    <ra> and <dec> are sexagesimal hours/deg or decimal deg/deg")
    print("    <size> is length of a side in arcmin (%s)" % def_radius)
    print("    -n name gives the source name for output files (%s)" % def_source)
    print("    -e epoch gives epoch, J2000 or B1950 (%s)" % def_epoch)

if __name__=="__main__":
    main()
    

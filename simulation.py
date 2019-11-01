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

        def checkexist(self,fname):
            ofile = ''

            if self.fname.endswith('.fits'):
                if os.path.exists(self.fname):
                    ofile=self.fname
                else:
                    print()

            return None


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
# Modifying This May Combine Later
###############################################################################

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
# Modifying This May Combine Later
###############################################################################

def check_exist(filename, status, clobber = globclob):

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
        print("Input filelist is not a string")
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
        pass # Will Change Later
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
#
# WILL BE CHANGING THIS BELOW MIGHT HAVE ASTROPY STUFF WILL WORK ON ASAP
#
###############################################################################

def shift_image(input,output,shift,border=0,bpmkey="BPM",bpmnew="",
                skysec="SKYSEC",clobber=globclob,verbose=globver):

    """ shifts input image.  WCS is preserved.

        input   name of input image
        output  name of output image
        shift   the shift:  [dx,dy] in pixels
        border  value to set in border regions [0]
        bpmkey  bad pixel mask header keyword [BPM]
        bpmnew  new bad pixel mask for shifted image [none]

        clobber clobber output files [yes]
        verbose print messages about actions [yes]
    """

    # Input checking
    if not os.path.exists(input) and not os.path.exists(input+".fits"):
        print "Couldn't open input image %s" % input
        return
    elif os.path.exists(input+".fits"):
        input+=".fits"
    if input==output:
        print "Can't shift to same filename, sorry"
        return
    check_exist(output,'w',clobber)

    # Open the input image as pyfits object
    fimg=fits.open(input)
    hdr=fimg[0].header
    D=fimg[0].data
    (iny,inx)=D.shape

    # The shift
    dx,dy=int(shift[0]),int(shift[1])
    szx,szy=inx-abs(dx),iny-abs(dy)

    # Coordinate ranges for zero/positive shifts
    x0,y0=0,0
    newx0,newy0=dx,dy

    # Change for negative shifts
    if dx<0:
        x0=-dx
        newx0=0
    if dy<0:
        y0=-dy
        newy0=0

    # Make the output data as a numarray
    DX=0*D+border

    # Insert data from the input image
    DX[newy0:newy0+szy,newx0:newx0+szx]=D[y0:y0+szy,x0:x0+szx]

    fimg[0].data=DX
    hdr.update('SHIFTED','Shifted from %s' % input)

    # Correct the WCS (CRPIX) settings
    if 'CRPIX1' in hdr.keys():
        crpix1=hdr['CRPIX1']
        hdr.update('CRPIX1',crpix1+dx)
    if 'CRPIX2' in hdr.keys():
        crpix2=hdr['CRPIX2']
        hdr.update('CRPIX2',crpix2+dy)

    # Reset the BPM keyword if requested
    if len(bpmkey)>0:
        if len(bpmnew)>0:
            hdr.update(bpmkey,bpmnew)
        else:
            del hdr['BPM']

    # Set or adjust the SKYSEC keyword, as appropriate
    if len(skysec)>0:
        # Default sky region = Full original image (IRAF style)
        skyreg="[%d:%d,%d:%d]" % (newx0+1,newx0+szx,newy0+1,newy0+szy)
        # Check if there was a sky region defined already
        if check_head(input,skysec):
            skyval=get_head(input,skysec)
            resky=re.search("\[(\d+):(\d+),(\d+):(\d+)\]",skyval)
            if resky:
                oldx0,oldx1,oldy0,oldy1=int(resky.group(1)), \
                                        int(resky.group(2)), \
                                        int(resky.group(3)), \
                                        int(resky.group(4))
                skyreg="[%d:%d,%d:%d]" % \
                        (oldx0+dx,oldx1+dx,
                         oldy0+dy,oldy1+dy)
        hdr.update(skysec,skyreg)

    # Write/close the output image
    fimg.writeto(output)

###############################################################################
# WILL CHANGE LATER
###############################################################################
def rotate():

    return None
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

    # read initial image input files and will gather crms

    # store initial crvs within a list

    '''
    Using CRAY class to store data? Not sure if this is a good Idea as of yet
    '''

    '''
    crv=crmfiles
    crm=crvfiles
    MIGHT MAKE THIS EASIER TO PREVENT FURTHER CONFUSION
    Use resultant name say
    fitsfile(listoffitsfiles)

    '''
    ###########################################################################

    for i in range(noutput):

        # pick base image

        n = randomImage(nfiles) # calls function to generate random value
        imageName = fitsImage[n] # Chooses name of fits Image


        # read in data, CRM

        '''
        call function that converts crm
        '''


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

        '''
        use following class to output fitsfile to reduce excess data
        loss using FitsFileImages() class for this process
        '''

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


###############################################################################
# This will be modified or removed
def usage():

    (xdir,xname)=os.path.split(sys.argv[0])

    print(" Usage: %s [-n name] [-e epoch] <ra> <dec> [size]" % xname)
    print("    <ra> and <dec> are sexagesimal hours/deg or decimal deg/deg")
    print("    <size> is length of a side in arcmin (%s)" % def_radius)
    print("    -n name gives the source name for output files (%s)" % \
          def_source)
    print("    -e epoch gives epoch, J2000 or B1950 (%s)" % def_epoch)

if __name__=="__main__":
    main()

#!/bin/env python
#
# $Id: geometry.py 3227 2007-06-20 16:15:52Z tc $
#
# This file is part of the Protool package
#
# (c) Jens Erik Nielsen, http://www.cmbi.kun.nl/gv/nielsen
#
from errors import *


try:
    import math
    from numpy import *
    from numpy.linalg import *
    linear_least_squares=lstsq
    inverse=inv
    Float=float
    eigenvectors=eig
except:
    from Numeric import *
    from LinearAlgebra import *

def length(vector):
    # This function returns the length of vector
    import math
    sum=0.0
    for value in vector:
        sum=sum+math.pow(value,2)
    return math.sqrt(sum)

def normalise (vector):
    # normalise vector
    newvector=[]
    l=length(vector)
    if l!=0:
        for item in vector:
            newvector.append(float(item)/l)
    else:
        newvector=vector
    return array(newvector)



class geometry:
    """Geometry class"""

    def distance(self,atom1,atom2):
        """Return the distance between two atoms"""
        return self.dist(atom1,atom2)

    #
    # ----
    #

    def dist(self,atom1,atom2):
        #
        # Calculates the distance between two atoms
        #
        vector=self.GetPosition(atom1)-self.GetPosition(atom2)
        return length(vector)

    #
    # ---
    #

    def superpose(self,reference_coords,fit_coords):
        """Superpose the fit_coords on the reference_coords"""
        import quatfit
        self.refcenter,self.fitcenter,self.rotation=quatfit.qfit(len(reference_coords),reference_coords,fit_coords)
        self.newcoords = quatfit.qtransform(len(reference_coords), fit_coords, self.refcenter, self.fitcenter, self.rotation)
        #
        # Calc rmsd
        #
        rmsd=0.0
        for i in range(len(self.newcoords)):
            disp = self.newcoords[i] - reference_coords[i]
            #disp=length(disp)*length(disp)
            dist_squared = dot(disp,disp)
            rmsd = rmsd + dist_squared
        rmsd=rmsd/float(len(self.newcoords))
        import math
        rmsd=math.sqrt(rmsd)
        return rmsd


    #
    # ----
    #
    

    def translate(self,ref,fit):
        #
        # Return the vector that translates the center of
        # gravity of ref into the center of gravity of ref
        #
        # Center of geometry of reference
        cogref = array([ 0., 0., 0.])
        for coord in ref:
            cogref = cogref + coord
        cogref = cogref/len(ref)

        # Center of geometry of fit
        cogfit = array([ 0., 0., 0.])
        for coord in fit:
            cogfit = cogfit + coord
        cogfit = cogfit/len(fit)

        # Find the vector to be added to fit to translate
        # COG to ref's COG
        tcogfit = cogref - cogfit
        return tcogfit

 
    #
    # ----
    #
    
    def center(self,ref,fit):
      #
      # Center both coordinate sets at 0,0,0
      #
      #
      # First check that we have some coordinates
      #
      if len(ref)==0 or len(fit)==0:
        raise "Either reference of fit coordinate set contained no coordinates"
      #
      # Center of coordinates of reference
      #
      cogref = array([0.,0.,0.])
      for coord in ref:
        cogref = cogref + coord
      cogref = cogref/len(ref)
      #
      # Produce new coordinates
      #
      center_ref=[]
      vect_ref=array([0.,0.,0.])-cogref
      for coord in ref:
        center_ref.append(coord+vect_ref)


      # Center of coordinates of fit
      cogfit = array([ 0., 0., 0.])
      for coord in fit:
        cogfit = cogfit + coord
      cogfit = cogfit/len(fit)
      #
      # Produce new coordinates
      #
      center_fit=[]
      vect_fit=array([0.,0.,0.])-cogfit
      for coord in fit:
        center_fit.append(coord+vect_fit)

      return center_ref,center_fit,vect_ref,vect_fit


#
# Validate the superpos function every time this module is imported
#


import copy
rv=array([[1,1,1],[0,-1,1],[1,-2,1],[1,-3,1],[-1,-3,2]])
fv=array([[1,1,0],[-1,2,0],[-2,1,0],[-3,1,0],[-3,3,1]])

import quatfit
refcenter,fitcenter,rotation=quatfit.qfit(len(rv),rv,fv)
newcoords = quatfit.qtransform(len(rv), fv, refcenter, fitcenter, rotation)
#
# Calc rmsd
#
rmsd=0.0
for i in range(len(newcoords)):
    disp = newcoords[i] - rv[i]
    dist = dot(disp,disp)
    rmsd = rmsd + dist
if rmsd>0.0000001:
    print 'Superpositioning test failed!!! in geometry.py' 
    print 'RMSD',rmsd
    print newcoords
    import sys
    sys.exit(1)
else:
    None
##    print 'superpositioning algorithm validated'

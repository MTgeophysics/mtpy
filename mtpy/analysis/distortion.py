#!/usr/bin/env python

"""
mtpy/mtpy/analysis/distortion.py

Contains functions for the determination of (galvanic) distortion of impedance
tensors.
The methods used follow Bibby et al 2005.
As it has been pointed out in that paper, there are various possibilities for
 constraining the solution, esp. in the 2D case.

Here we just implement the 'most basic' variety for the calculation of the 
distortion tensor.

Other mehtods can be implemented, but since the optimal assumtions and 
constraints depend on the application, the actual place for further functions
 is in an independent, personalised module.



    Functions:


@UofA, 2013
(LK)

"""

#=================================================================
import numpy as np

import mtpy.core.z as MTz 
import mtpy.analysis.geometry as MTge 
import mtpy.utils.exceptions as MTex
import mtpy.utils.calculator as MTcc

#reload(MTex)
#reload(MTz)
#reload(MTcc)
#reload(MTge)


#=================================================================
#Finding the distortion of a Z array. Using the phase tensor 
#(so, Z arrays are transformed into PTs first), following Bibby et al. 2005.
#

# First, try to find periods that indicate 1D. From them determine D incl.
# the g-factor by calculatiing a weighted mean. The g is assumed in order to 
# cater for the missing unknown in the system, it is here set to det(X)^0.5. 
# After that is found, the function no_distortion from the Z module can be 
# called to obtain the unperturbated regional impedance tensor.

#Second, if there are no 1D sections: Find the strike angle, then rotate the
# Z to the principal axis. In order to do that, use the rotate(-strike) method
# of the Z module. Then take the real part of the rotated Z. As in the 1D case,
# we need an assumption to get rid of the (2) unknowns:
# set det(D) = P and det(D) = T, where P,T can be chosen. Common choice is to
# set  one of P,T to an arbitrary value (e.g. 1). Then check, for which values
# of the other parameter  S^2 = T^2+4*P*X_12*X_21/det(X) > 0 holds.



def find_distortion(z_object, g = 'det', lo_dims = None):
    """
    find optimal distortion tensor from z object

    automatically determine the dimensionality over all frequencies, then find
    the appropriate distortion tensor D
    """

    z_obj = z_object

    if lo_dims is None :
        lo_dims = MTge.dimensionality(z_object = z_obj)
    try:
        if len(lo_dims) != len(z_obj.z):
            lo_dims = MTge.dimensionality(z_object = z_obj)
    except:
        pass
    
    #dictionary of values that should be no distortion in case distortion
    #cannot be calculated for that component
    dis_dict = {(0,0):1, (0,1):0, (1,0):0, (1,1):1}

    lo_dis = []
    lo_diserr = []

    if 1 in lo_dims:
        idx_1 = np.where(np.array(lo_dims) == 1)[0]

        for idx in idx_1:

            realz = np.real(z_obj.z[idx])
            imagz = np.imag(z_obj.z[idx])

            mat1 = np.matrix([[0, -1],[1, 0]])


            if g in ['01','10']:
                gr = np.abs(realz[int(g[0]),int(g[1])])
                gi = np.abs(imagz[int(g[0]),int(g[1])])
            else:
                gr = np.sqrt(np.linalg.det(realz))
                gi = np.sqrt(np.linalg.det(imagz))

            lo_dis.append(1./gr*np.dot(realz,mat1))  
            lo_dis.append(1./gi*np.dot(imagz,mat1))  

            if z_obj.zerr is not None:
                #find errors of entries for calculating weights

                lo_diserr.append(1./gr*\
                                np.array([[np.abs(z_obj.zerr[idx][0,1]),
                                           np.abs(z_obj.zerr[idx][0,0])],
                                           [np.abs(z_obj.zerr[idx][1,1]),
                                           np.abs(z_obj.zerr[idx][1,0])]])) 

                lo_diserr.append(1./gi*\
                                 np.array([[np.abs(z_obj.zerr[idx][0,1]),
                                            np.abs(z_obj.zerr[idx][0,0])],
                                            [np.abs(z_obj.zerr[idx][1,1]),
                                             np.abs(z_obj.zerr[idx][1,0])]])) 

            else:
                #otherwise go for evenly weighted average
                lo_diserr.append(np.ones((2, 2)))
                lo_diserr.append(np.ones((2, 2)))

        
        dis = np.identity(2)
        diserr = np.identity(2)
        for i in range(2):
            for j in range(2):
                try:
                    dis[i,j], dummy = np.average(np.array([k[i, j] 
                                                          for k in lo_dis]), 
                                                 weights=np.array([1./(k[i,j])**2 
                                                          for k in lo_diserr]),
                                                 returned=True)
                    diserr[i,j] = np.sqrt(1./dummy)
                    
                    #if the distortion came out as nan set it to an appropriate
                    #value 
                    if np.nan_to_num(dis[i,j]) == 0:
                        dis[i, j] = dis_dict[i, j]
                        diserr[i, j] = dis_dict[i, j]

                except ZeroDivisionError:
                    
                    print ('Could not get distortion for dis[{0}, {1}]'.format(
                           i, j)+' setting value to {0}'.format(dis_dict[i,j]))
                    dis[i, j] = dis_dict[i, j]
                    diserr[i, j] = dis_dict[i, j]*1e-6
    
        return dis, diserr

    if 2 in lo_dims:
        idx_2 = np.where(np.array(lo_dims) == 2)[0]
        #follow bibby et al. 2005 first alternative: P = 1
        P = 1

        lo_strikes = MTge.strike_angle(z_object = z_obj)
        lo_tetms = []
        lo_t = []
        lo_tetm_errs =[]

        for idx in idx_2:

            mat = z_obj.z[idx]
            ang = -lo_strikes[idx][0]
            if np.isnan(ang):
                ang = 0.

            errmat = None
            if z_obj.zerr is not None:
                errmat = z_obj.zerr[idx]
            tetm_mat, tetm_err = MTcc.rotatematrix_incl_errors(mat, 
                                                               ang, 
                                                               inmatrix_err=errmat)

            lo_tetms.append(tetm_mat)
         
            lo_tetm_errs.append(tetm_err)

            realz = np.real(tetm_mat)
            imagz = np.imag(tetm_mat)
            lo_t.append(-4*P*realz[0,1]*realz[1,0]/np.linalg.det(realz) )
            lo_t.append(-4*P*imagz[0,1]*imagz[1,0]/np.linalg.det(imagz) )

        #since there is no 'wrong' solution by a different value of T, no 
        #error is given/calculated for T !
        try:
            #just add 0.1% for avoiding numerical issues in the squareroots
            #later on
            T = np.sqrt(max(lo_t))+0.001
        except:
            T = 2


        for idx in range(len(lo_tetms)):

            realz = np.real(lo_tetms[idx])
            imagz = np.imag(lo_tetms[idx])
            errmat = lo_tetm_errs[idx]
            
            sr = np.sqrt(T**2+4*P*realz[0, 1]*realz[1, 0]/np.linalg.det(realz))
            si = np.sqrt(T**2+4*P*imagz[0, 1]*imagz[1, 0]/np.linalg.det(imagz))

            par_r = 2*realz[0, 1]/(T-sr)
            orth_r = 2*realz[1, 0]/(T+sr)
            par_i = 2*imagz[0, 1]/(T-si)
            orth_i = 2*imagz[1, 0]/(T+si)

            mat2_r = np.matrix([[0, 1./orth_r], [1./par_r, 0]])
            mat2_i = np.matrix([[0, 1./orth_i], [1./par_i ,0]])

            lo_dis.append(np.dot(realz,mat2_r))
            lo_dis.append(np.dot(imagz,mat2_i))

            if z_obj.zerr is not None:
                #find errors of entries for calculating weights
                sigma_sr = np.sqrt((-(2*P*realz[0,1]*realz[1,0]*\
                                      realz[1,1]*errmat[0,0])/\
                                      (np.linalg.det(realz)**2*sr))**2+\
                                    ((2*P*realz[0,0]*realz[1,0]*\
                                     realz[1,1]*errmat[0,1])/\
                                    (np.linalg.det(realz)**2*sr))**2+\
                                    ((2*P*realz[0,0]* realz[0,1]*\
                                      realz[1,1]*errmat[1,0])/\
                                      (np.linalg.det(realz)**2*sr))**2 +\
                                    (-(2*P*realz[0,1]* realz[1,0]*\
                                     realz[0,0]*errmat[1,1])/\
                                     (np.linalg.det(realz)**2*sr))**2)

                sigma_dr_11 = 0.5*sigma_sr
                sigma_dr_22 = 0.5*sigma_sr

                sigma_dr_12 = np.sqrt((mat2_r[0,1]/realz[0,0]*errmat[0,0])**2+\
                                      (mat2_r[0,1]/realz[1,0]*errmat[1,0])**2+\
                                      (0.5*realz[0,0]/realz[1,0]*sigma_sr)**2)
                sigma_dr_21 = np.sqrt((mat2_r[1,0]/realz[1,1]*errmat[1,1])**2+\
                                      (mat2_r[1,0]/realz[0,1]*errmat[0,1])**2+\
                                      (0.5*realz[1,1]/realz[0,1]*sigma_sr)**2)

                lo_diserr.append(np.array([[sigma_dr_11, sigma_dr_12],
                                           [sigma_dr_21, sigma_dr_22]]))

                sigma_si = np.sqrt((-(2*P*imagz[0,1]*imagz[1,0]*\
                                      imagz[1,1]*errmat[0,0])/\
                                      (np.linalg.det(imagz)**2*sr))**2+\
                                     ((2*P*imagz[0,0]*imagz[1,0]*\
                                      imagz[1,1]*errmat[0,1])/\
                                      (np.linalg.det(imagz)**2*sr))**2+\
                                     ((2*P*imagz[0,0]*imagz[0,1]*\
                                      imagz[1,1]*errmat[1,0])/\
                                      (np.linalg.det(imagz)**2*sr))**2+\
                                     (-(2*P*imagz[0,1]*imagz[1,0]*\
                                      imagz[0,0]*errmat[1,1])/\
                                      (np.linalg.det(imagz)**2*sr))**2)

                sigma_di_11 = 0.5*sigma_si
                sigma_di_22 = 0.5*sigma_si
                sigma_di_12 = np.sqrt((mat2_i[0,1]/imagz[0,0]*errmat[0,0])**2+\
                                      (mat2_i[0,1]/imagz[1,0]*errmat[1,0])**2+\
                                      (0.5*imagz[0,0]/imagz[1,0]*sigma_si)**2)
                sigma_di_21 = np.sqrt((mat2_i[1,0]/imagz[1,1]*errmat[1,1])**2+\
                                      (mat2_i[1,0]/imagz[0,1]*errmat[0,1])**2+\
                                      (0.5*imagz[1,1]/imagz[0,1]*sigma_si)**2)

                lo_diserr.append(np.array([[sigma_di_11, sigma_di_12],
                                           [sigma_di_21, sigma_di_22]]))

            else:
                #otherwise go for evenly weighted average
                lo_diserr.append(np.ones((2, 2)))
                lo_diserr.append(np.ones((2, 2)))


        dis = np.zeros((2, 2))
        diserr = np.zeros((2, 2))
        for i in range(2):
            for j in range(2):

                dis[i, j], dummy = np.average(np.array([k[i, j] 
                                                       for k in lo_dis]), 
                                              weights=np.array([1./(k[i,j])**2
                                                       for k in lo_diserr]),
                                              returned=True )
                diserr[i, j] = np.sqrt(1./dummy)

        return dis, diserr

    #if only 3D, use identity matrix - no distortion calculated
    dis = np.identity(2)
    diserr = diserr = np.zeros((2, 2))
    return dis, diserr



def find_1d_distortion(z_object, include_non1d = False):
    """
    find 1D distortion tensor from z object

    ONly use the 1D part of the Z to determine D. 
    Treat all frequencies as 1D, if  "include_non1d = True".

    
    """

    if not isinstance(z_object, MTz.Z):
        raise MTex.MTpyError_inputarguments('first argument must be an ' 
                                            'instance of the Z class')

    z_obj = z_object

    lo_dims = MTge.dimensionality(z_object=z_obj)

    if include_non1d is True:
        lo_dims = [1 for i in lo_dims]

    if len(list(np.where(np.array(lo_dims) == 1))) == 0:
        raise MTex.MTpyError_inputarguments('Z object does not have '
                                  'frequencies with spatial 1D characteristic')

    print lo_dims

    return  find_distortion(z_obj, lo_dims = lo_dims)



def find_2d_distortion(z_object, include_non2d=False):
    """
    find 2D distortion tensor from z object

    ONly use the 2D part of the Z to determine D. 
    Treat all frequencies as 2D, if  "include_non2d = True".
    
    """

    if not isinstance(z_object, MTz.Z):
        raise MTex.MTpyError_inputarguments('first argument must be an '
                                            'instance of the Z class')

    z_obj = z_object

    lo_dims = MTge.dimensionality(z_object = z_obj)

    #avoid the (standard) 1D distortion call -> remove all 1
    lo_dims = [ 4 if i == 1 else i for i in lo_dims ]

    if include_non2d is True:
        lo_dims = [2 for i in lo_dims]

    if len(list(np.where(np.array(lo_dims) == 2))) == 0:
        raise MTex.MTpyError_inputarguments('Z object does not have'
                                 ' frequencies with spatial 2D characteristic')


    return  find_distortion(z_obj, lo_dims = lo_dims)

def remove_distortion(z_array=None, z_object=None):

    if z_array is not None:
        z_obj = MTz.Z(z_array=z_array)
    
    elif z_object is not None:
        z_obj = z_object

    #0. generate a Z object
    #1. find distortion via function above, 
    #2. remove distortion via method of z object

    dis, diserr = find_distortion(z_obj)

    try:
        distortion_tensor, zd, zd_err = z_obj.no_distortion(dis, 
                                                distortion_err_tensor=diserr)
        distortion_z_obj = z_obj
        distortion_z_obj.z = zd
        distortion_z_obj.zerr = zd_err 

        return distortion_tensor, distortion_z_obj
        
    except MTex.MTpyError_Z:
        print 'Could not compute distortion tensor'
        
        return np.identity(2), z_obj
                                    






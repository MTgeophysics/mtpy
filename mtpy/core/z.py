#!/usr/bin/env python

"""
mtpy/mtpy/core/z.py

Contains classes and functions for handling impedance tensors (Z). 
 
    Class:
    "Z" contains information about an impedance tensor Z. 

        Methods:
        - set_edi_object
        - set_z
        - set_zerr
        - real
        - set_real
        - imag
        - set_imag
        - rho
        - phi
        - set_rho_phi
        - inverse
        - rotate
        - no_ss
        - no_distortion
        - no_ss_no_distortion
        - ellipticity
        - resistivity
        - invariants


    Class:
    "Tipper" contains information about the (complex valued) Tipper vector. 

        Methods:

        - set_edi_object
        - set_tipper
        - set_tippererr
        - real
        - set_real
        - imag
        - set_imag
        - rho
        - phi
        - set_rho_phi
        - rotate

    Functions:

     - rotate_z
     - remove_distortion
     - remove_ss
     - remove_ss_and_distortion
     - z2rhophi
     - rotate_tipper
     - _read_z_array
     - _read_tipper_array

@UofA, 2013
(LK)

"""

#=================================================================
import numpy as np
import math, cmath
import copy
import mtpy.utils.calculator as MTc


import mtpy.utils.exceptions as MTexceptions

#reload(MTexceptions)
#reload(MTformat)
#reload(MTc)


#=================================================================


#------------------------
class Z(object):
    """
        Z class - generates an impedance tensor (Z-) object.

        Methods  include reading and writing from and to edi-objects, rotations/combinations of Z instances, as well as 
        calculation of invariants, inverse, amplitude/phase,...

        
        Z is a complex array of the form (n_frequencies, 2, 2), 
        with indices in the following order: 
            Zxx: (0,0) - Zxy: (0,1) - Zyx: (1,0) - Zyy: (1,1)   

        Errors are given as standard deviations (sqrt(VAR))

    """

    def __init__(self, z_array = None, zerr_array = None, edi_object = None, Hscale = 1):

        if Hscale in ['b','B']:
            Hscale = 10e-3 * MTc.mu0
        try:
            float(Hscale)
        except:
            raise MTexceptions.MTpyError_edi_file('ERROR - H field scaling factor not undertood ')
    
    

        import mtpy.core.edi as MTedi 

        self.z = None
        self.zerr = None
              
        try:
            if len(z_array.shape) == 3 and z_array.shape[1:3] == (2,2):
                if z_array.dtype in ['complex', 'float']:
                    self.z = z_array
        except:
            pass

        try:
            if len(z_array.shape) == 2 and z_array.shape == (2,2):
                if z_array.dtype in ['complex', 'float']:
                    self.z = np.zeros((1,2,2),'complex')
                    self.z[0] = z_array            
        except:
            pass

        try:
            if len(zerr_array.shape) == 3 and zerr_array.shape[1:3] == (2,2):
                if zerr_array.dtype in ['float']:
                    self.zerr = zerr_array
        except:
            pass
        try:
            if len(zerr_array.shape) == 2 and zerr_array.shape == (2,2):
                if zerr_array.dtype in ['float']:
                    self.z = np.zeros((1,2,2))
                    self.zerr = zerr_array
        except:
            pass


            
        self.frequencies = None
        self.edi_object = None

        if isinstance(edi_object,MTedi.Edi):
            self.edi_object = edi_object
            self.frequencies = edi_object.frequencies
            self.z = edi_object.z
            if  edi_object.zerr  is not None:
                self.zerr = edi_object.zerr
        try:
            if len(self.z) != len(zerr):
                self.zerr = None
        except:
            pass


        #if (self.z is not None) and (self.zerr is None):
        #    self.zerr = np.zeros(self.z.shape)

        self.rotation_angle = 0.


    def set_edi_object(self, edi_object):

        if not isinstance(edi_object,MTedi.Edi):
            print 'Object is not a valid instance of the Edi class - Z object not updated'
            return


        self.edi_object = edi_object
        self.frequencies = edi_object.frequencies

        try:
            if edi_object.z is None :
                raise
            
            z_new = edi_object.z
            try:
                zerr_new = edi_object.zerr
            except:
                zerr_new = np.zeros(z_new.shape)

            if len(z_new) != len(zerr_new):
                raise
            self.z = z_new
            self.zerr = zerr_new

        except:
            print 'Edi object does not contain correct z information - z object not updated'


       
    def set_z(self, z_array):

        z_orig = self.z 

        if (self.z is not None) and (self.z.shape != z_array.shape):
            print 'Error - shape of "z" array does not match shape of Z array: %s ; %s'%(str(z_array.shape),str(self.z.shape))
            return

        self.z = z_array


    def set_zerr(self, zerr_array):

        if (self.zerr is not None) and (self.zerr.shape != zerr_array.shape):
            print 'Error - shape of "zerr" array does not match shape of Zerr array: %s ; %s'%(str(zerr_array.shape),str(self.zerr.shape))
            return

        self.zerr = zerr_array


    def real(self):
        if self.z is None:
            print 'z array is None - cannot calculate real'
            return

        return np.real(self.z)

        
    def set_real(self, real_array):
        

        if (self.z is not None) and (self.z.shape != real_array.shape):
            print 'shape of "real" array does not match shape of Z array: %s ; %s'%(str(real_array.shape),str(self.z.shape))
            return

        #assert real array:
        if np.linalg.norm(np.imag(real_array )) != 0 :
            print 'Error - array "real" is not real valued !'
            return

        i = np.complex(0,1)

        if self.z is not None:
            z_new = real_array + i * self.imag() 
        else:
            z_new = real_array 
                 

        self.z = z_new


    def imag(self):
        if self.z is None:
            print 'z array is None - cannot calculate imag'
            return


        return np.imag(self.z)

        
    def set_imag(self, imag_array):


        if (self.z is not None) and (self.z.shape != imag_array.shape):
            print 'Error - shape of "imag" array does not match shape of Z array: %s ; %s'%(str(imag_array.shape),str(self.z.shape))
            return
        
        #assert real array:
        if np.linalg.norm(np.imag(imag_array )) != 0 :
            print 'Error - array "imag" is not real valued !'
            return

        i = np.complex(0,1)
        
        if self.z is not None:
            z_new = self.real() + i * imag_array 
        else:
            z_new = i * imag_array 

        self.z = z_new


    def rho_phi(self):
        
        if self.z is None:
            print 'Z array is None - cannot calculate rho/phi'
            return
        rhoerr = None
        phierr = None
        if self.zerr is not None:
            rhoerr = np.zeros(self.zerr.shape)
            phierr = np.zeros(self.zerr.shape)

        rho = np.zeros(self.z.shape)
        phi = np.zeros(self.z.shape)


        for idx_f in range(len(self.z)): 
            omega = self.freq[idx_f] * 2 *math.pi

            for i in range(2):                        
                for j in range(2):
                    rho[idx_f,i,j] = np.abs(self.z[idx_f,i,j])**2 /omega/ MTc.mu0
                    phi[idx_f,i,j] = math.degrees(cmath.phase(self.z[idx_f,i,j]))%360
                
                    if self.zerr is not None:
                        r_err, phi_err = MTc.propagate_error_rect2polar( np.real(self.z[idx_f,i,j]), self.zerr[idx_f,i,j], np.imag(self.z[idx_f,i,j]), self.zerr[idx_f,i,j])
                        rhoerr[idx_f,i,j] = 2 * np.abs(self.z[idx_f,i,j])/omega/ MTc.mu0 * r_err
                        phierr[idx_f,i,j] = phi_err

        return rho, phi, rhoerr, phierr



    def set_rho_phi(self, rho_array, phi_array):

        if self.z is not None: 
            z_new = copy.copy(self.z) 

            if self.z.shape != rho_array.shape:
                print 'Error - shape of "rho" array does not match shape of Z array: %s ; %s'%(str(rho_array.shape),str(self.z.shape))
                return

            if self.z.shape != phi_array.shape:
                print 'Error - shape of "phi" array does not match shape of Z array: %s ; %s'%(str(phi_array.shape),str(self.z.shape))
                return
        else:
            z_new = p.zeros(rho_array.shape,'complex')
            if rho_array.shape != phi_array.shape:
                print 'Error - shape of "phi" array does not match shape of "rho" array: %s ; %s'%(str(phi_array.shape),str(rho_array.shape))
                return


        if (self.freq is None) or (len(self.freq) != len(rho_array)) :
            raise MTexceptions.MTpyError_EDI('ERROR -cannot set rho without proper frequency information - proper "freq" attribute must be defined ')

            
        #assert real array:
        if np.linalg.norm(np.imag(rho_array )) != 0 :
            print 'Error - array "rho" is not real valued !'
            return
        if np.linalg.norm(np.imag(phi_array )) != 0 :
            print 'Error - array "phi" is not real valued !'
            return

        for idx_f in range(len(z_new)):
            for i in range(2):
                for j in range(2):
                    abs_z = np.sqrt(2 * math.pi * self.freq[idx_f] * rho_array[idx_f,i,j])
                    z_new[idx_f,i,j] = cmath.rect( abs_z, math.radians(phi_array[idx_f,i,j] ))

        self.z = z_new


    def inverse(self):
        if self.z is None :
            print 'z array is "None" - I cannot invert that'
            return
        
        inverse = copy.copy(self.z)
        for idx_f in range(len(inverse)):
            try:
                inverse[idx_f,:,:] = np.array( (np.matrix(self.z[idx_f,:,:])).I )
            except:
                raise MTexceptions.MTpyError_Z('The %ith impedance tensor cannot be inverted'%(idx_f+1))

        return inverse


    def rotate(self, alpha):
        if self.z is None :
            print 'z array is "None" - I cannot rotate that'
            return

        #check for iterable list/set of angles - if so, it must have length 1 or same as len(tipper):
        if np.iterable(alpha) == 0:
            try:
                degreeangle = float(alpha%360)
            except:
                print '"Angle" must be a valid number (in degrees)'
                return

            #make an n long list of identical angles
            lo_angles = [degreeangle for i in self.z]
        else:
            if len(alpha) == 1:
                try:
                    degreeangle = float(alpha%360)
                except:
                    print '"Angle" must be a valid number (in degrees)'
                    return
                #make an n long list of identical angles
                lo_angles = [degreeangle for i in self.z]
            else:                    
                try:
                    lo_angles = [ float(i%360) for i in alpha]
                except:
                    print '"Angles" must be valid numbers (in degrees)'
                    return
            
        self.rotation_angle = lo_angles

        if len(lo_angles) != len(self.z):
            print 'Wrong number Number of "angles" - need %i '%(len(self.z))
            #self.rotation_angle = 0.
            return

        z_rot = copy.copy(self.z)
        zerr_rot = copy.copy(self.zerr)

        for idx_freq in range(len(self.z)):
                    
            angle = lo_angles[idx_freq]
            if np.isnan(angle):
                angle = 0.

            if self.zerr is not None:
                z_rot[idx_freq], zerr_rot[idx_freq] = MTc.rotatematrix_incl_errors(self.z[idx_freq,:,:], angle, self.zerr[idx_freq,:,:])
            else:
                z_rot[idx_freq], zerr_rot = MTc.rotatematrix_incl_errors(self.z[idx_freq,:,:], angle)


        self.z = z_rot
        self.zerr = zerr_rot
    


    def no_ss(self, reduce_rho_factor_x = 1., reduce_rho_factor_y = 1.):
        """
        Remove the static shift by providing the correction factors for x and y components.
        (Factor can be determined by using the "Analysis" module for the impedance tensor)

        Assume the original observed tensor Z is built by a static shift S and an unperturbated "correct" Z0 :
            Z = S * Z0

        returns:
            S, Z0   (over all frequencies)

        """
        
        #check for iterable list/set of reduce_rho_factor_x - if so, it must have length 1 or same as len(z):
        if np.iterable(reduce_rho_factor_x) == 0:
            try:
                x_factor = float(reduce_rho_factor_x)
            except:
                print '"reduce_rho_factor_x" must be a valid numbers'
                return

            lo_x_factors = [x_factor for i in self.z]
        else:
            if len(reduce_rho_factor_x) == 1:
                try:
                    x_factor = float(reduce_rho_factor_x)
                except:
                    print '"reduce_rho_factor_x" must be a valid numbers'
                    return
                lo_x_factors = [x_factor for i in self.z]
            else:                    
                try:
                    lo_x_factors = [x_factor for i in reduce_rho_factor_x]
                except:
                    print '"reduce_rho_factor_x" must be valid numbers'
                    return
            
        if len(lo_x_factors) != len(self.z):
            print 'Wrong number Number of "reduce_rho_factor_x" - need %i '%(len(self.z))
            return
  
        #check for iterable list/set of reduce_rho_factor_y - if so, it must have length 1 or same as len(z):
        if np.iterable(reduce_rho_factor_y) == 0:
            try:
                y_factor = float(reduce_rho_factor_y)
            except:
                print '"reduce_rho_factor_y" must be a valid numbers'
                return

            lo_y_factors = [y_factor for i in self.z]
        else:
            if len(reduce_rho_factor_y) == 1:
                try:
                    y_factor = float(reduce_rho_factor_y)
                except:
                    print '"reduce_rho_factor_y" must be a valid numbers'
                    return
                lo_y_factors = [y_factor for i in self.z]
            else:                    
                try:
                    lo_y_factors = [y_factor for i in reduce_rho_factor_y]
                except:
                    print '"reduce_rho_factor_y" must be valid numbers'
                    return
            
        if len(lo_y_factors) != len(self.z):
            print 'Wrong number Number of "reduce_rho_factor_y" - need %i '%(len(self.z))
            return
  

        z_corrected = copy.copy(self.z)
        static_shift = np.zeros((len(self.z),2,2))

        for idx_f in range(len(self.z)):
            z_corrected[idx_f,0,:] = self.z[idx_f,0,:]*np.sqrt(lo_x_factors[idx_f])
            z_corrected[idx_f,1,:] = self.z[idx_f,1,:]*np.sqrt(lo_y_factors[idx_f])
            static_shift[idx_f,0,0] = np.sqrt(lo_x_factors[idx_f])
            static_shift[idx_f,1,1] = np.sqrt(lo_y_factors[idx_f])

        return  static_shift, z_corrected,  z_corrected_err



    def no_distortion(self, distortion_tensor, distortion_err_tensor = None):
        """
            Remove distortion D form an observed impedance tensor Z to obtain the uperturbed "correct" Z0:
            Z = D * Z0

            Propagation of errors/uncertainties included
        """

        if distortion_err_tensor is None:
            distortion_err_tensor = np.zeros_like(distortion_tensor)
        #for all frequencies, calculate D.Inverse, then obtain Z0 = D.I * Z
        try:
            if not ( len(distortion_tensor.shape) in [2,3] ) and  (len(distortion_err_tensor.shape) in [2,3]):
                raise
            if len(distortion_tensor.shape) == 3 or len(distortion_err_tensor.shape) == 3:
                print 'Distortion is not time-dependent - take only first of given distortion tensors'
                try:
                    distortion_tensor = distortion_tensor[0]
                    distortion_err_tensor = distortion_err_tensor[0]
                except:
                    raise

            if not (distortion_tensor.shape == (2,2) ) and (distortion_err_tensor.shape == (2,2) ):
                raise

            distortion_tensor = np.matrix(np.real(distortion_tensor))

        except:
            raise MTexceptions.MTpyError_Z('The array provided is not a proper distortion tensor')

        try: 
            DI = distortion_tensor.I
        except:
            raise MTexceptions.MTpyError_Z('The provided distortion tensor is singular - I cannot invert that!')

        #propagation of errors (using 1-norm) - step 1 - inversion of D:
        DI_err = np.zeros_like(distortion_err_tensor)

        #todo :include error on  determinant!!
        D_det = np.linalg.det(distortion_tensor)

        dummy, DI_err = MTc.invertmatrix_incl_errors(distortion_tensor, distortion_err_tensor)

        #propagation of errors - step 2 - product of D.inverse and Z; D.I * Z, making it 4 summands for each component:
        z_corrected = np.zeros_like(self.z)
        z_corrected_err = np.zeros_like(self.zerr)
 
        for idx_f in range(len(self.z)):
            z_corrected[idx_f] = np.array(np.dot( DI, np.matrix(self.z[idx_f]) ))           
            for i in range(2):
                for j in range(2):
                    z_corrected_err[idx_f,i,j] = np.sum(np.abs( np.array( [ DI_err[i,0] * self.z[idx_f,0,j],\
                                                                            DI[i,0] * self.zerr[idx_f,0,j],\
                                                                            DI_err[i,1] * self.z[idx_f,1,j],\
                                                                            DI[i,1] * self.zerr[idx_f,1,j] ]))) 

        return distortion_tensor , z_corrected, z_corrected_err



    def no_ss_no_distortion(self, rho_x = 1., rho_y = 1.):

        pass


    def only1d(self):

        z1d = copy.copy(self.z)

        for i in range(len(z1d)):
            z1d[i,0,0] = 0
            z1d[i,1,1] = 0
            sign01 = z1d[i,0,1]/np.abs(z1d[i,0,1])
            sign10 = z1d[i,1,0]/np.abs(z1d[i,1,0])
            mean1d = 0.5* (z1d[i,1,0]+z1d[i,0,1])
            z1d[i,0,1] = sign01 * mean1d
            z1d[i,1,0] = sign10 * mean1d

        return z1d


    def only2d(self):

        z2d = copy.copy(self.z)

        for i in range(len(z2d)):
            z2d[i,0,0] = 0
            z2d[i,1,1] = 0
            
        return z2d



    def invariants(self):

        invariants_dict = {}

        z1 = (self.z[:,0,1] - self.z[:,1,0])/2.
        invariants_dict['z1'] = z1 

        det_z = np.array( [np.linalg.det(i) for i in self.z ])
        invariants_dict['det'] = det_z
        
        det_real = np.array( [np.linalg.det(i) for i in self.real() ])
        invariants_dict['det_real'] = det_real
        
        det_imag = np.array( [np.linalg.det(i) for i in self.imag() ])
        invariants_dict['det_imag'] = det_imag

        trace_z = np.array( [np.linalg.trace(i) for i in self.z ])
        invariants_dict['trace'] = trace_z
        
        skew_z = np.array( [np.abs(trace_z[i]/2.)/np.abs(z1[i]) for i in range(len(z1)) ])
        invariants_dict['skew'] = skew_z
        
        norm_z = np.array( [np.linalg.norm(i) for i in self.z ])
        invariants_dict['norm'] = norm_z
        
        lambda_plus = np.array( [ z1[i] + np.sqrt(z1[i] * z1[i] - det_z[i]) for i in range(len(z1)) ])
        invariants_dict['lambda_plus'] = lambda_plus
        
        lambda_minus = np.array( [ z1[i] - np.sqrt(z1[i] * z1[i] - det_z[i]) for i in range(len(z1)) ])
        invariants_dict['lambda_minus'] = lambda_minus
        
        sigma_plus = np.array( [ 0.5*norm_z[i]**2 + np.sqrt( 0.25*norm_z[i]**4 + np.abs(det_z[i])**2) for i in range(len(norm_z)) ])
        invariants_dict['sigma_plus'] = sigma_plus
        
        sigma_minus = np.array( [ 0.5*norm_z[i]**2 - np.sqrt( 0.25*norm_z[i]**4 + np.abs(det_z[i])**2) for i in range(len(norm_z)) ])
        invariants_dict['sigma_minus'] = sigma_minus

        return invariants_dict


#------------------------

class Tipper(object):
    """
        Tipper class - generates a Tipper-object.


        Errors are given as standard deviations (sqrt(VAR))

    """

    def __init__(self, tipper_array = None, tippererr_array = None, edi_object = None):
    

        self.tipper = None        
        self.tippererr = None
        try:
            if len(tipper_array.shape) == 3 and tipper_array.shape[1:3] == (1,2):
                if tipper_array.dtype in ['complex', 'float']:
                    self.tipper = tipper_array
        except:
            pass
        try:
            if len(tippererr_array.shape) == 3 and tippererr_array.shape[1:3] == (1,2):
                if tippererr_array.dtype in ['float']:
                    self.tippererr = tippererr_array
        except:
            pass

        self.frequencies = None
        self.edi_object = None

        if isinstance(edi_object,MTedi.Edi):
            self.edi_object = edi_object
            self.frequencies = edi_object.frequencies
            if edi_object.tipper is not None:
                self.tipper = edi_object.tipper
            if edi_object.tippererr is not None:
                self.tippererr = edi_object.tippererr

        try:
            if len(self.tipper) != len(tippererr):
                self.tippererr = None
        except:
            pass


        self.rotation_angle = 0.


    def set_edi_object(self, edi_object):

        if not isinstance(edi_object,MTedi.Edi):
            print 'Object is not a valid Edi instance - Tipper object not updated'
            return


        self.edi_object = edi_object
        self.frequencies = edi_object.frequencies
        
        try:
            if edi_object.tipper is None or edi_object.tippererr is None:
                raise
            
            tipper_new = edi_object.tipper
            tippererr_new = edi_object.tippererr

            if len(tipper_new) != len(tippererr_new):
                raise
            self.tipper = tipper_new
            self.tippererr = tippererr_new

        except:
            print 'Edi object does not contain correct tipper information - Tipper object not updated'

        
    def set_tipper(self, tipper_array):

        if (self.tipper is not None) and (self.tipper.shape != tipper_array.shape):
            print 'Error - shape of "tipper" array does not match shape of tipper-array: %s ; %s'%(str(tipper_array.shape),str(self.tipper.shape))
            return

        self.tipper = tipper_array


    def set_tippererr(self, tippererr_array):


        if (self.tippererr is not None) and (self.tippererr.shape != tippererr_array.shape):
            print 'Error - shape of "tippererr" array does not match shape of tippererr array: %s ; %s'%(str(tippererr_array.shape),str(self.tippererr.shape))
            return

        self.tippererr = tippererr_array


    def real(self):
        if self.tipper is None:
            print 'tipper array is None - cannot calculate real'
            return

        return np.real(self.tipper)

        
    def set_real(self, real_array):
        

        if (self.tipper is not None ) and (self.tipper.shape != real_array.shape):
            print 'shape of "real" array does not match shape of tipper array: %s ; %s'%(str(real_array.shape),str(self.tipper.shape))
            return

        #assert real array:
        if np.linalg.norm(np.imag(real_array )) != 0 :
            print 'Error - array "real" is not real valued !'
            return


        i = np.complex(0,1)

        if self.tipper is not None:
            tipper_new = real_array + i * self.imag() 
        else:
            tipper_new = real_array

        self.tipper = tipper_new


    def imag(self):
        if self.tipper is None:
            print 'tipper array is None - cannot calculate imag'
            return

        return np.imag(self.tipper)

        
    def set_imag(self, imag_array):


        if (self.tipper is not None) and (self.tipper.shape != imag_array.shape):
            print 'Error - shape of "imag" array does not match shape of tipper array: %s ; %s'%(str(imag_array.shape),str(self.tipper.shape))
            return
        
        #assert real array:
        if np.linalg.norm(np.imag(imag_array )) != 0 :
            print 'Error - array "imag" is not real valued !'
            return

        i = np.complex(0,1)
        if self.tipper is not None:
            tipper_new = self.real() + i * imag_array 
        else:
            tipper_new = i * imag_array 


        self.tipper = tipper_new


    def rho_phi(self):
        
        if self.tipper is None:
            print 'tipper array is None - cannot calculate rho/phi'
            return
        rhoerr = None
        phierr = None
        if self.tippererr is not None:
            rhoerr = np.zeros(self.tippererr.shape)
            phierr = np.zeros(self.tippererr.shape)

        rho = np.zeros(self.tipper.shape)
        phi = np.zeros(self.tipper.shape)


        for idx_f in range(len(tipper)):                         
            for j in range(2):
                rho[idx_f,0,j] = np.abs(self.tipper[idx_f,0,j])
                phi[idx_f,0,j] = math.degrees(cmath.phase(self.tipper[idx_f,0,j]))
                
                if self.tippererr is not None:
                    r_err, phi_err = MTc.propagate_error_rect2polar( np.real(self.tipper[idx_f,0,j]), self.tippererr[idx_f,0,j], np.imag(self.tipper[idx_f,0,j]), self.tippererr[idx_f,0,j])
                    rhoerr[idx_f,0,j] = r_err
                    phierr[idx_f,0,j] = phi_err

        return rho, phi, rhoerr, phierr


    def set_rho_phi(self, rho_array, phi_array):

        if self.tipper is not None: 
                
            tipper_new = copy.copy(self.tipper) 

            if self.tipper.shape != rho_array.shape:
                print 'Error - shape of "rho" array does not match shape of tipper array: %s ; %s'%(str(rho_array.shape),str(self.tipper.shape))
                return

            if self.tipper.shape != phi_array.shape:
                print 'Error - shape of "phi" array does not match shape of tipper array: %s ; %s'%(str(phi_array.shape),str(self.tipper.shape))
                return
        else:

            tipper_new = np.zeros(rho_array.shape,'complex')

            if rho_array.shape != phi_array.shape:
                print 'Error - shape of "phi" array does not match shape of "rho" array: %s ; %s'%(str(phi_array.shape),str(rho_array.shape))
                return
       
        #assert real array:
        if np.linalg.norm(np.imag(rho_array )) != 0 :
            print 'Error - array "rho" is not real valued !'
            return
        if np.linalg.norm(np.imag(phi_array )) != 0 :
            print 'Error - array "phi" is not real valued !'
            return

        for idx_f in range(len(rho_array)):
                for j in range(2):
                    tipper_new[idx_f,0,j] = cmath.rect( rho_array[idx_f,0,j], math.radians(phi_array[idx_f,0,j] ))

        self.tipper = tipper_new


    def rotate(self, alpha):

        if self.tipper is None :
            print 'tipper array is "None" - I cannot rotate that'
            return

        #check for iterable list/set of angles - if so, it must have length 1 or same as len(tipper):
        if np.iterable(alpha) == 0:
            try:
                degreeangle = float(alpha%360)
            except:
                print '"Angle" must be a valid number (in degrees)'
                return

            #make an n long list of identical angles
            lo_angles = [degreeangle for i in self.tipper]
        else:
            if len(lo_angles) == 1:
                try:
                    degreeangle = float(alpha%360)
                except:
                    print '"Angle" must be a valid number (in degrees)'
                    return
                #make an n long list of identical angles
                lo_angles = [degreeangle for i in self.z]
            else:                    
                try:
                    lo_angles = [ float(i%360) for i in alpha]
                except:
                    print '"Angles" must be valid numbers (in degrees)'
                    return
            
        self.rotation_angle = lo_angles

        if len(lo_angles) != len(self.tipper):
            print 'Wrong number Number of "angles" - need %i '%(len(self.tipper))
            self.rotation_angle = 0.
            return

        tipper_rot = copy.copy(self.tipper)
        tippererr_rot = copy.copy(self.tippererr)

        for idx_freq in range(len(tipper_rot)):
            angle = lo_angles[idx_freq]

            if self.tippererr is not None:
                tipper_rot[idx_freq], tippererr_rot[idx_freq] =  MTc.rotatevector_incl_errors(self.tipper[idx_freq,:,:], angle,self.tippererr[idx_freq,:,:] )
            else:
                tipper_rot[idx_freq], tippererr_rot = MTc.rotatevector_incl_errors(self.tipper[idx_freq,:,:], angle)


 
        self.tipper = tipper_rot
        self.tippererr = tippererr_rot


#------------------------


def rotate_z(z_array, alpha, zerr_array = None):

    z_object = _read_z_array(z_array, zerr_array)

    z_object.rotate(alpha)

    return z_object.z, z_object.zerr



def remove_distortion(z_array, distortion_tensor, distortion_err_tensor, zerr_array = None):
    
    z_object = _read_z_array(z_array, zerr_array)
    
    distortion, z_corrected, z_corrected_err = z_object.no_distortion(distortion_tensor, distortion_err_tensor)

    return distortion, z_corrected, z_corrected_err, z_object.z


def remove_ss(z_array, zerr_array = None, rho_x = 1., rho_y = 1.):


    z_object = _read_z_array(z_array, zerr_array)

    z_corrected, static_shift = z_object.no_ss()

    return static_shift, z_corrected, z_object.z


def remove_ss_and_distortion(z_array, zerr_array = None, rho_x = 1., rho_y = 1.):
    pass



def z2rhophi(z_array, zerr_array = None):
    
    z_object = _read_z_array(z_array,zerr_array )

    return z_object.rho_phi()


def rotate_tipper(tipper_array, alpha, tippererr_array = None):

    tipper_object = _read_tipper_array(tipper_array, tippererr_array)

    tipper_object.rotate(alpha)

    return tipper_object.tipper, tipper_object.tippererr


def tipper2rhophi(tipper_array, tippererr_array = None):
    
    tipper_object = _read_tipper_array(tipper_array, tippererr_array )

    return tipper_object.rho_phi()


def _read_z_array(z_array, zerr_array = None):

    try:
        z_object = Z( z_array=z_array, zerr_array=zerr_array )
    except:
        raise MTexceptions.MTpyError_Z('Cannot generate Z instance - check z-array dimensions/type: (N,2,2)/complex ; %s'%(str(z_array.shape)))

    return z_object


def _read_tipper_array(tipper_array, tippererr_array = None):

    try:
        tipper_object = tipper( tipper_array=tipper_array, tippererr_array=tippererr_array )
    except:
        raise MTexceptions.MTpyError_tipper('Cannot generate Tipper instance - check dimensions/type: (N,1,2)/complex ; %s'%(str(tipper_array.shape)))

    return tipper_object

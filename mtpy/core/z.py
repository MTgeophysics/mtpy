#!/usr/bin/env python

"""
mtpy/mtpy/core/z.py

Contains classes and functions for handling impedance tensors (Z). 
 
    Class:
    "Z" contains information about an impedance tensor Z. 

        Methods:
        - read_edi_object
        - set_z
        - set_zerr
        - rho
        - phi
        - set_res_phase
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

        - _edi_object
        - set_tipper
        - set_tippererr
        - set_r_phi
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

        All errors are given as standard deviations (sqrt(VAR))


    """

    def __init__(self, z_array = None, zerr_array = None, edi_object = None):
        """
            Initialise an instance of the Z class.

            Optional input:
            z_array : Numpy array containing Z values
            zerr_array : Numpy array containing Z-error values (NOT variance, but stddev!)
            edi_object : instance of the MTpy Edi class

            Initialise the attributes with None
        """    

        import mtpy.core.edi as MTedi 

        self.z = None
        self.zerr = None
              
        try:
            if len(z_array.shape) == 3 and z_array.shape[1:3] == (2,2):
                if z_array.dtype in ['complex', 'float', 'int']:
                    self.z = z_array
        except:
            pass

        try:
            if len(z_array.shape) == 2 and z_array.shape == (2,2):
                if z_array.dtype in ['complex', 'float','int']:
                    self.z = np.zeros((1,2,2),'complex')
                    self.z[0] = z_array            
        except:
            pass

        try:
            if len(zerr_array.shape) == 3 and zerr_array.shape[1:3] == (2,2):
                if zerr_array.dtype in ['float', 'int']:
                    self.zerr = zerr_array
        except:
            pass
        try:
            if len(zerr_array.shape) == 2 and zerr_array.shape == (2,2):
                if zerr_array.dtype in ['float', 'int']:
                    self.z = np.zeros((1,2,2))
                    self.zerr = zerr_array
        except:
            pass

            
        self.frequencies = None
        self.edi_object = None

        if isinstance(edi_object, MTedi.Edi):
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


    def read_edi_object(self, edi_object):
        """
            Read in an instance of the MTpy Edi class.

            Update attributes "z, zerr"

        """


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
        """
            Set the attribute 'z'.

            Input:
            Z array

            Test for shape, but no test for consistency!

        """         

        z_orig = self.z 

        if (self.z is not None) and (self.z.shape != z_array.shape):
            print 'Error - shape of "z" array does not match shape of Z array: %s ; %s'%(str(z_array.shape),str(self.z.shape))
            return

        self.z = z_array


    def set_zerr(self, zerr_array):
        """
            Set the attribute 'zerr'.

            Input:
            Zerror array

            Test for shape, but no test for consistency!

        """ 
        if (self.zerr is not None) and (self.zerr.shape != zerr_array.shape):
            print 'Error - shape of "zerr" array does not match shape of Zerr array: %s ; %s'%(str(zerr_array.shape),str(self.zerr.shape))
            return

        self.zerr = zerr_array


    def _get_real(self):
        """
            Return the real part of Z.


        """ 

        if self.z is None:
            print 'z array is None - cannot calculate real'
            return

        return np.real(self.z)

        
    def _set_real(self, real_array):
        """
            Set the real part of 'z'.

            Input:
            Z-shaped, real valued array

            Test for shape, but no test for consistency!

        """         

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

    real = property(_get_real, _set_real, doc='Real part of Z')

    def _get_imag(self):
        """
            Return the imaginary part of Z.

        """

        if self.z is None:
            print 'z array is None - cannot calculate imag'
            return


        return np.imag(self.z)

        
    def _set_imag(self, imag_array):
        """
            Set the imaginary part of 'z'.

            Input:
            Z-shaped, real valued array

            Test for shape, but no test for consistency!

        """         


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

    imag = property(_get_imag, _set_imag, doc='Imaginary part of Z ')

    def _get_res_phase(self):
        """
            Return values for resistivity (rho - in Ohm m) and phase (phase - in degrees).

            Output is a 4-tuple of arrays:
            (Res, Phase, ResError, PhaseError)
        """ 
        
        if self.z is None:
            print 'Z array is None - cannot calculate Res/Phase'
            return
        reserr = None
        phaseerr = None
        if self.zerr is not None:
            reserr = np.zeros(self.zerr.shape)
            phaseerr = np.zeros(self.zerr.shape)

        res= np.zeros(self.z.shape)
        phase = np.zeros(self.z.shape)


        for idx_f in range(len(self.z)): 

            for i in range(2):                        
                for j in range(2):
                    res[idx_f,i,j] = np.abs(self.z[idx_f,i,j])**2 /self.freq[idx_f] *0.2
                    phase[idx_f,i,j] = math.degrees(cmath.phase(self.z[idx_f,i,j]))%360
                
                    if self.zerr is not None:
                        r_err, phi_err = MTc.propagate_error_rect2polar( np.real(self.z[idx_f,i,j]), self.zerr[idx_f,i,j], np.imag(self.z[idx_f,i,j]), self.zerr[idx_f,i,j])
                        reserr[idx_f,i,j] = 0.4 * np.abs(self.z[idx_f,i,j])/self.freq[idx_f] * r_err
                        phaseerr[idx_f,i,j] = phi_err

        return res, phase, reserr, phaseerr


    res_phase= property(_get_res_phase, doc='Resistivity and Phase angle of Z')


    def set_res_phase(self, res_array, phase_array):
        """
            Set values for resistivity (res - in Ohm m) and phase (phase - in degrees).

            Updates the attributes "z, zerr".

        """ 


        if self.z is not None: 
            z_new = copy.copy(self.z) 

            if self.z.shape != res_array.shape:
                print 'Error - shape of "res" array does not match shape of Z array: %s ; %s'%(str(res_array.shape),str(self.z.shape))
                return

            if self.z.shape != phase_array.shape:
                print 'Error - shape of "phase" array does not match shape of Z array: %s ; %s'%(str(phase_array.shape),str(self.z.shape))
                return
        else:
            z_new = p.zeros(res_array.shape,'complex')
            if res_array.shape != phase_array.shape:
                print 'Error - shape of "phase" array does not match shape of "res" array: %s ; %s'%(str(phase_array.shape),str(res_array.shape))
                return


        if (self.freq is None) or (len(self.freq) != len(res_array)) :
            raise MTexceptions.MTpyError_EDI('ERROR - cannot set res without proper frequency information - proper "freq" attribute must be defined ')

            
        #assert real array:
        if np.linalg.norm(np.imag(res_array )) != 0 :
            print 'Error - array "res" is not real valued !'
            return
        if np.linalg.norm(np.imag(phase_array )) != 0 :
            print 'Error - array "phase" is not real valued !'
            return

        for idx_f in range(len(z_new)):
            for i in range(2):
                for j in range(2):
                    abs_z = np.sqrt(5 * self.freq[idx_f] * res_array[idx_f,i,j])
                    z_new[idx_f,i,j] = cmath.rect( abs_z, math.radians(phase_array[idx_f,i,j] ))

        self.z = z_new


    def _get_inverse(self):
        """
            Return the inverse of Z.

            (no error propagtaion included yet)

        """

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

    inverse = property(_get_inverse, doc='Inverse of Z')

    def rotate(self, alpha):
        """
            Rotate  Z array. Change the rotation angles in Zrot respectively.

            Rotation angle must be given in degrees. All angles are referenced to geographic North, positive in clockwise direction. (Mathematically negative!)

            In non-rotated state, X refs to North and Y to East direction.

            Updates the attributes "z, zerr, zrot".

        """



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
    


    def no_ss(self, reduce_res_factor_x = 1., reduce_res_factor_y = 1.):
        """
        Remove the static shift by providing the respective correction factors for the resistivity in the x and y components.
        (Factors can be determined by using the "Analysis" module for the impedance tensor)

        Assume the original observed tensor Z is built by a static shift S and an unperturbated "correct" Z0 :
            Z = S * Z0

        returns:
            S, Z0   (over all frequencies)

        Note:
        The factors are on the resistivity scale, so the entries of the matrix "S" are given by their square-roots! 

        """
        
        #check for iterable list/set of reduce_res_factor_x - if so, it must have length 1 or same as len(z):
        if np.iterable(reduce_res_factor_x) == 0:
            try:
                x_factor = float(reduce_res_factor_x)
            except:
                print '"reduce_res_factor_x" must be a valid numbers'
                return

            lo_x_factors = [x_factor for i in self.z]
        else:
            if len(reduce_res_factor_x) == 1:
                try:
                    x_factor = float(reduce_res_factor_x)
                except:
                    print '"reduce_res_factor_x" must be a valid numbers'
                    return
                lo_x_factors = [x_factor for i in self.z]
            else:                    
                try:
                    lo_x_factors = [x_factor for i in reduce_res_factor_x]
                except:
                    print '"reduce_res_factor_x" must be valid numbers'
                    return
            
        if len(lo_x_factors) != len(self.z):
            print 'Wrong number Number of "reduce_res_factor_x" - need %i '%(len(self.z))
            return
  
        #check for iterable list/set of reduce_res_factor_y - if so, it must have length 1 or same as len(z):
        if np.iterable(reduce_res_factor_y) == 0:
            try:
                y_factor = float(reduce_res_factor_y)
            except:
                print '"reduce_res_factor_y" must be a valid numbers'
                return

            lo_y_factors = [y_factor for i in self.z]
        else:
            if len(reduce_res_factor_y) == 1:
                try:
                    y_factor = float(reduce_res_factor_y)
                except:
                    print '"reduce_res_factor_y" must be a valid numbers'
                    return
                lo_y_factors = [y_factor for i in self.z]
            else:                    
                try:
                    lo_y_factors = [y_factor for i in reduce_res_factor_y]
                except:
                    print '"reduce_res_factor_y" must be valid numbers'
                    return
            
        if len(lo_y_factors) != len(self.z):
            print 'Wrong number Number of "reduce_res_factor_y" - need %i '%(len(self.z))
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
        """
            Not implemented yet!!
        """

        pass


    def _get_only1d(self):
        """
            Return Z in 1D form.

            If Z is not 1D per se, the diagonal elements are set to zero, the off-diagonal elements keep their signs, but their absolute is set to the mean of the original Z off-diagonal absolutes.
        """

        z1d = copy.copy(self.z)

        for i in range(len(z1d)):
            z1d[i,0,0] = 0
            z1d[i,1,1] = 0
            sign01 = np.sign(z1d[i,0,1])
            sign10 = np.sign(z1d[i,1,0])
            mean1d = 0.5* (z1d[i,1,0]+z1d[i,0,1])
            z1d[i,0,1] = sign01 * mean1d
            z1d[i,1,0] = sign10 * mean1d

        return z1d

    only1d = property(_get_only1d, doc=""" Return Z in 1D form. If Z is not 1D per se, the diagonal elements are set to zero, the off-diagonal elements keep their signs, but their absolute is set to the mean of the original Z off-diagonal absolutes.""")

    def _get_only2d(self):
        """
            Return Z in 2D form.

            If Z is not 2D per se, the diagonal elements are set to zero.
        """

        z2d = copy.copy(self.z)

        for i in range(len(z2d)):
            z2d[i,0,0] = 0
            z2d[i,1,1] = 0
            
        return z2d
    
    only2d = property(_get_only2d, doc="""Return Z in 2D form. If Z is not 2D per se, the diagonal elements are set to zero. """)


    def _get_trace(self):
        """
            Return the trace of Z (incl. uncertainties).

            Output:
            - Trace(Z) - Numpy array
            - Error of Trace(Z) - Numpy array

        """

        tr = np.array( [np.trace(i) for i in self.z])

        tr_err = None
        if self.zerr is not None:
            tr_err = np.zeros_like(tr)
            tr_err[:] = self.zerr[:,0,0] + self.zerr[:,1,1]


        return tr, tr_err

    trace = property(_get_trace, doc='Trace of Z, incl. error')

    def _get_skew(self):
        """
            Return the skew of Z (incl. uncertainties).

            Output:
            - Skew(Z) - Numpy array
            - Error of Skew(Z) - Numpy array

        """
        
        skew =  np.array( [ i[0,1] - i[1,0] for i in self.z ] )
        
        skewerr = None
        if self.zerr is not None:
            skewerr = np.zeros_like(skew)
            skewerr[:] = self.zerr[:,0,1] + self.zerr[:,1,0]

        return skew, skewerr
    skew = property(_get_skew, doc='Skew of Z, incl. error')

    def _get_det(self):
        """
            Return the determinant of Z (incl. uncertainties).

            Output:
            - det(Z) - Numpy array
            - Error of det(Z) - Numpy array

        """

        det_Z = np.array( [np.linalg.det(i) for i in self.z])
        
        det_Z_err = None
        if self.zerr is not None:
            det_Z_err = np.zeros_like(det_phi)
            det_Z_err[:] = np.abs(self.z[:,1,1] * self.zerr[:,0,0]) + np.abs(self.z[:,0,0] * self.zerr[:,1,1]) + np.abs(self.z[:,0,1] * self.zerr[:,1,0]) + np.abs(self.z[:,1,0] * self.zerr[:,0,1])

        return det_Z, det_Z_err
    det = property(_get_det, doc='Determinant of Z, incl. error')


    def _get_norm(self):
        """
            Return the 2-/Frobenius-norm of Z (incl. uncertainties).

            Output:
            - Norm(Z) - Numpy array
            - Error of Norm(Z) - Numpy array

        """

        znormerr = None
        norm_z = np.array( [np.linalg.norm(i) for i in self.z ])

        return znorm, znormerr
    norm = property(_get_norm, doc='Norm of Z, incl. error')


    def _get_invariants(self):
        """
            Return a dictionary of Z-invariants.

            Contains:
            z1, det, det_real, det_imag, trace, skew, norm, lambda_plus/minus, sigma_plus/minus
        """


        invariants_dict = {}

        z1 = (self.z[:,0,1] - self.z[:,1,0])/2.
        invariants_dict['z1'] = z1 

        invariants_dict['det'] = self.det()[0]
        
        det_real = np.array( [np.linalg.det(i) for i in self.real() ])
        invariants_dict['det_real'] = det_real
        
        det_imag = np.array( [np.linalg.det(i) for i in self.imag() ])
        invariants_dict['det_imag'] = det_imag

        invariants_dict['trace'] = trace()[0]
        
        invariants_dict['skew'] = skew()[0]
        
        invariants_dict['norm'] = norm()[0]
        
        lambda_plus = np.array( [ z1[i] + np.sqrt(z1[i] * z1[i] - self.det()[i]) for i in range(len(z1)) ])
        invariants_dict['lambda_plus'] = lambda_plus
        
        lambda_minus = np.array( [ z1[i] - np.sqrt(z1[i] * z1[i] - self.det()[i]) for i in range(len(z1)) ])
        invariants_dict['lambda_minus'] = lambda_minus
        
        sigma_plus = np.array( [ 0.5*norm_z[i]**2 + np.sqrt( 0.25*norm_z[i]**4 + np.abs(self.det()[i])**2) for i in range(len(norm_z)) ])
        invariants_dict['sigma_plus'] = sigma_plus
        
        sigma_minus = np.array( [ 0.5*norm_z[i]**2 - np.sqrt( 0.25*norm_z[i]**4 + np.abs(self.det()[i])**2) for i in range(len(norm_z)) ])
        invariants_dict['sigma_minus'] = sigma_minus

        return invariants_dict
    invariants = property(_get_invariants, doc='Invariants of Z: z1, det, det_real, det_imag, trace, skew, norm, lambda_plus/minus, sigma_plus/minus')

#------------------------

class Tipper(object):
    """
        Tipper class - generates a Tipper-object.


        Errors are given as standard deviations (sqrt(VAR))

    """

    def __init__(self, tipper_array = None, tippererr_array = None, edi_object = None):
        """
            Initialise an instance of the Tipper class.

            Optional input:
            tipper_array : Numpy array containing Tipper values
            tippererr_array : Numpy array containing Tipper-error values (NOT variance, but stddev!)
            edi_object : instance of the MTpy Edi class

            Initialise the attributes with None
        """    

        self.tipper = None        
        self.tippererr = None
        try:
            if len(tipper_array.shape) == 3 and tipper_array.shape[1:3] == (1,2):
                if tipper_array.dtype in ['complex', 'float','int']:
                    self.tipper = tipper_array
        except:
            pass
        try:
            if len(tippererr_array.shape) == 3 and tippererr_array.shape[1:3] == (1,2):
                if tippererr_array.dtype in ['float','int']:
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


    def read_edi_object(self, edi_object):
        """
            Read in an instance of the MTpy Edi class.

            Update attributes "tipper, tippererr"

        """

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
        """
            Set the attribute 'tipper'.

            Input:
            Tipper array

            Test for shape, but no test for consistency!

        """         

        if (self.tipper is not None) and (self.tipper.shape != tipper_array.shape):
            print 'Error - shape of "tipper" array does not match shape of tipper-array: %s ; %s'%(str(tipper_array.shape),str(self.tipper.shape))
            return

        self.tipper = tipper_array


    def set_tippererr(self, tippererr_array):
        """
            Set the attribute 'tippererr'.

            Input:
            TipperError array

            Test for shape, but no test for consistency!

        """         


        if (self.tippererr is not None) and (self.tippererr.shape != tippererr_array.shape):
            print 'Error - shape of "tippererr" array does not match shape of tippererr array: %s ; %s'%(str(tippererr_array.shape),str(self.tippererr.shape))
            return

        self.tippererr = tippererr_array


    def _get_real(self):
        """
            Return the real part of the Tipper.

        """ 
        if self.tipper is None:
            print 'tipper array is None - cannot calculate real'
            return

        return np.real(self.tipper)

        
    def _set_real(self, real_array):
        """
            Set the real part of 'tipper'.

            Input:
            Tipper-shaped, real valued array

            Test for shape, but no test for consistency!

        """         
        

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

    real = property(_get_real, _set_real, doc='Real part of the Tipper')

    def _get_imag(self):
        """
            Return the imaginary part of the Tipper.

        """ 

        if self.tipper is None:
            print 'tipper array is None - cannot calculate imag'
            return

        return np.imag(self.tipper)

        
    def _set_imag(self, imag_array):
        """
            Set the imaginary part of 'tipper'.

            Input:
            Tipper-shaped, real valued array

            Test for shape, but no test for consistency!

        """         


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

    imag = property(_get_imag, _set_imag, doc='Imaginary part of the Tipper')

    def _get_rho_phi(self):
        """
            Return values for amplitude (rho) and argument (phi - in degrees).

            Output is a 4-tuple of arrays:
            (Rho, Phi, RhoError, PhiError)
        """ 
 

        
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

    rho_phi = property(_get_rho_phi, doc='Amplitude and Phase angle of the Tipper')

    def set_rho_phi(self, r_array, phi_array):
        """
            Set values for amplitude(r) and argument (phi - in degrees).

            Updates the attributes "tipper, tippererr".

        """ 

        if self.tipper is not None: 
                
            tipper_new = copy.copy(self.tipper) 

            if self.tipper.shape != r_array.shape:
                print 'Error - shape of "r" array does not match shape of tipper array: %s ; %s'%(str(r_array.shape),str(self.tipper.shape))
                return

            if self.tipper.shape != phi_array.shape:
                print 'Error - shape of "phi" array does not match shape of tipper array: %s ; %s'%(str(phi_array.shape),str(self.tipper.shape))
                return
        else:

            tipper_new = np.zeros(r_array.shape,'complex')

            if r_array.shape != phi_array.shape:
                print 'Error - shape of "phi" array does not match shape of "r" array: %s ; %s'%(str(phi_array.shape),str(r_array.shape))
                return
       
        #assert real array:
        if np.linalg.norm(np.imag(r_array )) != 0 :
            print 'Error - array "r" is not real valued !'
            return
        if np.linalg.norm(np.imag(phi_array )) != 0 :
            print 'Error - array "phi" is not real valued !'
            return

        for idx_f in range(len(r_array)):
                for j in range(2):
                    tipper_new[idx_f,0,j] = cmath.rect( r_array[idx_f,0,j], math.radians(phi_array[idx_f,0,j] ))

        self.tipper = tipper_new


    def rotate(self, alpha):
        """
            Rotate  Tipper array. Change the rotation angles in Zrot respectively.

            Rotation angle must be given in degrees. All angles are referenced to geographic North, positive in clockwise direction. (Mathematically negative!)

            In non-rotated state, X refs to North and Y to East direction.

            Updates the attributes "tipper, tippererr, zrot".

        """

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
    """
        Rotate a Z array

        Input:
        - Z array : (1,2,2) or (2,2) shaped Numpy array
        - rotation angle (in degrees) for clockwise rotation

        Optional:
        - Zerror : (1,2,2) or (2,2) shaped Numpy array

        Output:
        - rotated Z array
        - rotated Zerror array (or None, if no error given)

    """

    z_object = _read_z_array(z_array, zerr_array)

    z_object.rotate(alpha)

    return z_object.z, z_object.zerr



def remove_distortion(z_array, distortion_tensor, distortion_err_tensor = None, zerr_array = None):
    """
        Remove the distortion from a given Z array.

        Inputs:
        - Z array : (1,2,2) or (2,2) shaped Numpy array
        - distortion_tensor : (1,2,2) or (2,2) shaped Numpy array

        Optional:
        - Zerror array : (1,2,2) or (2,2) shaped Numpy array
        - distortion_error_tensor : (1,2,2) or (2,2) shaped Numpy array

        Output:
        - corrected Z array
        - Error of corrected Z array (or None)
        - original Z array

    """
    
    z_object = _read_z_array(z_array, zerr_array)
    
    distortion, z_corrected, z_corrected_err = z_object.no_distortion(distortion_tensor, distortion_err_tensor)

    return  z_corrected, z_corrected_err, z_array


def remove_ss(z_array, zerr_array = None, res_x = 1., res_y = 1.):
    """
        Remove the static shift from a given Z array.

        Inputs:
        - Z array : (1,2,2) or (2,2) shaped Numpy array

        Optional:
        - Zerror array : (1,2,2) or (2,2) shaped Numpy array
        - res_x : factor, by which the X component of the Resistivity is off
        - res_y : factor, by which the Y component of the Resistivity is off
    
        Output:
        - corrected Z array
        - static shift tensor
        - original Z array

    """
 

    z_object = _read_z_array(z_array, zerr_array)

    z_corrected, static_shift = z_object.no_ss()

    return z_corrected, static_shift, z_array


def remove_ss_and_distortion(z_array, zerr_array = None, res_x = 1., res_y = 1.):
    """
        Not implemented yet !!

    """
    pass



def z2resphi(z_array, zerr_array = None):
    """
        Return the resistivity/phase information for Z 
        (Z must be given in units "km/s" = "mu m / nT"!!).

        Input:
        - Z array

        Optional:
        - Zerror array

        Output:
        - Resistivity array 
        - Phase array (in degrees)
        - Resistivity uncertainties array
        - Phase uncertainties array
    """
    
    z_object = _read_z_array(z_array,zerr_array)

    return z_object.res_phase()


def rotate_tipper(tipper_array, alpha, tippererr_array = None):
    """
        Rotate a Tipper array

        Input:
        Tipper array : (1,2,2) or (2,2) shaped Numpy array
        rotation angle (in degrees) for clockwise rotation

        Optional:
        TipperError : (1,2,2) or (2,2) shaped Numpy array

        Output:
        - rotated Tipper array
        - rotated TipperError array (or None, if no error given)

    """

    tipper_object = _read_tipper_array(tipper_array, tippererr_array)

    tipper_object.rotate(alpha)

    return tipper_object.tipper, tipper_object.tippererr


def tipper2rhophi(tipper_array, tippererr_array = None):
    """
        Return values for amplitude (rho) and argument (phi - in degrees).

        Output is a 4-tuple of arrays:
        (Rho, Phi, RhoError, PhiError)
    """ 
     
    tipper_object = _read_tipper_array(tipper_array, tippererr_array )

    return tipper_object.rho_phi()


def _read_z_array(z_array, zerr_array = None):
    """
        Read a Z array and return an instance of the Z class.


        Input:
        - Z array

        Optional:
        - Zerror array
    """


    try:


        z_object = Z( z_array=z_array, zerr_array=zerr_array)
    except:
        raise MTexceptions.MTpyError_Z('Cannot generate Z instance - check z-array dimensions/type: (N,2,2)/complex ; %s'%(str(z_array.shape)))

    return z_object


def _read_tipper_array(tipper_array, tippererr_array = None):
    """
        Read a Tipper array and return an instance of the Tipper class.


        Input:
        - Tipper array

        Optional:
        - TipperError array

    """

    try:
        tipper_object = tipper( tipper_array=tipper_array, tippererr_array=tippererr_array )
    except:
        raise MTexceptions.MTpyError_tipper('Cannot generate Tipper instance - check dimensions/type: (N,1,2)/complex ; %s'%(str(tipper_array.shape)))

    return tipper_object


def correct4sensor_orientation(Z_prime, Bx=0, By=90, Ex=0, Ey=90, Z_prime_error = None):
    """
        Correct a Z-array for wrong orientation of the sensors.

        Assume, E' is measured by sensors orientated with the angles
            E'x: a
            E'y: b
        
        Assume, B' is measured by sensors orientated with the angles
            B'x: c
            B'y: d
         
        With those data, one obtained the impedance tensor Z':
            E' = Z' * B'

        Now we define change-of-basis matrices T,U so that
            E = T * E'
            B = U * B'

        =>   T contains the expression of the E'-basis in terms of E (the standard basis) 
        and  U contains the expression of the B'-basis in terms of B (the standard basis)
        The respective expressions for E'x-basis vector and E'y-basis vector are the columns of T.
        The respective expressions for B'x-basis vector and B'y-basis vector are the columns of U.

        We obtain the impedance tensor in default coordinates as: 

        E' = Z' * B' => T^(-1) * E = Z' * U^(-1) * B => E = T * Z' * U^(-1) * B => Z = T * Z' * U^(-1)


        Input:
        - Z_prime - 2x2 complex valued Numpy array (impedance tensor)
        - orientation of Bx-sensor (degrees) - default = 0
        - orientation of By-sensor (degrees) - default = 90
        - orientation of Ex-sensor (degrees) - default = 0
        - orientation of Ey-sensor (degrees) - default = 90

        Optional:
        - Z_prime_error - 2x2 real valued Numpy array (impedance tensor standard deviation)


        Output:
        - Z - 2x2 complex valued Numpy array (impedance tensor) in default orientation
        - Zerr - 2x2 real valued Numpy array (impedance tensor standard deviation) in default orientation

    """
    try:
        if len(Z_prime.shape) != 2:
            raise
        if Z_prime.shape != (2,2):
            raise
        
        if Z_prime.dtype not in ['complex', 'float', 'int']:
            raise

        Z_prime = np.matrix(Z_prime)

    except:
        raise MTexceptions.MTpyError_inputarguments('ERROR - Z array not valid! Must be 2x2 complex array')

    if Z_prime_error is not None:
        try:
            if len(Z_prime_error.shape) != 2:
                raise
            if Z_prime_error.shape != (2,2):
                raise
        
            if Z_prime_error.dtype not in [ 'float', 'int']:
                raise

        except:
            raise MTexceptions.MTpyError_inputarguments('ERROR - Z-error array not valid! Must be 2x2 real array')


    T = np.matrix(np.zeros((2,2)))
    U = np.matrix(np.zeros((2,2)))

    dummy1 = cmath.rect(1,math.radians(Ex))


    T[0,0] = np.real(dummy1)
    T[1,0] = np.imag(dummy1)
    dummy2 = cmath.rect(1,math.radians(Ey))
    T[0,1] = np.real(dummy2)
    T[1,1] = np.imag(dummy2)

    dummy3 = cmath.rect(1,math.radians(Bx))
    U[0,0] = np.real(dummy3)
    U[1,0] = np.imag(dummy3)
    dummy4 = cmath.rect(1,math.radians(By))
    U[0,1] = np.real(dummy4)
    U[1,1] = np.imag(dummy4)

    try:
        Z = np.array(np.dot(T,np.dot(Z_prime, U.I)))
    except:
        raise MTexceptions.MTpyError_inputarguments("ERROR - Given angles do not define basis for 2 dimensions - cannot convert Z'")

    Zerr = copy.copy(Z_prime_error)

    #TODO: calculate error propagation

    return Z, Zerr

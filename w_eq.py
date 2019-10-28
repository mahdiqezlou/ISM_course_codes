### a file to plot log W vs log N and find Doppler parameter from observation of equivalent width
import numpy as np
from scipy.optimize import brentq

class W_eq(object):
    
    def __init__(self, f, wave_length, gamma):

        self.f = f
        self.wave_length= wave_length
        self.gamma = gamma
        self.c = 3*10**5

    def get_tau_0(self, N, b):
    
        # draine Eq 9.10 
        # N in cm^-2
        # wave_length in Angstrom
        # b in kms^-1
       return 0.7580*(N/(1.0*10**13))*(self.f/0.4164)*(self.wave_length/1215.7)*(10/b)

    def get_W_lambda(self, N, b, part=None):
    
        # draine Eq 9.27
        # N in cm^-2
        # wave_length in Angstrom
        # b in kms^-1
        # gamma in s^-1

        tau_0 = self.get_tau_0(N=N, b=b)
        ind_1 = np.where(tau_0 < 1.25393)
        ind_2 = np.where(tau_0 >= 1.25393)

        w_eq_lambda = np.zeros_like(tau_0)
        
        w_eq_lambda[ind_1]= np.sqrt(np.pi)*(b/self.c)*(tau_0[ind_1]/(1+(tau_0[ind_1]/(2*np.sqrt(2)))))*self.wave_length

        
        w_eq_lambda[ind_2] = np.sqrt((4*b*b/(self.c*self.c))*np.log(tau_0[ind_2]/(np.log(2)))+ (b/self.c)*(self.gamma*self.wave_length/(self.c*10*13))*((tau_0[ind_2]-1.25393)/np.sqrt(np.pi)))*self.wave_length
        
        #w_eq_lambda= np.sqrt(np.pi)*(b/self.c)*(tau_0/(1+(tau_0/(2*np.sqrt(2)))))*self.wave_length
        #w_eq_lambda = np.sqrt((4*b*b/(self.c*self.c))*np.log(tau_0/(np.log(2)))+ (b/self.c)*(self.gamma*self.wave_length/(self.c*10*13))*((tau_0-1.25393)/np.sqrt(np.pi)))*self.wave_length



        if part == 'linear':
            print('linear part of C.O.G returned, first argiment is N, second is W_lambda')
            return (N[ind_1], w_eq_lambda[ind_1])
        else :
            return w_eq_lambda

    def _W_for_solving_b(self, N, b):

        """
        The same argument as get_W_lambda, exept that N is a float number not an array. It is just used
        by self.find_b().
        """
        tau_0 = self.get_tau_0(N=N, b=b)
        
        if tau_0 < 1.25393:

            w_eq_lambda= np.sqrt(np.pi)*(b/self.c)*(tau_0/(1+(tau_0/(2*np.sqrt(2)))))*self.wave_length

        else:
            
            w_eq_lambda = np.sqrt((4*b*b/(self.c*self.c))*np.log(tau_0/(np.log(2)))+ (b/self.c)*(self.gamma*self.wave_length/(self.c*10*13))*((tau_0-1.25393)/np.sqrt(np.pi)))*self.wave_length

        return w_eq_lambda

    
    def find_b(self, N, W):
        
        """
        N: colden found from a weak absorption of the same elemnt
        W: specific eqivalent width observed for the sronger absorption od the same element, Notice: specific eq width = dimensionless eq width * wave_length
        
        Returns : corrsponding Doppler parameter(b)
        """

        return brentq(lambda b: self._W_for_solving_b(b=b, N=N)-W, 1, 10)
    





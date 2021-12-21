import numpy as np
import sys as sys
import time
import scipy

#---------------------------------------------------------------------------------------------

def scale6npass10(eps10):
    """ 
    For a given eps for 10th order filter, this returns the eps value need for a 6-pass
    6th order filter to have the same R(eps10) = 0.5 filter response
    """
    return .125 * 4**(np.log10(eps10))

#---------------------------------------------------------------------------------------------

class TimerError(Exception):
    """A custom exception used to report errors in use of Timer class"""

class Timer:
    def __init__(self):
        self._start_time = None

    def start(self):
        """Start a new timer"""
        if self._start_time is not None:
            raise TimerError(f"Timer is running. Use .stop() to stop it")

        self._start_time = time.perf_counter()

    def stop(self):
        """Stop the timer, and report the elapsed time"""
        if self._start_time is None:
            raise TimerError(f"Timer is not running. Use .start() to start it")

        elapsed_time = time.perf_counter() - self._start_time
        self._start_time = None
        print(f"Elapsed time: {elapsed_time:0.4f} seconds")

#---------------------------------------------------------------------------------------------
            
def RaymondFilter6(xy2d, eps, npass = 1, **kwargs):
    """
    Pure python driver for Raymond's filters for 1 and 2D arrays
    
    Adapted from Raymond's original code and from example code on Program Creek.
       
    https://www.programcreek.com/python/example/97027/scipy.linalg.solve_banded [example 3]
    
    See RAYMOND, 1988, MWR, 116, 2132-2141 for details on parameter EPS
    
    Lou Wicker, Oct 2021 
    """
    
    try:
        from scipy.sparse import spdiags
        from scipy.sparse.linalg import spsolve
        from scipy import linalg        
    except ImportError: 
        raise ImportError("RaymondFilter6:  Filter1D requires scipy.sparse libraries.")

    #---------------------------------------------------------------------------------------------
    def Filter6_Init(N, EPS, bc_reflect=False, **kwargs):
        """
        Compute diagonal matrix values for 1D column Raymond's 1988 6th order implicit tangent filter.
        
        For a given 1D grid size, these are only needed once, so for 2D or 3D computations, this should
        increase the speed of the calculation.

        Adapted from Raymond's original code and from example code on Program Creek.

        https://www.programcreek.com/python/example/97027/scipy.linalg.solve_banded [example 3]

        See RAYMOND, 1988, MWR, 116, 2132-2141 for details on parameter EPS

        Lou Wicker, Oct 2021 
        
        Inputs:  N   - length of 1D array
                 EPS - filter scale value
                 
        Returns:  A - a special matrix of the form CRC that is used by the linear solver. 
        """
        
        if 'bc_reflect' in kwargs:
            bc_reflect = kwargs.get('bc_reflect')
            print("Filter6_Init:  bc_reflect = ", bc_reflect)

        NM1 = N-1
        NM2 = N-2
        NM3 = N-3
        NM4 = N-4
        NM5 = N-5
        
        EP  = 1.0 + EPS
        EM  = 1.0 - EPS

        # Initialize the matrix 
          
        U3 = np.zeros((N,), dtype=np.float64)
        U2 = np.zeros((N,), dtype=np.float64)
        U1 = np.zeros((N,), dtype=np.float64)

        D  = np.zeros((N,), dtype=np.float64)
        
        L1 = np.zeros((N,), dtype=np.float64)
        L2 = np.zeros((N,), dtype=np.float64)
        L3 = np.zeros((N,), dtype=np.float64)
        
        # compute the internal diagonals
        
        U3[:] =       EM
        U2[:] =   6.0*EP
        U1[:] =  15.0*EM
   
        D [:] =  20.0*EP
        
        L1[:] =  15.0*EM
        L2[:] =   6.0*EP
        L3[:] =       EM

        if bc_reflect: # use reflective condition at boundary
            
            U3[0:3] = (  2.0*EM,       EM,       EM)
            U2[0:3] = ( 12.0*EP,   6.0*EP,   6.0*EP)
            U1[0:3] = ( 30.0*EM,  16.0*EM,  15.0*EM)
              
            D [0:3] = ( 20.0*EP,  26.0*EP,  20.0*EP)  
            
            L1[0:3] = (  0.0,     15.0*EM,  16.0*EM)        
            L2[0:3] = (  0.0,      0.0,      6.0*EP)
            L3[0:3] = (  0.0,      0.0,        0.0 )


            U3[NM3:N] = (  0.0,      0.0,        0.0 )
            U2[NM3:N] = (  6.0*EP,   0.0,        0.0 )            
            U1[NM3:N] = ( 16.0*EM,  15.0*EM,     0.0 )
            
            D [NM3:N] = ( 20.0*EP,  26.0*EP,  20.0*EP)  
            
            L1[NM3:N] = ( 15.0*EM,  16.0*EM,  30.0*EM)
            L2[NM3:N] = (  6.0*EP,   6.0*EP,  12.0*EP)
            L3[NM3:N] = (      EM,       EM,   2.0*EM)
            
            # Convert diagonals to sparse matrix format (the tricky part as far as I am concerned!)
        
            diagonals = [D, L1, U1, L2, U2, L3, U3]
        
            return spdiags(diagonals, [0, -1, 1, -2, 2, -3, 3], N, N, format='csc')

        else:  # reduce order of approximation at boundary assuming Dirichet at i=0
        
            U3[0:3] = ( 0.0,     0.0,     0.0 )
            U2[0:3] = ( 0.0,     0.0,      EP )
            U1[0:3] = ( 0.0,      EM,  4.0*EM )

            D [0:3] = ( 1.0,  2.0*EP,  6.0*EP )

            L1[0:3] = ( 0.0,      EM,  4.0*EM )
            L2[0:3] = ( 0.0,     0.0,      EP )
            L3[0:3] = ( 0.0,     0.0,     0.0 )

            U3[NM3:N] = (    0.0,    0.0,     0.0 )
            U2[NM3:N] = (     EP,    0.0,     0.0 )
            U1[NM3:N] = ( 4.0*EM,     EM,     0.0 )

            D [NM3:N] = ( 6.0*EP, 2.0*EP,     1.0 )

            L1[NM3:N] = ( 4.0*EM,     EM,     0.0 )
            L2[NM3:N] = (     EP,    0.0,     0.0 )
            L3[NM3:N] = (    0.0,    0.0,     0.0 )
        
            # Convert diagonals to sparse matrix format (the tricky part as far as I am concerned!)
        
            diagonals = [D[1:-1], L1[1:-1], U1[1:-1], L2[1:-1], U2[1:-1], L3[1:-1], U3[1:-1]]
 
            # The use of Dirichet BCs means that the first and last values of solution == 0.  Rewrite system to be N-2.

            return spdiags(diagonals, [0, -1, 1, -2, 2, -3, 3], N-2, N-2, format='csc')
    
    #---------------------------------------------------------------------------------------------
    def Filter1D(XY, EPS, A, bc_reflect=False, **kwargs):
        """
        Compute solution for 1D column Raymond's 1988 6th order implicit tangent filter.

        Adapted from Raymond's original code and from example code on Program Creek.

        https://www.programcreek.com/python/example/97027/scipy.linalg.solve_banded [example 3]

        See RAYMOND, 1988, MWR, 116, 2132-2141 for details on parameter EPS

        Lou Wicker, Oct 2021 
        """
                  
        if 'bc_reflect' in kwargs:
            bc_reflect = kwargs.get('bc_reflect')
            print("Filter1D:  bc_reflect = ", bc_reflect)

        N  = len(XY)

        RHS = np.zeros((N,), dtype=np.float64)
        XF  = XY.copy()
        
        # Construct RHS (0-based indexing, not 1 like in fortran)
        
        NM1 = N-1
        NM2 = N-2
        NM3 = N-3
        NM4 = N-4
        NM5 = N-5
        NM6 = N-6
        
        # Compute inner RHS
        
        RHS[3:NM3] = EPS*((XY[0:NM3-3]+XY[6:NM3+3])
                    - 6.0*(XY[1:NM3-2]+XY[5:NM3+2])
                    +15.0*(XY[2:NM3-1]+XY[4:NM3+1])
                    -20.0* XY[3:NM3]               )

        if bc_reflect: 

            RHS[  0] =  EPS*(-20.0*XY[0] + 30.0*XY[1] - 12.0*XY[2] +  2.0*XY[3])
            RHS[  1] =  EPS*( 15.0*XY[0] - 26.0*XY[1] + 16.0*XY[2] -  6.0*XY[3] +     XY[4]) 
            RHS[  2] =  EPS*( -6.0*XY[0] + 16.0*XY[1] - 20.0*XY[2] + 15.0*XY[3] - 6.0*XY[4] + XY[5])
        

            RHS[NM3] = EPS*( -6.0*XY[NM1] + 16.0*XY[NM2] - 20.0*XY[NM3] + 15.0*XY[NM4] - 6.0*XY[NM5] + XY[NM6])
            RHS[NM2] = EPS*( 15.0*XY[NM1] - 26.0*XY[NM2] + 16.0*XY[NM3] -  6.0*XY[NM4] +     XY[NM5])   
            RHS[NM1] = EPS*(-20.0*XY[NM1] + 30.0*XY[NM2] - 12.0*XY[NM3] +  2.0*XY[NM4])
       
            XF = XF + spsolve(A, RHS)

            return XF

        else:
            
            RHS[  0] = 0.0
            RHS[  1] =  0.1*EPS*( XY[0] - 2.0*XY[1] +     XY[2] ) 
            RHS[  2] = -0.5*EPS*( XY[0] - 4.0*XY[1] + 6.0*XY[2] - 4.0*XY[3] + XY[4] )
        

            RHS[NM3] = -0.5*EPS*( XY[NM1] - 4.0*XY[NM2] + 6.0*XY[NM3] - 4.0*XY[NM4] + XY[NM5] )
            RHS[NM2] =  0.1*EPS*( XY[NM1] - 2.0*XY[NM2] +     XY[NM3] )   
            RHS[NM1] = 0.0 

            B = scipy.linalg.lu_factor(A.toarray())

            XF[1:-1] = XF[1:-1] + scipy.linalg.lu_solve(B, RHS[1:-1])
            
            print('X')


#            XF[1:-1] = XF[1:-1] + spsolve(A, RHS[1:-1])

            return XF

    #---------------------------------------------------------------------------------------------
    # Code to do 1D or 2D input

    if len(xy2d.shape) < 2:
        
        A = Filter6_Init(xy2d.shape[0], eps, **kwargs)
        
        print(f'\nCOND:  \n',np.linalg.cond(A.toarray()))
        
        x_copy = xy2d[:].copy()
        for n in np.arange(npass):
            x_copy = Filter1D(x_copy.copy(), eps, A, **kwargs)
        return x_copy
    
    elif len(xy2d.shape) > 2:
        print("RaymondFilter6:  3D filtering not implemented as of yet, exiting\n")
        sys.exit(-1)
        
    else:
    
        ny, nx = xy2d.shape
        
        print("RaymondFilter6 (Python version) called:  Shape of array:  NY: %d  NX:  %d, NPASS:  %d\n" % (ny, nx, npass))
        
        x1d = np.zeros((nx,))
        y1d = np.zeros((ny,))
        
        ytmp    = xy2d.copy()
        xy_copy = xy2d.copy()

        tic = time.perf_counter()

        for n in np.arange(npass):  # multiple pass capability
            
            tic = time.perf_counter()
            
            A = Filter6_Init(ny, eps, **kwargs)
            
            for i in np.arange(nx):
                y1d[:]    = xy_copy[:,i]
                ytmp[:,i] = Filter1D(y1d, eps, A, **kwargs)
        
            A = Filter6_Init(nx, eps, **kwargs)
            
            xtmp = ytmp.copy()
            
            for j in np.arange(ny):
                xy_copy[j,:] = Filter1D(xtmp[j,:], eps, A, **kwargs)
        
            toc = time.perf_counter()
            
            print(f"J-loop {tic - toc:0.4f} seconds")
        
        return xy_copy

#---------------------------------------------------------------------------------------------
            
def RaymondFilter6F(xy2d, eps, npass=1, **kwargs):
    """
    Combined Fortran and Python driver for Raymond's filters for 1 and 2D arrays
    
    Adapted from Raymond's original code 
    
    Wraps Raymond's 6th order fortran code in python using f2py
    
    See RAYMOND, 1988, MWR, 116, 2132-2141 for details on parameter EPS
    
    Lou Wicker, Oct 2021 
    """
    
    try:
        from filter.raymond_lowpass import raymond1d_lowpass, raymond2d_lowpass 
    except ImportError: 
        raise ImportError("RaymondFilter6F:  Requires compiled shared library called raymond_lowpass.")
            #---------------------------------------------------------------------------------------------
    # Code to do 1D or 2D input

    if len(xy2d.shape) ==1:
        
        toc = time.perf_counter()

        x_copy = xy2d.copy()
        for n in np.arange(npass):
            x_copy = raymond1d_lowpass(x_copy.copy(), eps)
            
        tic = time.perf_counter()
            
        print(f"NPass Loop for 1D array took {tic - toc:0.4f} seconds\n")

        return x_copy
    
    elif len(xy2d.shape) == 2:
            
        xy_copy = xy2d.copy()

        toc = time.perf_counter()

        for n in np.arange(npass):  # multiple pass capability
            
            xy_copy = raymond2d_lowpass(xy_copy.transpose(),eps)
                      
        tic = time.perf_counter()
            
        print(f"NPass Loop for 2D array took {tic - toc:0.4f} seconds\n")
        
        return xy_copy
        
    elif len(xy2d.shape) == 3:
        
        print("RaymondFilter6F:  Input array is 3D, 2D filtering implemented on outer two dimensions\n")
        
        xy_copy = xy2d.copy()
        
        for n in np.arange(npass):  # multiple pass capability
            
            for k in np.arange(xy2d.shape[0]):
                xy_copy[k] = raymond2d_lowpass(xy_copy[k].transpose().copy,eps)
                      
        tic = time.perf_counter()
            
        print(f"NPass Loop for 3D array took {tic - toc:0.4f} seconds\n")

        return xy_copy

#---------------------------------------------------------------------------------------------
    
def RaymondFilter10(xy2d, eps, **kwargs):
    """
    Driver for Raymond's filters for 1 and 2D arrays
    
    Adapted from Raymond's original code and from example code on Program Creek.
       
    https://www.programcreek.com/python/example/97027/scipy.linalg.solve_banded [example 3]
    
    See RAYMOND, 1988, MWR, 116, 2132-2141 for details on parameter EPS
    
    Lou Wicker, Oct 2021 
    """

    try:
        from scipy.sparse import spdiags
        from scipy.sparse.linalg import spsolve
    except ImportError: 
        raise ImportError("Raymond_10_Filter1D requires scipy.sparse libraries.")
        
        
    #---------------------------------------------------------------------------------------------
    def Filter10_Init(N, EPS, **kwargs):
        """
        Compute diagonal matrix values for 1D column Raymond's 1988 6th order implicit tangent filter.
        
        For a given 1D grid size, these are only needed once, so for 2D or 3D computations, this should
        increase the speed of the calculation.

        Adapted from Raymond's original code and from example code on Program Creek.

        https://www.programcreek.com/python/example/97027/scipy.linalg.solve_banded [example 3]

        See RAYMOND, 1988, MWR, 116, 2132-2141 for details on parameter EPS

        Lou Wicker, Oct 2021 
        
        Inputs:  N   - length of 1D array
                 EPS - filter scale value
                 
        Returns:  A - a special matrix of the form CRC that is used by the linear solver. 
        """

        # Construct RHS (0-based indexing, not 1 like in fortran)
        NM1 = N-1
        NM2 = N-2
        NM3 = N-3
        NM4 = N-4
        NM5 = N-5
        NM6 = N-6
        NM7 = N-7
        NM8 = N-7
        NM9 = N-8
        
        EP  = 1.0 + EPS
        EM  = 1.0 - EPS
        
        # Initialize the diagonal arrays 
          
        U5 = np.zeros((N,), dtype=np.float64)
        U4 = np.zeros((N,), dtype=np.float64)
        U3 = np.zeros((N,), dtype=np.float64)
        U2 = np.zeros((N,), dtype=np.float64)
        U1 = np.zeros((N,), dtype=np.float64)

        D  = np.zeros((N,), dtype=np.float64)
        
        L1 = np.zeros((N,), dtype=np.float64)
        L2 = np.zeros((N,), dtype=np.float64)
        L3 = np.zeros((N,), dtype=np.float64)
        L4 = np.zeros((N,), dtype=np.float64)
        L5 = np.zeros((N,), dtype=np.float64)

        # set internal diagonals

        U5[:] =          EM
        U4[:] =   10.0 * EP 
        U3[:] =   45.0 * EM 
        U2[:] =  120.0 * EP 
        U1[:] =  210.0 * EM
        
        D [:] =  252.0 * EP
        
        L1[:] =  210.0 * EM
        L2[:] =  120.0 * EP 
        L3[:] =   45.0 * EM 
        L4[:] =   10.0 * EP 
        L5[:] =          EM

        # Set boundary values - overwrite outer diagonal entries
                
        U5[0:5] = (0.0,       0.0,        0.0,        0.0,        0.0)
        U4[0:5] = (0.0,       0.0,        0.0,        0.0,         EP)
        U3[0:5] = (0.0,       0.0,        0.0,         EM,     8.0*EM)
        U2[0:5] = (0.0,       0.0,         EP,     6.0*EP,    28.0*EP)
        U1[0:5] = (0.0,        EM,     4.0*EM,    15.0*EM,    56.0*EM)
        
        D [0:5] = (1.0,    2.0*EP,     6.0*EP,    20.0*EP,    70.0*EP)

        L1[0:5] = (0.0,        EM,     4.0*EM,    15.0*EM,    56.0*EM)
        L2[0:5] = (0.0,       0.0,         EP,     6.0*EP,    28.0*EP)
        L3[0:5] = (0.0,       0.0,        0.0,         EM,     8.0*EM)
        L4[0:5] = (0.0,       0.0,        0.0,        0.0,         EP)
        L5[0:5] = (0.0,       0.0,        0.0,        0.0,        0.0)

        U5[NM5:N] = (     0.0,       0.0,       0.0,     0.0,     0.0)
        U4[NM5:N] = (      EP,       0.0,       0.0,     0.0,     0.0)
        U3[NM5:N] = (  8.0*EM,        EM,       0.0,     0.0,     0.0)
        U2[NM5:N] = ( 28.0*EP,    6.0*EP,        EP,     0.0,     0.0)
        U1[NM5:N] = ( 56.0*EM,   15.0*EM,    4.0*EM,      EM,     0.0)
        
        D [NM5:N] = ( 70.0*EP,   20.0*EP,    6.0*EP,  2.0*EP,     1.0)

        L1[NM5:N] = ( 56.0*EM,   15.0*EM,    4.0*EM,      EM,     0.0)
        L2[NM5:N] = ( 28.0*EP,    6.0*EP,        EP,     0.0,     0.0)
        L3[NM5:N] = (  8.0*EM,        EM,       0.0,     1.0,     0.0)
        L4[NM5:N] = (      EP,       0.0,       0.0,     0.0,     0.0)
        L5[NM5:N] = (     0.0,       0.0,       0.0,     0.0,     0.0)


#         U5[0:5] = (0.0,       0.0,        0.0,        0.0,        0.0)
#         U4[0:5] = (0.0,       0.0,        0.0,        0.0,        0.0)
#         U3[0:5] = (0.0,       0.0,        0.0,        0.0,        0.0)
#         U2[0:5] = (0.0,       0.0,        0.0,        0.0,        0.0)
#         U1[0:5] = (0.0,        EM,         EM,         EM,         EM)
        
#         D [0:5] = (1.0,    2.0*EP,     2.0*EP,     2.0*EP,     2.0*EP)

#         L1[0:5] = (0.0,        EM,         EM,         EM,         EM)
#         L2[0:5] = (0.0,       0.0,        0.0,        0.0,        0.0)
#         L3[0:5] = (0.0,       0.0,        0.0,        0.0,        0.0)
#         L4[0:5] = (0.0,       0.0,        0.0,        0.0,        0.0)
#         L5[0:5] = (0.0,       0.0,        0.0,        0.0,        0.0)

#         U5[NM5:N] = (     0.0,       0.0,       0.0,     0.0,     0.0)
#         U4[NM5:N] = (     0.0,       0.0,       0.0,     0.0,     0.0)
#         U3[NM5:N] = (     0.0,       0.0,       0.0,     0.0,     0.0)
#         U2[NM5:N] = (     0.0,       0.0,       0.0,     0.0,     0.0)
#         U1[NM5:N] = (      EM,        EM,        EM,      EM,     0.0)
        
#         D [NM5:N] = (  2.0*EP,    2.0*EP,    2.0*EP,  2.0*EP,     1.0)

#         L1[NM5:N] = (      EM,        EM,        EM,      EM,     0.0)
#         L2[NM5:N] = (     0.0,       0.0,       0.0,     0.0,     0.0)
#         L3[NM5:N] = (     0.0,       0.0,       0.0,     0.0,     0.0)
#         L4[NM5:N] = (     0.0,       0.0,       0.0,     0.0,     0.0)
#         L5[NM5:N] = (     0.0,       0.0,       0.0,     0.0,     0.0)

        # Convert diagonals to sparse matrix format (the tricky part as far as I am concerned!)

        diagonals = [D[1:-1],L1[1:-1],U1[1:-1],L2[1:-1],U2[1:-1],
                             L3[1:-1],U3[1:-1],L4[1:-1],U4[1:-1],L5[1:-1],U5[1:-1]]


        return spdiags(diagonals, [0, -1, 1, -2, 2, -3, 3, -4, 4, -5, 5], N-2, N-2, format='csc')
    
    #---------------------------------------------------------------------------------------------
    def Filter1D(XY, EPS, A, **kwargs):
        """
        Compute solution for 1D column Raymond's 1988 10th order implicit tangent filter.

        Adapted from Raymond's original code and from example code on Program Creek.

        https://www.programcreek.com/python/example/97027/scipy.linalg.solve_banded [example 3]

        See RAYMOND, 1988, MWR, 116, 2132-2141 for details on parameter EPS

        Equation

        (1-eps)*(p(i+2)+p(i-2)) + 10*(1+eps)(p(i+4)+p(i-4)) + 45*(1-eps)(p(i+3)+p(i-3)) 
        + 120*(1+eps)(p(i+2)+p(i-2)) + 210*(1-eps)(p(i+1)+p(i-1)) -252*(eps)(p(i))  

        =

        eps * ((u(i+2)+p(i-2)) - 010*u(i+4)+u(i-4)) + 45*(u(i+3)+u(i-3)) - 120*(u(i+2)+u(i-2))
         + 210*(u(i+2)+u(i-2)) - 252*u(i))

        [-252.  210.  210. -120. -120.   45.   45.  -10.  -10.    1.    1.]

        Lou Wicker, Oct 2021 
        """

        N = len(XY)
        RHS = np.zeros((N,),   dtype=np.float64)
        XF  = XY.copy()

        # Construct RHS (0-based indexing, not 1 like in fortran)
        NM1 = N-1
        NM2 = N-2
        NM3 = N-3
        NM4 = N-4
        NM5 = N-5
        NM6 = N-6
        NM7 = N-7
        NM8 = N-8
        NM9 = N-9

        # Compute RHS filter
        
        RHS[  0] =  0.0
        RHS[  1] =  0.01*EPS*( XY[0] - 2.0*XY[1] +      XY[2] ) 
        RHS[  2] = -0.05*EPS*( XY[0] - 4.0*XY[1] +  6.0*XY[2] -  4.0*XY[3] +      XY[4] )
        RHS[  3] =  0.1*EPS*( XY[0] - 6.0*XY[1] + 15.0*XY[2] - 20.0*XY[3] + 15.0*XY[4] -  6.0*XY[5] +      XY[6] )
        RHS[  4] = -0.5*EPS*( XY[0] - 8.0*XY[1] + 28.8*XY[2] - 56.0*XY[3] + 70.0*XY[4] - 56.0*XY[5] + 28.0*XY[6] - 8.0*XY[7] + XY[8] )

        RHS[NM5] = -0.5*EPS*( XY[NM1] - 8.0*XY[NM2] + 28.8*XY[NM3] - 56.0*XY[NM4] + 70.0*XY[NM5] - 56.0*XY[NM6] + 28.0*XY[NM7] - 8.0*XY[NM8] + XY[NM9] )        
        RHS[NM4] =  0.1*EPS*( XY[NM1] - 6.0*XY[NM2] + 15.0*XY[NM3] - 20.0*XY[NM4] + 15.0*XY[NM5] -  6.0*XY[NM6] +      XY[NM7] )
        RHS[NM3] = -0.05*EPS*( XY[NM1] - 4.0*XY[NM2] +  6.0*XY[NM3] -  4.0*XY[NM4] +      XY[NM5] )
        RHS[NM2] =  0.01*EPS*( XY[NM1] - 2.0*XY[NM2] +      XY[NM3] )  
        RHS[NM1] =  0.0 
                     
        RHS[5:NM5] = EPS*((XY[0:NM5-5]+XY[10:NM5+5])
                   - 10.0*(XY[1:NM5-4]+XY[ 9:NM5+4])
                   + 45.0*(XY[2:NM5-3]+XY[ 8:NM5+3])
                   -120.0*(XY[3:NM5-2]+XY[ 7:NM5+2])
                   +210.0*(XY[4:NM5-1]+XY[ 6:NM5+1])
                   -252.0* XY[5:NM5]                )
        
        B = scipy.linalg.lu_factor(A.toarray())

        XF[1:-1] = XF[1:-1] + scipy.linalg.lu_solve(B, RHS[1:-1])

#        XF[1:-1] = XF[1:-1] + spsolve(A, RHS[1:-1])

        
        return XF

    #---------------------------------------------------------------------------------------------
    # Code to do 1D or 2D input

    if len(xy2d.shape) < 2:
        A = Filter10_Init(xy2d.shape[0], eps)
        
        print(f'\nCOND:  \n',np.linalg.cond(A.toarray()))


        return Filter1D(xy2d[:], eps, A, **kwargs)
    
    elif len(xy2d.shape) > 2:
        print("RaymondFilter6:  3D filtering not implemented as of yet, exiting\n")
        sys.exit(-1)
        
    else:

        ny, nx = xy2d.shape
        
        print("RaymondFilter10 called:  Shape of array:  NY: %d  NX:  %d" % (ny, nx))

        x1d = np.zeros((nx,))
        y1d = np.zeros((ny,))
        
        XYRES = xy2d.copy()

        tic = time.perf_counter()
        
        A = Filter10_Init(ny, eps)
        for i in np.arange(nx):
            y1d[:]     = xy2d[:,i]
            XYRES[:,i] = Filter1D(y1d, eps, A, **kwargs)
        
        toc = time.perf_counter()
        print(f"I-loop {toc - tic:0.4f} seconds")
    
        tic = time.perf_counter()
        XYRES2 = XYRES.copy()
        A = Filter10_Init(nx, eps)
        for j in np.arange(ny):
            XYRES[j,:] = Filter1D(XYRES2[j,:], eps, A, **kwargs)
        
        toc = time.perf_counter()
        print(f"J-loop {toc - tic:0.4f} seconds")
        
        return XYRES

#---------------------------------------------------------------------------------------------
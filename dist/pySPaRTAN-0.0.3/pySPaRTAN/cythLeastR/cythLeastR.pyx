from __future__ import division
import numpy.matlib as np
import numpy.linalg as npl
cimport numpy as np
include "l1lq.pxi"

cpdef tuple LeastR(double[:,:] _A, double[:,:] _y, double z, dict opts):

    # Verify and initialize the parameters
    cdef np.ndarray A = np.asmatrix(_A)
    cdef np.ndarray y = np.asmatrix(_y)
    cdef np.ndarray x, funVal, ValueL
    cdef np.ndarray ind
    cdef np.ndarray mu, nu
    cdef int m, n
    cdef double q
    cdef double lamda
    cdef long i, j
    

    m = int(A.shape[0])
    n = int(A.shape[1])
 
    if y.size != m:
        raise ValueError('Check the length of y!')

    if z <= 0:
        raise ValueError('z should be positive!')
    
    sll_opts(opts)
    

    # Normalization
    cdef np.ndarray ind_zero
    if opts['nFlag'] != 0:
        if 'mu' in opts:
            mu = opts['mu']
            if mu.shape[1] != n:
                raise ValueError('Check the input .mu')
        else:
            mu = A.mean(axis=0)

        if opts['nFlag'] == 1:
            if 'nu' in opts:
                nu = opts['nu']
                if nu.shape[0] != n:
                    raise ValueError('Check the input .nu')
            else:
                nu = np.power(np.power(A, 2).sum(axis=0)/m, 0.5)
                nu = nu.T
        elif opts['nFlag'] == 2:
            if 'nu' in opts:
                nu = opts['nu']
                if nu.shape[0] != m:
                    raise ValueError('Check the input .nu')
            else:
                nu = np.power(np.power(A, 2).sum(axis=1)/n, 0.5)
        ind_zero = np.absolute(nu) <= 1e-10
        nu[ind_zero] = 1

    # Starting point initialization
    cdef np.ndarray ATy = np.zeros(n)
    cdef np.ndarray ind_i_m
    cdef np.ndarray ind_i
    cdef long m_i
    cdef np.ndarray tt
    cdef np.ndarray invNu
    cdef np.ndarray mu_invNu
    cdef double lamda_max = 0
    cdef np.ndarray Ax = np.zeros(m)
    cdef double x_norm = 0
    cdef double ratio
    
    
    if opts['nFlag'] == 0:
        ATy = A.T * y
    elif opts['nFlag'] == 1:
        ATy = A.T * y - y.sum() * mu.T
        ATy=ATy / nu
    else:
        invNu = y / nu
        ATy= A.T * invNu - invNu.sum() * mu.T
     
    if 'rsL2' in opts:
        rsL2 = opts['rsL2']
        if (rsL2 < 0):
            raise ValueError('\n opts.rsL2 should be nonnegative!')
    else:
        rsL2 = 0

    

    if opts['rFlag'] == 0:
        lamda = z
    else:
        if (z < 0) or (z > 1):
            raise ValueError('opts.rFlag=1, and z should be in [0,1]')

        lambda_max = np.maximum(np.absolute(ATy))
        lamda = z * lambda_max;
    
        rsL2 = rsL2 * lambda_max  # the input rsL2 is a ratio of lambda_max


    if opts['init'] == 2:
        x = np.zeros(n)
    else:
        if 'x0' in opts:
            x = opts['x0']

            if (x.shape[0] != n):
                raise ValueError('Check the input .x0')
        else:
            x = ATy

    
    if opts['nFlag'] == 0:
        Ax = A * x
    elif opts['nFlag'] == 1:
        invNu = x / nu
        mu_invNu = mu * invNu
        Ax = A * invNu - np.repmat(mu_invNu, m, 1)
    else:
        Ax = A * x - np.repmat(mu * x, m, 1)
        Ax = Ax / nu            
    
            
    if opts['init'] == 0:
        x_norm = np.absolute(x).sum() 
        x_2norm = x.T * x
        if x_norm >= 1e-6:
            ratio=initFactor(x_norm, Ax, y, lamda,'LeastR', rsL2, x_2norm);
            x = ratio * x
            Ax = ratio * Ax

    # The main program
    # The Armijo Goldstein line search schemes + accelerated gradient descent
    cdef long bFlag
    cdef double L, beta
    cdef np.ndarray xp
    cdef np.ndarray Axp
    cdef np.ndarray xxp
    cdef long iterStep
    cdef double alphap = 0, alpha = 1
    cdef np.ndarray s
    cdef np.ndarray As
    cdef np.ndarray ATAs
    cdef np.ndarray g
    cdef np.ndarray v
    cdef np.ndarray Av
    cdef double r_sum, l_sum
    cdef np.ndarray Axy
    cdef double norm_xxp, norm_xp
    cdef double gamma
    cdef np.ndarray t
    cdef np.ndarray tp
    cdef np.ndarray ATAx
    cdef np.ndarray ATAxp
    cdef np.ndarray s_t
    cdef np.ndarray xnew, tnew
    cdef np.ndarray Axnew
    cdef np.ndarray u, v_t
    cdef double tao
 
    if (opts['mFlag'] == 0) and (opts['lFlag'] == 0):
        bFlag = 0
        L = 1 + rsL2
        xp = x.copy()
        Axp = Ax.copy()

        xxp = np.zeros((n, 1))
        
        funVal = np.squeeze(np.asarray(np.zeros(opts['maxIter'])))
        ValueL = np.squeeze(np.asarray(np.zeros(opts['maxIter'])))
        for iterStep in xrange(opts['maxIter']):
            # step 1
            beta = (alphap - 1) / alpha
            s = x + beta * xxp

            # step 2
            As = Ax + beta * (Ax - Axp)

            ATAs = np.zeros((n, 1))
            
            if opts['nFlag'] == 0:
                ATAs = A.T * As
            elif opts['nFlag'] == 1:
                ATAs = A.T * As - As.sum() * mu.T  
                ATAs = ATAs / nu
            else:
                invNu = As / nu
                ATAs = A.T * invNu - invNu.sum() * mu.T
                
            g = ATAs - ATy + rsL2 * s

            xp = x.copy()
            Axp = Ax.copy()
            while 1:
                v = s - g / L				
                x = np.multiply(np.sign(v), np.maximum((np.absolute(v)-lamda / L),0))
                v = x - s
              
                if opts['nFlag'] == 0:
                    Ax = A * x
                elif opts['nFlag'] == 1:  
                    invNu = x / nu
                    mu_invNu = mu * invNu
                    Ax = A * invNu - np.repmat(mu_invNu, m, 1)
                else:
                    Ax = A * x - np.repmat(mu * x, m, 1)
                    Ax = Ax / nu
                

                Av = Ax - As

                r_sum = npl.norm(v, 'fro') ** 2
                l_sum = npl.norm(Av, 'fro') ** 2

                if r_sum <= 1e-20:
                    bFlag = 1
                    break

                if(l_sum <= r_sum * (L - rsL2)):
                    break
                else:
                    L = np.maximum(2 * L, l_sum / r_sum  + rsL2)
            ValueL[iterStep] = L
 
            
            # step 3
            alphap = alpha
            alpha = (1 + np.sqrt(4 * alpha * alpha + 1)) / 2
            

            xxp = x - xp
            Axy = Ax - y

            funVal[iterStep] = Axy.T * Axy / 2 + rsL2 / 2 * x.T @ x + np.absolute(x).sum() * lamda

            if bFlag != 0:
                break

            if opts['tFlag'] == 0:

                if iterStep >= 2:
                    if np.absolute(funVal[iterStep] - funVal[iterStep-1]) <= opts['tol']:
                        break
            elif opts['tFlag'] == 1:

                if iterStep >= 2:
                    if np.absolute(funVal[iterStep] - funVal[iterStep-1]) <= (opts['tol'] * funVal[iterStep-1]):
                        break
            elif opts['tFlag'] == 2:
                if funVal[iterStep] <= opts['tol']:
                    break
            elif opts['tFlag'] == 3:
                norm_xxp = npl.norm(xxp, 'fro')
                if norm_xxp <= opts['tol']:
                    break
            elif opts['tFlag'] == 4:
                norm_xp = npl.norm(xp, 'fro')
                norm_xxp = npl.norm(xxp, 'fro')
                if norm_xxp <= (opts['tol'] * np.maximum(norm_xp, 1)):
                    break
            elif opts['tFlag'] == 5:
                if iterStep >= opts['maxIter']:
                    break

    # Reformulated problem + adaptive line search

    if (opts['mFlag'] == 1) and (opts['lFlag'] == 1):

        bFlag = 0
        L = 1 + rsL2
        gamma = 1
        xp = x.copy()
        Axp = Ax.copy()
        xxp = np.zeros((n, 1))
        t = np.absolute(x)
        tp = t.copy()
        As = Ax.copy()

        if opts['nFlag'] == 0:
            ATAx = A.T * As
        elif opts['nFlag'] == 1: 
            ATAx =  A.T * As - As.sum() * mu.T
            ATAx = ATAx / nu
        else:
            invNu = Ax / nu
            ATAx=A.T * invNu - invNu.sum() * mu.T
            
        
        for iterStep in xrange(opts['maxIter']):
            ATAxp = ATAx.copy()
            if iterStep != 1:
                
               # compute AT Ax
                if opts['nFlag'] == 0:
                    ATAx = A.T * Ax;
                elif opts['nFlag'] == 1:
                    ATAx = A.T * Ax - Ax.sum() * mu.T
                    ATAx = ATAx / nu
                else:
                    invNu = Ax / nu
                    ATAx = A.T * invNu - invNu.sum() * mu.T
                

            #--------Line Search for L begins
            Axnew = np.zeros((m, 1))
            while 1:
                if iterStep != 1:
                    alpha = (-gamma + np.sqrt(gamma ** 2 + 4 * L * gamma)) / (2 * L)
                    beta = (gamma - gamma * alphap) / (alphap * gamma + alphap * L * alpha)
                    s = x + beta * xxp
                    s_t = t + beta * (t - tp)
                    As = Ax + beta * (Ax - Axp)
                    ATAs = ATAx + beta * (ATAx - ATAxp)
                else:
                    alpha = (-1 + np.sqrt(5)) / 2
                    beta = 0
                    s = x.copy()
                    s_t = t.copy()
                    As = Ax.copy()
                    ATAs = ATAx.copy()
                    
                g = ATAs - ATy + rsL2 * s
                u = s - g / L
                v = s_t - lamda / L
                xnew, tnew = ep21R(np.asarray(u, order='F'), np.asarray(v, order='F'), n, 1)
                v = xnew - s
                v_t = tnew - s_t

                if opts['nFlag'] == 0:
                    Axnew=A * xnew
                elif opts['nFlag'] == 1:
                    invNu = xnew / nu
                    mu_invNu = mu * invNu
                    Axnew = A * invNu - np.repmat(mu_invNu, m, 1)
                else:
                    Axnew = A * xnew - np.repmat(mu * xnew, m, 1)     
                    Axnew = Axnew / nu
                                
                
                Av = Axnew - As
                r_sum = npl.norm(v, 'fro') ** 2 + np.double(v_t.T * v_t)
                l_sum = npl.norm(Av, 'fro') ** 2 + npl.norm(v, 'fro') ** 2 * rsL2
                

                if r_sum <= 1e-20:
                    bFlag = 1
                    break

                if l_sum <= r_sum * L:
                    break
                else:
                    L = np.maximum(2 * L, l_sum / r_sum)
            #--------Line Search for L ends
            gamma = L * alpha * alpha
            alphap = alpha
            ValueL[iterStep] = L

            tao = L * r_sum / l_sum
            if tao >= 5:
                L = L * 0.8

            xp = x.copy()
            x = xnew.copy()
            xxp = x - xp
            Axp = Ax.copy()
            Ax = Axnew.copy()
            tp = t.copy()
            t = tnew.copy()

            Axy = Ax - y
            funVal[iterStep] = Axy.T * Axy / 2 + rsL2 / 2 * x.T * x + lamda * t.sum()

            if bFlag != 0:
                break

            if opts['tFlag'] == 0:
                 if iterStep >= 2:
                    if np.absolute(funVal[iterStep] - funVal[iterStep-1]) <= opts['tol']:
                        break
            elif opts['tFlag'] == 1:
                if iterStep >= 2:
                    if np.absolute(funVal[iterStep] - funVal[iterStep-1]) <= (opts['tol'] * funVal[iterStep-1]):
                        break
            elif opts['tFlag'] == 2:
                if funVal[iterStep] <= opts['tol']:
                    break
            elif opts['tFlag'] == 3:
                norm_xxp = npl.norm(xxp, 'fro')
                if norm_xxp <= opts['tol']:
                    break
            elif opts['tFlag'] == 4:
                norm_xp = npl.norm(xp, 'fro')
                norm_xxp = npl.norm(xxp, 'fro')
                if norm_xxp <= (opts['tol'] * np.maximum(norm_xp, 1)):
                    break
            elif opts['tFlag'] == 5:
                if iterStep >= opts['maxIter']:
                    break

    # Reformulated problem + Nemirovski's line search
    if (opts['mFlag'] == 1) and (opts['lFlag'] == 0):
        L = 1 + rsL2
        bFlag = 0
        xp = x.copy()
        Axp = Ax.copy()
        xxp = np.zeros((n, 1))

        alphap = 0
        alpha = 1
        
        t = np.absolute(x)
        tp = t.copy()

        for iterStep in xrange(opts['maxIter']):
            # step 1
            beta = (alphap - 1) / alpha
            s = x + beta * xxp
            s_t = t + beta * (t - tp)

            # step 2
            As = Ax + beta * (Ax - Axp)
            ATAs = np.zeros((n, 1))
            
            
            if opts['nFlag'] == 0:
                ATAs = A.T * As
            elif opts['nFlag'] == 1:
                ATAs = A.T * As - As.sum() * mu.T
                ATAs = ATAs / nu
            else:
                invNu = As / nu
                ATAs = A.T * invNu - invNu.sum() * mu.T
           
            g = ATAs - ATy + rsL2 * s
            
            xp = x.copy()
            Axp = Ax.copy()
            tp = t.copy()

            while 1:
                u = s - g / L
                v = s_t - lamda / L
                
                xnew, tnew = ep21R(np.asarray(u ,order='F'), np.asarray(v ,order='F'), n, 1)
                v = xnew - s
                v_t = tnew-s_t
                if opts['nFlag'] == 0:
                    Ax = A * xnew
                elif opts['nFlag'] == 1:
                    invNu = xnew / nu
                    mu_invNu = mu * invNu
                    Ax = A * invNu - np.repmat(mu_invNu, m, 1)
                else:
                    Ax = A * xnew - np.repmat(mu * xnew, m, 1)
                    Ax = Ax / nu
            
                Av = Ax - As
                r_sum = npl.norm(v, 'fro') ** 2 + npl.norm(v_t, 'fro') ** 2
                l_sum = np.double(Av.T * Av) + npl.norm(v, 'fro') ** 2 * rsL2

                if r_sum <= 1e-20:
                    bFlag = 1
                    break

                if l_sum <= r_sum * L:
                    break
                else:
                    L = np.maximum(2 * L, l_sum / r_sum)

            ValueL[iterStep] = L

            # step 3
            alphap = alpha
            alpha = (1 + np.sqrt(4 * alpha ** 2)) / 2

            xxp = x - xp
            Axy = Ax - y
            funVal[iterStep] = Axy.T * Axy / 2 + lamda * t.sum()

            if bFlag != 0:
                break

            if opts['tFlag'] == 0:
                if iterStep >= 1:
                    if np.absolute(funVal[iterStep] - funVal[iterStep-1]) <= opts['tol']:
                        break
            elif opts['tFlag'] == 1:
                if iterStep >= 1:
                    if np.absolute(funVal[iterStep] - funVal[iterStep-1]) <= (opts['tol'] * funVal[iterStep-1]):
                        break
            elif opts['tFlag'] == 2:
                if funVal[iterStep] <= opts['tol']:
                    break
            elif opts['tFlag'] == 3:
                norm_xxp = npl.norm(xxp, 'fro')
                if norm_xxp <= opts['tol']:
                    break
            elif opts['tFlag'] == 4:
                norm_xp = npl.norm(xp, 'fro')
                norm_xxp = npl.norm(xxp, 'fro')
                if norm_xxp <= (opts['tol'] * np.maximum(norm_xp, 1)):
                    break
            elif opts['tFlag'] == 5:
                if iterStep >= opts['maxIter']:
                    break
    if (opts['mFlag'] == 0) and (opts['lFlag'] == 1):
        raise ValueError('The function does not support opts.mFlag=0 & opts.lFlag=1!')

    return (np.asarray(x), np.asarray(funVal))


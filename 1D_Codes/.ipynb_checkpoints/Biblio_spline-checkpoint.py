# -*- coding: utf-8 -*-
from numpy import ones, zeros,linspace,array,meshgrid
import matplotlib.pyplot as plt



def find_span(knots,t,p):
    low  = p
    high = len(knots)-p-1
    mid  =(low+high)//2
    if t <= knots[low ]:return low
    if t >= knots[high ]:return high-1
    assert t > knots[low ], 'out of your starts'
    assert t < knots[high ],'out of your end'
    while t< knots[mid] or t>=knots[mid+1]:
        if t< knots[mid]:
            high=mid
        else:
            low=mid
        mid=(high+low)//2
    return mid
#########################################

def make_knots(grid, p):
    knots = list(grid[0]*ones(p)) + list(grid) + list(grid[-1]*ones(p))
    return knots

#############################################

def basis_funct(t,p,T,i):
    left = [0] * (p+1)
    right = [0] * (p+1)
    N = zeros(p+1) 
    N[0] = 1
    for j in range(1,p+1):
        left[j]=t-T[i+1-j]
        right[j]=T[i+j]-t
        saved=0
        for r in range(j):
            temp=N[r]/(right[r+1]+left[j-r])
            N[r]=saved+right[r+1]*temp
            saved=left[j-r]*temp
        N[j] = saved
    return N
##################################################
def nth_Derv_basis_functs(t, p, n, knots):
    i = find_span(knots,t,p)
    left = [0] * (p+1)
    right = [0] * (p+1)
    ndu  = zeros((p+1,p+1))
    a  = zeros((2,p+1))
    ders = zeros((n+1,p+1))
    ndu[0,0] = 1.0
    #minelem = min(p,n)
    for j in range(1,p+1):
        left[j]=t-knots[i+1-j]# en commence par 0
        right[j]=knots[i+j]-t
        saved = 0.0    
        for r in range(0,j):
            ndu[j,r]=right[r+1]+left[j-r]
            temp=ndu[r,j-1] / ndu[j,r]            
            ndu[r,j]=saved+right[r+1]*temp
            saved=left[j-r]*temp
        ndu[j,j]=saved

    for j in range(p+1):   # Laod the basis fct
        ders[0,j] = ndu[j,p]
        
  # This section Compute  the derivative 
    for r in range(p+1):
        s1 = 0 ; s2 = 1; a[0,0] = 1
    # loop to computz K th derivative 
        for k in range(1,n+1):
            d = 0.0; rk = r-k; pk = p-k
            if(r>=k):
                a[s2,0] = a[s1,0]/ndu[pk+1,rk]
                d       = a[s2,0]*ndu[rk,pk]
            if(rk>=-1):
                j1 = 1
            else:
                j1 = -rk
            if(r-1<=pk):
                j2 = k-1
            else:
                j2 =p-r
            for j in range(j1,j2+1): # PAS CLAIR!!!!!!!  
                a[s2,j]=(a[s1,j]-a[s1,j-1])/ndu[pk+1,rk+j]
                d+= a[s2,j]*ndu[rk+j,pk]
            if(r<=pk):
                a[s2,k] = -a[s1,k-1]/ndu[pk+1,r]
                d += a[s2,k]*ndu[r,pk]
            ders[k,r] = d
            j = s1; s1 = s2; s2 = j
        r = p
    for k in range(1,n+1):
        for j in range(p+1):
            ders[k,j] *= r
        r *= (p-k)
    return ders
######################################## B spline basis 
plt.figure(figsize=(6,3.5))
def plot_bassis_function(knots,degree,nx=101):
    x_min=knots[degree]
    x_max=knots[-degree-1]
    x_s=linspace(x_min,x_max,nx)
    nspline=len(knots)-degree-1# NOMBRE DE FONCTION DE BASE Bspline 
    P=zeros((nx,nspline))
    for i,xi in enumerate(x_s):
        ispan=find_span(knots,xi,degree)
        values_xi=basis_funct(xi,degree,knots,ispan)
        P[i,ispan-degree:ispan+1]+=values_xi
    
    for j in range(nspline):
        plt.plot(x_s,P[:,j],label=r'$N_{}^{}$'.format(j,degree))
        plt.legend(loc='best')
    plt.grid()
    if os.path.exists('Bspline_p=0.png'):
        os.remove('Bspline_p=0.png')
    plt.savefig('Bspline_p=0.png')
    plt.show()


######################################## B spline CUrve 
def plot(knots,degree,alpha,nx):
    t=linspace(knots[degree],knots[-degree-1],nx)
    P=zeros((len(alpha),1))
    P[:,0]= alpha[:]
    Q = zeros((nx,1))
    for i,xi in enumerate(t):
        C = zeros(P.shape[-1])
        ispan=find_span(knots,xi,degree)
        values_xi=basis_funct(xi,degree,knots,ispan)
        for jk in range(degree+1):
            C[:] += values_xi[jk]*P[ispan-degree+jk,:]
        Q[i,:]=C[:]
    return Q[:,0]    









########################################


def Point_on_Rational_Bspline_curve(knots, degree, alpha,  nx, w):
    t = np.linspace(knots[degree], knots[-degree-1], nx)
    P = np.zeros((len(alpha), 1))
    P[:, 0] = alpha[:]
    Q = zeros((nx, 1))
    for i, xi in enumerate(t):
        C = zeros(P.shape[-1])
        ispan = find_span(knots, xi, degree)
        values_xi = basis_funct(xi, degree, knots, ispan)
        ratio = 0.0
        for k in range(0, degree+1):  
            values_xi[k] = values_xi[k] * w[ispan-degree+k] 
            ratio += values_xi[k]
        R = values_xi / ratio
        for jk in range(degree+1):
            C[:] += R[jk] * P[ispan-degree+jk, :]
        Q[i, :] = C[:]
    return Q[:, 0]
    
########################################
    


def plot_surface_2D(knots_1, knots_2, alpha, degree_1, degree_2, nx = 100, ny = 100):
    xs = linspace(knots_1[degree_1],knots_1[-degree_1-1], nx)
    ys = linspace(knots_2[degree_2],knots_2[-degree_2-1], ny)
    n,m,d = alpha.shape
    P = zeros((n,m,d,1))
    Q = zeros((nx,ny,d,1))
    P[:,:,:,0] = alpha.copy()
    saved = zeros(d)
    for i ,xi in enumerate(xs):
        for j,yj in enumerate(ys):
            ispan = find_span(knots_1,xi,degree_1)     
            jspan = find_span(knots_2,yj,degree_2)
            values_xi = basis_funct(xi,degree_1,knots_1,ispan)
            values_yj = basis_funct(yj,degree_2,knots_2,jspan)  
            saved[:]  = 0.0
            for ik in range(degree_1 + 1):
                 for jk in range(degree_2 + 1):
                     index_x = ispan-degree_1+ik
                     index_y = jspan-degree_2+jk
                     saved[:]   += values_xi[ik]*values_yj[jk]*P[index_x, index_y,:,0]
            Q[i,j,:,0] = saved[:]
    return  Q[:,:,:,0]     



#########################################################"


def plot_R_surface_2D(knots_1, knots_2, alpha, degree_1, degree_2,w1,w2, nx = 100, ny = 100):
    xs = linspace(knots_1[degree_1],knots_1[-degree_1-1], nx)
    ys = linspace(knots_2[degree_2],knots_2[-degree_2-1], ny)
    n,m, d = alpha.shape
    P = zeros((n,m,d,1))
    Q = zeros((nx,ny,d,1))
    P[:,:,:,0] = alpha.copy()
    saved = zeros(d)
    for i ,xi in enumerate(xs):
        for j,yj in enumerate(ys):
            ispan = find_span(knots_1,xi,degree_1)     
            jspan = find_span(knots_2,yj,degree_2)
            values_xi = basis_funct(xi,degree_1,knots_1,ispan)
            values_yj = basis_funct(yj,degree_2,knots_2,jspan)  
            saved[:]  = 0.0
            ratio = 0.0
            for ik in range(degree_1 + 1):
                for jk in range(degree_2 + 1):
                    R  = values_xi[ik] * w1[ispan-degree_1+ik] 
                    R *= values_yj[jk] * w2[jspan-degree_2+jk]
                    ratio += R
                    saved[:]   += R*P[ispan-degree_1+ik,jspan-degree_2+jk,:,0]
                     # saved[:]   += values_xi[ik]*values_yj[jk]*P[ispan-degree_1+ik,jspan-degree_2+jk,:,0]

            
            Q[i,j,:,0] = saved[:]/ratio
    return  Q[:,:,:,0]























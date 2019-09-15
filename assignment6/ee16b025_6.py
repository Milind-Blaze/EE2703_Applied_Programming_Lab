
# coding: utf-8

# # Abstract 
# 
# Tubelights work by the emission of electrons at the cathode that undergo accelecration due to an applied electric field and excite atoms which emit photons during subsequent relaxation. Thus, models of tubelights can be built efficiently and accurately using Python with the aid of libraries such as numpy and scipy that enable one to visualise the working of this device using intensity plots and phase spaces.
# 
# # Introduction 
# 
# The tubelight is modeled as a one dimensional two terminal device with the ends being the cathode and the anode. Electrons are emitted at the anode and accelerate under an applied uniform electric field $E_0$ with an acceleration fo $1ms^{-2}$. When they reach a critical velocity $u_0$ they can undergo fruitful collisions with atoms that get excited and relax immediately to emit photons. Electrons that reach the anode are lost. This model is simulated for $nk$ turns, each beginning with an injection of electrons. The injection of electrons is modeled as follows as a normally distributed random variable.
# \begin{equation}
# m = N + Msig \cdot X
# \end{equation}
# where X is a normally distributed random variable. Thus 
# \begin{align}
# E(m) &=N + Msig \cdot E(X) = N\\
# Var(m) &= Msig^2Var(X)\\
# \end{align}
# Also $E(X)$ is chosen to be $0$. The integer part of $m$ is chosen to be the number of electrons injected. Msig determines the variance of the random variable m.
# 
# # Code and results
# 
# The necessary libraries are imported and the default size for images is set.

# In[36]:




from __future__ import division
import matplotlib.pyplot as plt
from matplotlib import cm, colors
import numpy as np
import sys
size=(10,8)


# In[15]:


# defining the necessary functions


# The default parameters are set as given below. Alternately the same can be taken from commandline arguments. The default length of the tube is taken to be $n=100$ units. Further, the probability that a sufficiently energitic electron undergoes a collision in a turn is given by $p$. 

# In[16]:


# defining the constants 

# defaults

n= 100
M= 5
nk= 500
u0= 7
p= 0.5
Msig=2


# In[4]:


# # considering the command line arguments
if len(sys.argv)==7:
    n= float(sys.argv[1])
    M= float(sys.argv[2])
    nk= int(sys.argv[3])
    u0= float(sys.argv[4])
    p= float(sys.argv[5])
  #  if p==1:
   #     print("This is a pathological case, try 0.99")
    #    sys.exit()
    Msig= float(sys.argv[6])

if (len(sys.argv)>1 and len(sys.argv)!=7):
    print ("\nCommand line arguments to be given as follows:\n  $ python ee16b025_6.py n M nk u0 p Msig \nwhere n, M nk, u0, p, Msig represent the same things as given in the question.\n")
    sys.exit()

# # Edit this
# The length list, in which each element represents one electron, is made sufficiently large. Other arrays are defined as follows
# \begin{array}{c c c}
# xx &-& electron position
# u &-& electron velocity
# dx &-& electron displacement in one turn
# I &-& list of positions of every photon ever emitted
# X &-& list of positions of electrons that existed at the end of every turn
# V &-& the electron velocities corresponding to X
# 
# \end{array}

# In[23]:


length= n*M # shouldn't it be nk*M?

# electron information

xx= np.zeros((length))   # electron position
u= np.zeros((length))    # electron velocity
dx= np.zeros((length))   # disp. in current turn

# extra info i don't really understand why i'm defining

I=[]
X=[]
V=[]
            


# The code does the following-
#  - ii finds the positions of electrons in the array xx which exist. $0$ in the xx array indicates the non-existence of an electron. Each existing electron undergoes motion according to the equations
# \begin{align}
# dx &= u_0t + \frac{1}{2}at^2\\
# u &= u_0 + at\\
# \end{align}
#  - Due to this motion, some electrons might reach the anode and get absorbed. Correspondingly, their xx, u, dx values are set to zero. Which basically translates to their disappearance from our simulation universe.
#  - Some electons obtain critical velocity $u_0$. Of these a few undergo collisions leading to immediate photon emission. This is modeled as a uniform random variable. Each energitic electron is assigned a random number (following a uniform distribution) and each electron collides if this number is lesser than $p$ which means the probability of a collsiion is $p$.
#  - The velocities of these colliding electrons are set to zero as they lose all their energy in the collision. Since the collision might occur at any position betweet $x_i$ and $x_{i+1}$ (where $x_i$ is the position of some colliding electron corresponding to the $i^{th}$ turn) the new position is updated by subtracting some random fraction of the displacement that had been added to the position of the colliding electron in that turn.
#  - The positions where the collisions occur are appended to the list $I$. 
#  - New electrons are generated as explained in Section 1 and fill the empty positionsin the xx lists.

# In[28]:

ii= np.where(xx>0)
for _ in range(0,nk):
# to find where electrons are active
    dx[ii]= u[ii]+0.5
    xx[ii]= xx[ii]+dx[ii]
    u[ii]= u[ii]+1
    anode= np.where(xx>n)
    xx[anode]= 0
    dx[anode]= 0
    u[anode]= 0
    kk= np.where(u>u0) # electrons with energy
    ll= np.where(np.random.random(len(kk[0]))<=p)
    kl= np.array(kk)[0][ll]  # which electrons will ionize
    u[kl]= 0
    rho= np.random.random(1)[0]
    xx[kl]= xx[kl]- dx[kl]*rho # minus cos xx has already been updated with an addition of dx in the loop
    I.extend(xx[kl])
    m= int(M + Msig*np.random.randn()) # mean  5 standard deviation 2
    new= np.where(xx==0)
    xx[new[0][:m]]=1
    ii= np.where(xx>0)
    X.extend(xx[ii])
    V.extend(u[ii])
  
    
    

    


# In[51]:


plt.figure(0,figsize=size)
plt.title("Electron density",fontsize=18)
plt.xlabel("Position",fontsize=18)
plt.ylabel("Number of electrons to have ever occupied the bin",fontsize=18)
plt.grid(True)
plt.hist(X,np.arange(0,n+1,1),color="white")
plt.show()    
plt.close()


# The following histogram shows I vs position. This is a representation of the average value of the number of photons emitted between any two consecutive integral positions along hte tubelight and thus a representation of the average intensity.

# In[50]:

if p!=0:
    plt.figure(1,figsize=size)
    plt.title("Intensity (average) of the tubelight vs position", fontsize=18)
    plt.xlabel("Position", fontsize=18)
    plt.ylabel("Intensity/ Number of photons emitted",fontsize=18)
    plt.grid(True)
    ret=plt.hist(I,np.arange(0,n+1,1),color="white")
    vals=ret[0]
    vals=1-(vals/max(vals))
    norm = colors.Normalize(vals.min(), vals.max())
    for thisfrac, thispatch in zip(vals, ret[2]):
        color = cm.Greys(norm(thisfrac))
        thispatch.set_facecolor(color)
    plt.show()    
    plt.close()  


# We plot the phase space of the electrons and also make a table of the centers of the bins and the intensity (number of photons) values.

# In[53]:


plt.figure(2,figsize=size)
plt.xlabel("X",fontsize=18)
plt.title("Electron phase space",fontsize=18)
plt.ylabel("V",fontsize=18)
plt.grid(True)
plt.plot(X,V,"bx")
plt.show() 
plt.close()


# The following is a 2D extension of the one dimensional model wherein the intensity value along $y$ axis has been kept the same at a given $x$.

# In[87]:

if p!=0:
    intensity=np.array([ret[0]]*int(n))
    y,x=np.meshgrid(ret[1][:-1],ret[1][:-1])
    plt.figure(figsize=(14,4))
    plt.contourf(y,x,-intensity,cmap=cm.Blues)
    plt.title("Tubelight, extension from 1D to 2D", fontsize=18)
    plt.xlabel("x", fontsize=18)
    plt.show()
    plt.close()


# In[9]:

if p!=0:
    bins=ret[1]
    bins= 0.5*(bins[0:-1]+bins[1:])
    print ("Intensity data \n")
    print ("xpos     count")
    for i in range(len(bins)):
        if bins[i]<10:
            print str(bins[i])+"       "+str(ret[0][i])
        else:
            print str(bins[i])+"      "+str(ret[0][i])

if p==0:
    print "Nothing to show for the intensity plots and table as no collisions occur with p=0"

# # Discussion and Conclusion-
# 
# The population plot of Figure 1 indicates the number of electrons that have been at a given position during the simulation. Since all the electrons start at the same point the first bin has the highest peak. This is followed by $x\approx25$ which is where most electrons have the energy necessary for ionisation. 
# 
# However the plot of Figure 2 is of far greater interest. It is effectively a plot of the number of electron atom collisions or the number of photons emitted at a given point or in a given bin. Since most of the electrons have insufficent energy for the first 20 - 25 units or so, there are no bars in this region. However most of the electrons gain sufficient energy by 25 (for the given defaults), and thus cause photon emission. Thus there is a huge peak at this point. Electrons that don't suffer collisions move on and perhaps do so at later stages contributing to the further bars. The next set of peaks is noticed at around 60, which incidentally happens to be nearly twice the initial value of position of the first peak. This could probably be due to the face that  a majority of the electrons that suffered  collisions at the position of the first peak have gained enought energy for a second set of collisions.
# 
# The phase space further adds insight. The phase space is the set of all possible states of a system. In this case, it is evident that no high velocity, low position states are present. This arises from the fact that for all low valus of position, until the first peak in intensity, the velocity of the electrons is only increasing and there are no collisions. Once the electrons have sufficient energy, they begin to suffer collisions which means they end up with zero velocity at multiple values of X and similarly at higher values of velocity for other positions. This means that a multitude of combinations of $X_i$ and $V_i$ are possible which is made obvious by the phase space plot.
# 
# To conclude, the various scientific python libraries have been made use of effectively to model what is a fairly complex situation which inherently involves some randomness. Further the intensity distribution of the tube is studied and found to match expectations.

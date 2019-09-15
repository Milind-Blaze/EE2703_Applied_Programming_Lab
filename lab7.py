
# coding: utf-8

# # Abstract 
# 
# This work pursues the use of scientific python, particularly the use of the scipy library to study the transfer functions of a few systems. Bode plots and time domain responses of the system are obtained and studied using inbuilt tools.

# # Introduction
# 
# The report pursues broadly the following three problems.
# - The use of bode plots to understand the behaviour of systems whose responses are known in the laplace domain. Here, the scipy.signal.impulse function is used.
# - The use of the laplace transforms to analyse the given systems by determining their transfer functions and plotting the corresponding bode plots. This discussion primarily revolves around a second order low pass filter. By applying inputs of multiple frequencies, the behaviour of this system is studied and understood.

# # Method and code-
# 
# The necessary libraries are imported in the following piece of code. A helper function is also defined. This is essentially the unit step function $u(t)$.
# 

# In[77]:


from __future__ import division
get_ipython().magic(u'matplotlib inline')
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sp
size=(10,8)


# In[78]:


# defining useful functions
def u(t):
    return 1*(t>=0)


# ## Problem1
# 
# A system modeled by the following differential equation is studied here. 
# \begin{equation}\label{eqn.diff}
# \ddot{x} + 2.25x = f(t)
# \end{equation}
# where $f(t)$ is given by 
# \begin{equation}
# f(t) = cos(bt)e^{-at}u(t)
# \end{equation} whose laplace transform is given by
# \begin{equation}
# F(s) = \frac{s+a}{(s+a)^2 + b^2}
# \end{equation}. Here we choose $a = 0.5$ and $b = 1.5$. $x$ is the time domain response of the system. The following initial conditions are given to us $\dot{x}(0)=0$ and $x(0)=0$. Switching to the laplace domain the differential equation leads to the following
# \begin{equation}
# X(s) = \frac{s + 0.5}{s^4 + s^3 + 4.75s^2 + 2.25s + 5.625}
# \end{equation}. This transfer function is passed to the $sp.impulse$ function whose job is to essential convert this and the initial conditions into the time domain response. 

# In[79]:


# Problem 1
den=[1,1,4.75,2.25,5.625]
num=[1,0.5]
X= sp.lti(num,den)
T= np.linspace(0,50,301)
t,y=sp.impulse(X,None,T)
plt.plot(t,y)


# In[53]:


w1,S1,phi1=X.bode()
plt.semilogx(w1,S1)


# In[54]:


num= [1,0.05]
den=[1,0.1,4.502,0.225,5.068]
X= sp.lti(num,den)
T= np.linspace(0,50,101)
t,y=sp.impulse(X,None,T)
plt.plot(t,y)


# In[55]:


w2,S2,phi2=X.bode()
plt.semilogx(w2,S2,"b-")
plt.semilogx(w1,S1,"g-")
# note the higher quality factor


# In[56]:


H=sp.lti([1],[1,0,2.25])
o=[]
t=np.linspace(0,100,201)

for i in np.arange(1.4,1.65,0.05):
    ip= (np.cos(i*t))*(np.exp(-0.05*t))*u(t)
    T,op,svec= sp.lsim(H,ip,t)
    plt.figure(figsize=size)   
    plt.plot(T,op)
    plt.show()
    plt.close()


# In[57]:


num= [1,0,2]
den=[1,0,3,0]
X= sp.lti(num,den)
T= np.linspace(0,50,201)
t,y=sp.impulse2(X,None,T)
plt.plot(t,y)


# In[61]:


num= [2]
den=[1,0,3,0]
X= sp.lti(num,den)
T= np.linspace(0,50,201)
t,y=sp.impulse2(X,None,T)
plt.plot(t,y)


# In[63]:


den=[10**(-12),10**(-4),1]
num=[1]
X= sp.lti(num,den)
w,S,phi=X.bode()
plt.subplot(211)
plt.semilogx(w,S,"b-")
plt.subplot(212)
plt.semilogx(w,phi,"g-")
# note the higher quality factor


# In[69]:


t=np.arange(0,10**(-1), 10**(-7))
f= (np.cos(1000*t)-np.cos((10**(6))*t))*u(t)
T,op,svec= sp.lsim(X,f,t)


# In[73]:


plt.figure(figsize=size)   
plt.semilogx(T,op)
#plt.plot(T[0:300],op[0:300])
plt.show()
plt.close()


# In[80]:


conda install notebook jupyterlab # for updating jupyter notebook and lab


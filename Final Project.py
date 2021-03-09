#!/usr/bin/env python
# coding: utf-8

# Hermite's Equation (1) $$y''-2xy'+2py=0$$
# Legendre's Equation (2) $$(1-x^2)y''-2xy'+p(p+1)y=0$$
# Chebyshev's Equation (3) $$(1-x^2)y''-xy'+p^2y=0$$

# All of these equations have ordinary points at $x_o=0$ and have power series of the form $$y=\Sigma_0^\infty{a_nx^n}$$
# 
# To show this, the following are the needed derivatives of y $$y'=\Sigma_1^\infty{na_nx^{n-1}}=\Sigma_0^\infty{(n+1)a_{n+1}x^n}$$ 
# $$y''=\Sigma_2^\infty{n(n-1)a_nx^{n-2}}=\Sigma_1^\infty{n(n+1)a_{n+1}x^{n-1}}=\Sigma_0^\infty{(n+2)(n-1)a_{n+2}x^{n}}$$ 
# 
# We will designate the difference in these by using yd1, yd0, ydd2, ydd1, ydd0 based on their starting points, respectively.

# 1. By first setting $a_0=1$ and $a_1=0$ and then setting $a_0=0$ and $a_1=1$, and a pair of independent solutions to equations (1), (2), and (3) for $p=0,1,2,3$.

# In[201]:


x=var('x')
an,n=var('an,n')
an1=var('an1') #this is a_n+1
an2=var('an2')# this is a_n+2

y=an*x^n
yd1=n*an*x^(n-1)
yd0=(n+1)*an1*x^n
ydd2=n*(n-1)*an*x^(n-2)
ydd1=n*(n+1)*an1*x^(n-1)
ydd0=(n+2)*(n+1)*an2*x^n


# Equation (1): With $a_0=1$, $a_1=0$, and $p=0$. So $$y''-2xy'=0$$

# In[240]:


solve([ydd0-2*x*yd1==0],an2)


# In[109]:


#a0=1,a1=0,p=0
def coef(an,n):
    return(2*an*n/(n^2 + 3*n + 2))
a_nlist=[1,0] #a_0=1, a_1=0
for i in range (0,10):
    new=coef(a_nlist[i],i)
    a_nlist.append(new)
    
show(a_nlist)


#    So based on the a_nlist above, the solution is $y=1$ when $a_0=1$, $a_1=0$.

# Equation (1): With $a_0=0$, $a_1=1$, and $p=0$. So $y''-2xy'=0$

# In[111]:


#a0=0,a1=1,p=0
def coef(an,n):
    return((2*an*n/(n^2 + 3*n + 2)))
a_nlist=[0,1] #a_0=0, a_1=1
for i in range (0,10):
    new=coef(a_nlist[i],i)
    a_nlist.append(new)
    
show(a_nlist)


# To find how the series was made I completely factored each denominator and found the commonalities between them.
# Based on the above, the series solution was found to be $$y=\Sigma_0^\infty{\frac{x^{2n+1}}{n!(2n+1)}}$$

# Equation (1): With $a_0=1$, $a_1=0$, and $p=1$. So $$y''-2xy'+2(1)y=0$$

# In[124]:


n=var('n')
an2=var('an2')

solve([ydd0-2*x*yd1+2*1*y],an2)


# In[118]:


#a0=1,a1=0,p=1
def coef(an,n):
    return(2*(an*n - an)/(n^2 + 3*n + 2))
a_nlist=[1,0] #a_0=1, a_1=0
for i in range (0,10):
    new=coef(a_nlist[i],i)
    a_nlist.append(new)
    
show(a_nlist)


# To find how the series was made I completely factored each denominator and found the commonalities between them.
# 
# Based on the above the series solution was found to be $$1+\Sigma_1^\infty{\frac{-(x^{2n})}{n!(2n-1)}}$$

# In[105]:


#a0=0,a1=1,p=1
def coef(an,n):
    return(2*(an*n - an)/(n^2 + 3*n + 2))
a_nlist=[0,1] #a_0=0, a_1=1
for i in range (0,10):
    new=coef(a_nlist[i],i)
    a_nlist.append(new)
    
show(a_nlist)


# Based on the above the series solution was found to be $$y=x$$

# Equation (1): With $a_0=1$, $a_1=0$, and $p=2$. So $$y''-2xy'+2(2)y=0$$

# In[125]:


n=var('n')
an2=var('an2')

solve([ydd0-2*x*yd1+2*2*y],an2)


# In[128]:


#a0=1,a1=0,p=2
def coef(an,n):
    return(2*(an*n - 2*an)/(n^2 + 3*n + 2))
a_nlist=[1,0] #a_0=1, a_1=0
for i in range (0,10):
    new=coef(a_nlist[i],i)
    a_nlist.append(new)
    
show(a_nlist)


# Based on the above the series solution was found to be $$y=-2x^2+1$$

# Equation (1): With $a_0=0$, $a_1=1$, and $p=2$. So $$y''-2xy'+2(2)y=0$$

# In[130]:


#a0=0,a1=1,p=2
def coef(an,n):
    return(2*(an*n - 2*an)/(n^2 + 3*n + 2))
a_nlist=[0,1] #a_0=0, a_1=1
for i in range (0,10):
    new=coef(a_nlist[i],i)
    a_nlist.append(new)
    
show(a_nlist)


# To find how the series was made I completely factored each denominator and found the commonalities between them.
# 
# Based on the above, the series solution is $$y=\Sigma_0^\infty{\frac{-(x^{2n+1})}{(2n+1)(2n-1)n!}}$$

# Equation (1): With $a_0=1$, $a_1=0$, and $p=3$. So $$y''-2xy'2(3)y=0$$

# In[126]:


n=var('n')
an2=var('an2')

solve([ydd0-2*x*yd1+2*3*y],an2)


# In[137]:


#a0=1,a1=0,p=3
def coef(an,n):
    return(2*(an*n - 3*an)/(n^2 + 3*n + 2))
a_nlist=[1,0] #a_0=1, a_1=0
for i in range (0,10):
    new=coef(a_nlist[i],i)
    a_nlist.append(new)
    
show(a_nlist)


# To find how the series was made I completely factored each denominator and found the commonalities between them.
# 
# Based on the above, the series solution is $$y=\Sigma_0^\infty{\frac{3(x^{2n})}{(2n-3)(2n-1)n!}}$$

# In[138]:


#a0=0,a1=1,p=3
def coef(an,n):
    return(2*(an*n - 3*an)/(n^2 + 3*n + 2))
a_nlist=[0,1] #a_0=0, a_1=1
for i in range (0,10):
    new=coef(a_nlist[i],i)
    a_nlist.append(new)
    
show(a_nlist)


# This solution is found to be $$y=x-\frac{2}{3}x^3$$
# 
# This answer concludes Equation (1)

# Equation (2): $$(1-x^2)y''-2xy'+p(p+1)y=0$$ With $a_0=1$, $a_1=0$, and $p=0$. So $$(1-x^2)y''-2xy'=0$$ This is also: $$y''-x^2y''-2xy'=0$$

# In[176]:


n=var('n')
an2=var('an2')

solve([(ydd0-x^2*ydd2-2*x*yd1)],an2)


# In[140]:


#a0=1,a1=0,p=0
def coef(an,n):
    return(an*n/(n + 2))
a_nlist=[1,0] #a_0=1, a_1=0
for i in range (0,10):
    new=coef(a_nlist[i],i)
    a_nlist.append(new)
    
show(a_nlist)


# This solution is $$y=1$$

# In[153]:


#a0=0,a1=1,p=0
def coef(an,n):
    return(an*n/(n + 2))
a_nlist=[0,1] #a_0=1, a_1=0
for i in range (0,10):
    new=coef(a_nlist[i],i)
    a_nlist.append(new)
    
show(a_nlist)


# This solution is straight forward: $$y=\Sigma_0^\infty{\frac{x^{2n+1}}{2n+1}}$$

# With $a_0=1$, $a_1=0$, and $p=1$. So $$y''-x^2y''-2xy'+1(1+1)y=0$$

# In[149]:


solve([(ydd0-x^2*ydd2-2*x*yd1+1*(1+1)*y)],an2)


# In[148]:


#a0=1,a1=0,p=1
def coef(an,n):
    return((an*n - an)/(n + 1))
a_nlist=[1,0] #a_0=1, a_1=0
for i in range (0,10):
    new=coef(a_nlist[i],i)
    a_nlist.append(new)
    
show(a_nlist)


# This solution is again, straight forward: $$y=\Sigma_0^\infty{\frac{-(x^{2n})}{2n-1}}$$

# In[146]:


#a0=0,a1=1,p=1
def coef(an,n):
    return((an*n - an)/(n + 1))
a_nlist=[0,1] #a_0=1, a_1=0
for i in range (0,10):
    new=coef(a_nlist[i],i)
    a_nlist.append(new)
    
show(a_nlist)


# This solution is $$y=x$$

# With $a_0=1$, $a_1=0$, and $p=2$. So $$y''-x^2y''-2xy'+2(2+1)y=0$$

# In[150]:


solve([(ydd0-x^2*ydd2-2*x*yd1+2*(2+1)*y)],an2)


# In[151]:


#a0=1,a1=0,p=2
def coef(an,n):
    return((an*n^2 + an*n - 6*an)/(n^2 + 3*n + 2))
a_nlist=[1,0] #a_0=1, a_1=0
for i in range (0,10):
    new=coef(a_nlist[i],i)
    a_nlist.append(new)
    
show(a_nlist)


# This solution is $$y=-3x^2+1$$

# In[154]:


#a0=0,a1=1,p=2
def coef(an,n):
    return((an*n^2 + an*n - 6*an)/(n^2 + 3*n + 2))
a_nlist=[0,1] #a_0=1, a_1=0
for i in range (0,10):
    new=coef(a_nlist[i],i)
    a_nlist.append(new)
    
show(a_nlist)


# This solution is more complex but when each piece is written out the solution was found to be:$$y=\Sigma {x+\frac{-4}{3!}x^3+\frac{-24}{5!}x^5+\frac{-576}{7!}x^7+\frac{-28800}{9!}x^9+\frac{-2419200}{11!}x^{11}+...}$$
# 

# With $a_0=1$, $a_1=0$, and $p=3$. So $$y''-x^2y''-2xy'+3(3+1)y=0$$

# In[155]:


solve([(ydd0-x^2*ydd2-2*x*yd1+3*(3+1)*y)],an2)


# In[159]:


#a0=1,a1=0,p=3
def coef(an,n):
    return((an*n^2 + an*n - 12*an)/(n^2 + 3*n + 2))
a_nlist=[1,0] #a_0=1, a_1=0
for i in range (0,10):
    new=coef(a_nlist[i],i)
    a_nlist.append(new)
    
show(a_nlist)


# This solution is more complex but when each piece is written out the solution was found to be: $$y=\sum_0^\infty {1+\frac{-12}{2!}x^2+\frac{72}{4!}x^4+\frac{576}{6!}x^6+\frac{17280}{8!}x^8+\frac{1036800}{10!}x^{10}+...}$$
# 

# In[158]:


#a0=0,a1=1,p=3
def coef(an,n):
    return((an*n^2 + an*n - 12*an)/(n^2 + 3*n + 2))
a_nlist=[0,1] #a_0=1, a_1=0
for i in range (0,10):
    new=coef(a_nlist[i],i)
    a_nlist.append(new)
    
show(a_nlist)


# This is another straight forward solution with $$y=x-\frac{5}{3}x^3$$
# 
# This answer concludes Equation (2).

# Equation (3): $$(1-x^2)y''-xy'+p^2y=0$$ With $a_0=1$, $a_1=0$, and $p=0$. So $$(1-x^2)y''-xy'=0$$ This is also: $$y''-x^2y''-xy'=0$$

# In[160]:


solve([ydd0-x^2*ydd2-x*yd1],an2)


# In[162]:


#a0=1,a1=0,p=0
def coef(an,n):
    return(an*n^2/(n^2 + 3*n + 2))
a_nlist=[1,0] #a_0=1, a_1=0
for i in range (0,10):
    new=coef(a_nlist[i],i)
    a_nlist.append(new)
    
show(a_nlist)


# This solution is $$y=1$$

# In[163]:


#a0=0,a1=1,p=0
def coef(an,n):
    return(an*n^2/(n^2 + 3*n + 2))
a_nlist=[0,1] #a_0=1, a_1=0
for i in range (0,10):
    new=coef(a_nlist[i],i)
    a_nlist.append(new)
    
show(a_nlist)


# This solution was found by factoring all of the parts as well as solving for each a without simplifying: $$y=\Sigma_0^\infty{\frac{(2n-1^{2})!!}{(2n-1)!}x^{2n+1}}$$

# With $a_0=1$, $a_1=0$, and $p=1$: $$y''-x^2y''-xy'+1^2y=0$$

# In[164]:


solve([ydd0-x^2*ydd2-x*yd1+1^2*y],an2)


# In[167]:


#a0=1,a1=0,p=1
def coef(an,n):
    return((an*n - an)/(n + 2))
a_nlist=[1,0] #a_0=1, a_1=0
for i in range (0,10):
    new=coef(a_nlist[i],i)
    a_nlist.append(new)
    
show(a_nlist)


# This solution is: $$y=1+\sum_1^\infty{\frac{-(2n-3)!!}{(2n)!!}x^{2n}}$$

# In[166]:


#a0=0,a1=1,p=1
def coef(an,n):
    return((an*n - an)/(n + 2))
a_nlist=[0,1] #a_0=1, a_1=0
for i in range (0,10):
    new=coef(a_nlist[i],i)
    a_nlist.append(new)
    
show(a_nlist)


# This is again, a straight forward solution with $$y=x$$

# With $a_0=1$, $a_1=0$, and $p=2$: $$y''-x^2y''-xy'+2^2y=0$$

# In[168]:


solve([ydd0-x^2*ydd2-x*yd1+2^2*y],an2)


# In[172]:


#a0=1,a1=0,p=2
def coef(an,n):
    return((an*n - 2*an)/(n + 1))
a_nlist=[1,0] #a_0=1, a_1=0
for i in range (0,10):
    new=coef(a_nlist[i],i)
    a_nlist.append(new)
    
show(a_nlist)


# $$y=-2x^2+1$$

# In[174]:


#a0=0,a1=1,p=2
def coef(an,n):
    return((an*n - 2*an)/(n + 1))
a_nlist=[0,1] #a_0=1, a_1=0
for i in range (0,10):
    new=coef(a_nlist[i],i)
    a_nlist.append(new)
    
show(a_nlist)


# $$y=x+\sum_1^\infty{\frac{-(2n-3)!!}{(2n)!!}x^{2n+1}}$$

# With $a_0=1$, $a_1=0$, and $p=3$: $$y''-x^2y''-xy'+3^2y=0$$

# In[177]:


solve([ydd0-x^2*ydd2-x*yd1+3^2*y],an2)


# In[179]:


#a0=1,a1=0,p=3
def coef(an,n):
    return((an*n^2 - 9*an)/(n^2 + 3*n + 2))
a_nlist=[1,0] #a_0=1, a_1=0
for i in range (0,10):
    new=coef(a_nlist[i],i)
    a_nlist.append(new)
    
show(a_nlist)


# $$y=1+\Sigma_1^\infty{-\frac{9}{2!}x^2+\frac{45}{4!}x^4+\frac{315}{6!}x^6+\frac{8505}{8!}x^8+...}$$

# In[180]:


#a0=0,a1=1,p=3
def coef(an,n):
    return((an*n^2 - 9*an)/(n^2 + 3*n + 2))
a_nlist=[0,1] #a_0=1, a_1=0
for i in range (0,10):
    new=coef(a_nlist[i],i)
    a_nlist.append(new)
    
show(a_nlist)


# $$y=x-\frac{4}{3}x^3$$
# 
# This solution concludes Problem #1.

# 2. For $p=0,1,2,3,4,5,6$, Find polynomial solutions for Equation (2) with $y(1)=1$.Argue that, for each, $p$ there is only one such solution.  These polynomials are called the Legendre polynomials.

# Legendre's Equation (2) $$(1-x^2)y''-2xy'+p(p+1)y=0$$

# When 
# 1. $p=0$, $y''-x^2y''-2xy'=0$
# 2. $p=1$, $y''-x^2y''-2xy'+2y=0$
# 3. $p=2$, $y''-x^2y''-2xy'+2(3)y=0$
# 4. $p=3$, $y''-x^2y''-2xy'+3(4)y=0$
# 5. $p=4$, $y''-x^2y''-2xy'+4(5)y=0$
# 6. $p=5$, $y''-x^2y''-2xy'+5(6)y=0$
# 7. $p=6$, $y''-x^2y''-2xy'+6(7)y=0$

# In[203]:


p=var('p')
solve([ydd0-x^2*ydd2-2*x*yd1+p*(p+1)*y],an2)


# In[200]:


show(an2==(an*n^2 - an*p^2 + an*n - an*p)/(n^2 + 3*n + 2))


# This is also $$a_{n+2}=\frac{a_n(n^2-p^2+n-p)}{(n+1)(n+2)}$$
# Which will simplify to $$a_{n+2}=\frac{a_n(n-p)(n+p+1)}{(n+1)(n+2)}$$

# Because we are given $y(1)=1$, that will make it easier to find the solutions that are correct from the #1-#6 with $a_0=1,a_1=0$ and $a_0=0,a_1=1$.
# 
# There will only be only solution equation for each p value because 
# 
# The polynomials for $p=0,1,2,3$ have already been computed above and we will reference those as we complete this question.

# From above we show that for $a_0=1,a_1=0$ and $p=0$, $y=a_0=1$. We consider this to be the Legendre polynomial because it terminates. If the series doesn't terminate, it can't be used as the wanted polynomial. 
# 
# So, based on the parts that have been completed above and we are going to but them into a set thing with $P_p(x)=solution$.
# 
# Looking at the pattern above we can see that when p is even, the even series terminates, and when p is odd, the odd series terminates.
# 
# The full solution will be listed together when all have been found at the end of this section.

# Due to the fact that we have the initial value $y(1)=1$ 
# we are going to impliment an $a_0$ to show this. 

# This will be visible at the end of this section

# In[205]:


#a0=1,a1=0,p=4
solve([(ydd0-x^2*ydd2-2*x*yd1+4*(4+1)*y)],an2)


# In[207]:


def coef(an,n):
    return((an*n^2 + an*n - 20*an)/(n^2 + 3*n + 2))
a_nlist=[1,0] #a_0=1, a_1=0
for i in range (0,10):
    new=coef(a_nlist[i],i)
    a_nlist.append(new)
    
show(a_nlist)


# $$y=\frac{35}{3}x^4-10x^2+1$$

# In[209]:


#a0=0,a1=1,p=5
solve([(ydd0-x^2*ydd2-2*x*yd1+5*(5+1)*y)],an2)


# In[210]:


def coef(an,n):
    return((an*n^2 + an*n - 30*an)/(n^2 + 3*n + 2))
a_nlist=[0,1] 
for i in range (0,10):
    new=coef(a_nlist[i],i)
    a_nlist.append(new)
    
show(a_nlist)


# $$y=\frac{21}{5}x^5-\frac{14}{3}x^3+x$$

# In[211]:


#a0=1,a1=0,p=6
solve([(ydd0-x^2*ydd2-2*x*yd1+6*(6+1)*y)],an2)


# In[213]:


def coef(an,n):
    return((an*n^2 + an*n - 42*an)/(n^2 + 3*n + 2))
a_nlist=[1,0] #a_0=1, a_1=0
for i in range (0,10):
    new=coef(a_nlist[i],i)
    a_nlist.append(new)
    
show(a_nlist)


# $$y=\frac{-231}{5}x^6+63x^4-21x^2+1$$

# 1. $p=0, a_0=1,a_1=0 \rightarrow P_0=1a_0$
# 2. $p=1, a_0=1,a_1=1 \rightarrow P_1=xa_0$
# 3. $p=2, a_0=1,a_1=0 \rightarrow P_2=(1-3x^2)a_0$
# 4. $p=3, a_0=1,a_1=1 \rightarrow P_3=(x-\frac{5}{3}x^3)a_0$
# 5. $p=4, a_0=1,a_1=0 \rightarrow P_4=(\frac{35}{3}x^4-10x^2+1)a_0$
# 6. $p=5, a_0=1,a_1=1 \rightarrow P_5=(\frac{21}{5}x^5-\frac{14}{3}x^3+x)a_0$
# 7. $p=6, a_0=1,a_1=0 \rightarrow P_6=(\frac{-231}{5}x^6+63x^4-21x^2+1)a_0$
# 
# These are the Legendre Polynomials before the initial value requirement has been applied.
# We will also for $a_0$ for each problem based on $y(1)=1$.
# 
# Some are very easy to see, for $p=0,1$ $a_0=1$

# In[270]:


a=var('a')
P2=a*(1-3*x^2)
P3=a*(x-5/3*x^3)
P4=a*(35/3*x^4-10*x^2+1) 
P5=a*(21/5*x^5-14/3*x^3+x)
P6=a*(-231/5*x^6+63*x^4-21*x^2+1)
x=1
p2=solve([P2==1],a)
show(p2)
p3=solve([P3==1],a)
show(p3)
p4=solve([P4==1],a)
show(p4)
p5=solve([P5==1],a)
show(p5)
p6=solve([P6==1],a)
show(p6)


# Now the final equations are:
# 1. $p=0, a_0=1,a_1=0 \rightarrow P_0=1$
# 2. $p=1, a_0=1,a_1=1 \rightarrow P_1=x$
# 3. $p=2, a_0=1,a_1=0 \rightarrow P_2=\frac{1}{2}(3x^2-1)$
# 4. $p=3, a_0=1,a_1=1 \rightarrow P_3=\frac{1}{2}(5x^3-3x)$
# 5. $p=4, a_0=1,a_1=0 \rightarrow P_4=\frac{1}{8}(35x^4-30x^2+3)$
# 6. $p=5, a_0=1,a_1=1 \rightarrow P_5=\frac{1}{8}(63x^5-70x^3+15x)$
# 7. $p=6, a_0=1,a_1=0 \rightarrow P_6=\frac{1}{16}(231x^6-315x^4+105x^2-5)$
# 
# These are the Legendre Polynomials!!

# 3. Plot the first seven Legendre polynomials on a single plot and conjecture a property for these polynomials that depends on whether p is even or odd.

# In[271]:


x=var('x')
P0=1 #green
P1=x #red
P2=(-1/2)*(1-3*x^2) #yellow
P3=(-3/2)*(x-5/3*x^3) #purple
P4=(3/8)*(35/3*x^4-10*x^2+1) #orange
P5=(15/8)*(21/5*x^5-14/3*x^3+x) #aqua
P6=(-5/16)*(-231/5*x^6+63*x^4-21*x^2+1) #blue

show(sum([plot(P0,color='green')
          +plot(P1,color='red')
          +plot(P2,color='yellow')
          +plot(P3,color='purple')
          +plot(P4,color='orange')
          +plot(P5,color='aqua')
          +plot(P6)]))


# Based on the plot we can conjecture that even polynomials (not including zero) will have an odd number of minimum and maximum points added together and odd polynomials (not including one or zero) will have an even number of minumum and maximum points added together.

# 4. Consider Hermite’s equation (1).  For a positive integer p, how should you choose values in {0,1} for $a_0$ and $a_1$ so that the cooresponding power series solution $h_p(x)$is actually a polynomial? Compute $h_p(x)$for p ∈ {0,1,2,3,4,5,6}.

# Similar to question 2 above, some of these polynomials have already been calculated. And again, the be considered a polynomials solution, the series must terminate. So we we only be looking at the series that terminates from the solutions above and finding a pattern to calculating the polynomials. All of the solutions will be listed together at the end of this section.

# Hermite's Equation (1) $$y''-2xy'+2py=0$$
# 
# When:
# 1. $p=0$, $y''-2xy'=0$
# 2. $p=1$, $y''-2xy'+2(1)y=0$
# 3. $p=2$, $y''-2xy'+2(2)y=0$
# 4. $p=3$, $y''-2xy'+2(3)y=0$
# 5. $p=4$, $y''-2xy'+2(4)y=0$
# 6. $p=5$, $y''-2xy'+2(5)y=0$
# 7. $p=6$, $y''-2xy'+2(6)y=0$
# 
# From the work done in problem 1 we can see that 

# In[247]:


p=var('p')
solve([ydd0-2*x*yd1+2*p*y],an2)


# In[248]:


show(an2 == 2*(an*n - an*p)/(n^2 + 3*n + 2))


# This is $$a_{n+2}=a_n\frac{2(n-p)}{(n+1)(n+2)}$$

# In[250]:


#a0=1,a1=0,p=4
p=4
def coef(an,n):
    return(2*(an*n - an*p)/(n^2 + 3*n + 2))
a_nlist=[1,0] #a_0=1, a_1=0
for i in range (0,10):
    new=coef(a_nlist[i],i)
    a_nlist.append(new)
    
show(a_nlist)


# $$y_4=\frac{4}{3}x^4-4x^2+1$$

# In[252]:


#a0=0,a1=1,p=5
p=5
def coef(an,n):
    return(2*(an*n - an*p)/(n^2 + 3*n + 2))
a_nlist=[0,1] #a_0=0, a_1=1
for i in range (0,10):
    new=coef(a_nlist[i],i)
    a_nlist.append(new)
    
show(a_nlist)


# $$y_5=\frac{4}{15}x^5-\frac{4}{3}x^3+x$$

# In[253]:


#a0=1,a1=0,p=4
p=6
def coef(an,n):
    return(2*(an*n - an*p)/(n^2 + 3*n + 2))
a_nlist=[1,0] #a_0=1, a_1=0
for i in range (0,10):
    new=coef(a_nlist[i],i)
    a_nlist.append(new)
    
show(a_nlist)


# $$y_6=-\frac{8}{15}x^6+4x^4-6x^2+1$$

# 1. $p=0, a_0=1,a_1=0 \rightarrow h_0(x)=1$
# 2. $p=1, a_0=0,a_1=1 \rightarrow h_1(x)=x$
# 3. $p=2, a_0=1,a_1=0 \rightarrow h_2(x)=1-2x^2$
# 4. $p=3, a_0=0,a_1=1 \rightarrow h_3(x)=x-\frac{2}{3}x^3$
# 5. $p=4, a_0=1,a_1=0 \rightarrow h_4(x)=1-4x^2+\frac{4}{3}x^4$
# 6. $p=5, a_0=0,a_1=1 \rightarrow h_5(x)=\frac{4}{5}x^5-\frac{4}{3}x^3+x$
# 7. $p=6, a_0=1,a_1=0 \rightarrow h_6(x)=\frac{-8}{15}x^6+4x^4-6x^2+1$
# 
# There are the Hermite's polynomials!!

# 5. Let $G(x,t)=exp(2xt−t^2)$ and, for a positive integer $j$, define $$H_j(x)=\frac{\delta^jG}{\delta t^j}|_{t=0}$$ Then, using Taylor's Formula, we get $$\sum_{n=0}^\infty\frac{H_n(x)}{n!}t^n$$ which is a formal taylor series for G. Show that $$\frac{\delta^2G}{\delta{x^2}}-2x\frac{\delta{G}}{\delta{x}}+2(2xt-2t^2)G=0$$ Assume that G can be expressed as the power series formula and that it can be differentiated term wise with respect to x and t to argue that,  for any fixed non-negative integer, $H_n$ satisfies Hermite’s equation. Finally, notice that $H_n$ is a scalar multiple of $h_n$ for every such n. The polynomials $H_n$ are called Hermite polynomials.  Make a conjecture for a formula that gives the scalars relating $h_n$ and $H_n$.

# So here we see that $G(x,t)=e^{2xt-t^2}$
# 
# That first equation is holding x constant for making the derivative. So we want to find $H_j(x)$ by taking the derivative of the above equation while treating x as a constant.
# 
# The first thing to do is make sure the ODE above can be evaluated correctly.

# In[331]:


x,t,n=var('x,t,n')
G=e^(2*x*t-t^2)
dgx=G.diff(x)
dgx2=G.diff(x,2)

ode=dgx2-2*x*dgx+2*(2*x*t-2*t^2)*G
ode.full_simplify()


# This proves the equation from the question above.

# In[332]:


Glist=[]
for i in range (0,5):
    new=G.diff(t,i)
    Glist.append(new)
    print('dG',i,'=',new)
    
show(Glist)


# In[333]:


t=0
G=e^(-t^2 + 2*t*x)
dG=-2*(t - x)*e^(-t^2 + 2*t*x)
d2G=4*(t - x)^2*e^(-t^2 + 2*t*x) - 2*e^(-t^2 + 2*t*x)
d3G=-8*(t - x)^3*e^(-t^2 + 2*t*x) + 12*(t - x)*e^(-t^2 + 2*t*x)
d4G=16*(t - x)^4*e^(-t^2 + 2*t*x) - 48*(t - x)^2*e^(-t^2 + 2*t*x) + 12*e^(-t^2 + 2*t*x)
show('H0=G=',G)
show('H1=dG=',dG)
show('H2=d2G=',d2G)
show('H3=d3G=',d3G)
show('H4=d4G=',d4G)
H0=G
H1=dG
H2=d2G
H3=d3G
H4=d4G


# So from question 4 and this question we want to see a relationship between $h_n$ and $H_n$.
# 
# $$h_0(x)=1$$
# $$h_1(x)=x$$
# $$h_2(x)=1-2x^2$$
# $$h_3(x)=x-\frac{2}{3}x^3$$
# $$h_4(x)=1-4x^2+\frac{4}{3}x^4$$

# Here you can see that these equations are related based on a scalar.$$H_0=h_0$$ $$H_1=2h_1$$ $$H_2=-2h_2$$ $$H_3=-12h_3$$ $$H_4=12h_4$$ 

# 6. Plot $h_j$ for $j=0,...,4$ on a single plot. Do the same for $H_j$ for $j=0,...,4$.
# 
# The equations for $h_j$ are pulled from problem number 4.

# In[264]:


h0=1
h1=x
h2=1-x^2
h3=x-2/3*x^3
h4=1-4*x^2+4/3*x^4

show(sum([plot(h0,color='green')
          +plot(h1,color='red')
          +plot(h2,color='yellow')
          +plot(h3,color='purple')
          +plot(h4,color='orange')]))


# In[334]:


show(sum([plot(H0,color='green')
          +plot(H1,color='red')
          +plot(H2,color='yellow')
          +plot(H3,color='purple')
          +plot(H4,color='orange')]))


# 7. Suppose that $\alpha$, $\beta$, and $\gamma$ are fixed real numbers. The hypergeometic equation is $$x(1-x)y''+(\gamma-(\alpha+\beta+1)x)y'-\alpha\beta y=0$$ This equation has regular singular points at 0,1, and $\infty$.  (It turns out that any second order linear difierential equation with three regular singular points is equivalent to the hypergeometric equation, so it is quite important.)
# 
# Suppose that $\gamma$ is not a non-positive integer. The hypergeometric series is the series $$F(\alpha,\beta,\gamma,x)=1+\frac{\alpha\beta}{1!\gamma}x+\frac{\alpha(\alpha+1)\beta(\beta+1)}{2!\gamma(\gamma+1)}x^2+\frac{\alpha(\alpha+1)(\alpha+2)\beta(\beta+1)()\beta+2}{3!\gamma(\gamma+1)(\gamma+2)}x^3+...$$ 
# 
# Show that $F(\alpha,\beta,\gamma,x)$ converges on $(-1,1)$ and is a solution for the ODE above corresponding to the indical root zero.

# So based on the function above, $$F(\alpha,\beta,\gamma,x)=\sum_0^\infty\frac{(\alpha)_n(\beta)_n}{n!(\gamma)_n}x^n$$
# 
# To do this question we are going to take the ratio test $$\lim{n\to\infty}\left\lvert\frac{a_{n+1}}{a_n}\right\rvert$$ To make this easier $\alpha=a$, $\beta=b$, and $\gamma=c$ 
# 
# Because the limit is multiplied by x we are going to take the reciprocal R. The interval of convergence will be (-R, R). Also, because most of the terms can cancel we will end up with $$\lim{n\to\infty}\left\lvert\frac{(c+n)(n+1)}{(a+n)(b+n)}\right\rvert$$ From this we can see that the limits has an n-squared term over another one and this shows that the it converges on (-1, 1).

# An indical root of zero is when m=0 from this equation: $$y=x^m\sum_{n=0}^\infty d_nx^n$$

# In[392]:


x,m,n=var('x,m,n')
dn,dn1,dn2=var('dn,dn1,dn2')
y=dn*x^(n+m)
yd1=(n+m)*dn*x^(n+m-1)
yd0=(n+m+1)*dn1*x^(n+m)
ydd2=(n+m)*(n+m-1)*dn*x^(n+m-2)
ydd1=(n+m)*(n+m+1)*dn1*x^(n+m-1)
ydd0=(n+m+2)*(n+m+1)*dn2*x^(n+m)


# Now we are going to use these and plug them into the original hypergeometric equation. $$x(1-x)y''+(\gamma-(\alpha+\beta+1)x)y'-\alpha\beta y=0$$ 
# 
# In order for the equation to be useable we want all of the exponents of x to be equal to x+m-1. And we are going to assume that $d_0$ is not equal to zero.

# In[445]:


a,b,c=var('a,b,c')
d0=var('d0')
m=var('m')
F=dn*(n+m)*(n+m-1)*x^(n+m-1)-dn1*(m+m-1)*(n+m-2)*x^(n+m-1)+c*dn*(n+m)*x^(n+m-1)-(1+a+b)*dn1*(n+m-1)*x^(n+m-1)-a*b*dn1*x^(n+m-1)
#Now this equation is not perfect and we want to be able to take the first term of the sums
#(they can't be shown in this previous equation) and show that we have a recursion relation
rr=solve([(m*(m-1)+c*m)],m)
show(rr)


# Because we took out the first term, now the equation becomes:

# In[446]:


dn=var('dn')
F=((n+m)*(n+m-1)+c*(n+m))*dn+(-1*(n+m-1)*(n+m-1)-(1+a+b)*(n+m-1)-a*b)*dn1
print(F)


# In[447]:


rr2=solve([(c*(m + n) + (m + n)*(m + n - 1))*dn - (a*b + (a + b + 1)*(m + n - 1) + (m + n - 1)^2)*dn1],dn)
print(rr2)
show(rr2)


# In[434]:


m=0
n=var('n')
dn = ((a + b - 1)*dn1*m + dn1*m^2 + dn1*n^2 + ((a - 1)*b - a)*dn1 + ((a + b - 1)*dn1 + 2*dn1*m)*n)/((c - 1)*m + m^2 + (c + 2*m - 1)*n + n^2)
dn.full_simplify().show()


# In[443]:


n=1
d0=var('d0')
d1=((a + b)*d0*n + ((a - 1)*b - a + 1)*d0)/((c - 1)*n + n^2)
print('d1=',d1)
d1.full_simplify().show()


# So when we review the solution we found in the first part of this question we can see that when n=1 we get the same solution as we did to find the radius of convergence.

# 8. Suppose that $\gamma$ is not an integer and express a solution to the ODE in problem 7 corresponding to the indical root $1-\gamma$ in terms of F. Give the general solution to the ODE in terms of F. 
# 
# Knowing that the indical root (m) is equal to $1-\gamma$ we can use the following:

# In[450]:


m=1-c
n=var('n')
dn = ((a + b - 1)*dn1*m + dn1*m^2 + dn1*n^2 + ((a - 1)*b - a)*dn1 + ((a + b - 1)*dn1 + 2*dn1*m)*n)/((c - 1)*m + m^2 + (c + 2*m - 1)*n + n^2)
dn.full_simplify()


# In[451]:


n=1
d0=var('d0')
d1=-((a + b - 2*c + 1)*dn1*n + dn1*n^2 + (a*b - (a + b + 1)*c + c^2)*dn1)/((c - 1)*n - n^2)
print('d1=',d1)
d1.full_simplify().show()


# This gives us a recurssion relation and we can make sure that this will evaluate into the the original equation (kind of) because the m is equal to $1-\gamma$ our answers will be slightly different. But it will still evaluate to $$F(\alpha-\gamma+1,\beta-\gamma+1,2-\gamma,x)$$

# 9. Use the substitution $\frac{1}{2}(1-x)=t$ to show that a particular solution to Legendre's equation (2) is given by $$F\left(p+1,-p,1,\frac{1}{2}(1-x)\right)$$. So based on this $$a=p+1$$ $$b=-p$$$$c=1$$$$t=\frac{1}{2}(1-x)$$
# 
# Now putting this into the original equation will give us a useful equation. This will be shown in the code below.

# In[388]:


p,x=var('p,x')
a=p+1
b=-p
c=1
t=1/2*(1-x)
f=(1+a*b/c*t
   +a*(a+1)*b*(b+1)/(factorial(2)*c*(c+1))*t^2
   +a*(a+1)*(a+2)*b*(b+1)*(b+2)/(factorial(3)*c*(c+1)*(c+2))*t^3
   +a*(a+1)*(a+2)*(a+3)*b*(b+1)*(b+2)*(b+3)/(factorial(4)*c*(c+1)*(c+2)*(c+3))*t^4
   +a*(a+1)*(a+2)*(a+3)*(a+4)*b*(b+1)*(b+2)*(b+3)*(b+4)/(factorial(5)*c*(c+1)*(c+2)*(c+3)*(c+4))*t^5
   +a*(a+1)*(a+2)*(a+3)*(a+4)*(a+5)*b*(b+1)*(b+2)*(b+3)*(b+4)*(b+5)/(factorial(6)*c*(c+1)*(c+2)*(c+3)*(c+4)*(c+5))*t^6
   +a*(a+1)*(a+2)*(a+3)*(a+4)*(a+5)*(a+6)*b*(b+1)*(b+2)*(b+3)*(b+4)*(b+5)*(b+6)/(factorial(7)*c*(c+1)*(c+2)*(c+3)*(c+4)*(c+5)*(c+6))*t^7)
f.full_simplify()


# Now we are going to plug in different values of p to show that this equation is a particular solution to equation 2. $$(1-x^2)y''-2xy'+p(p+1)y=0$$

# In[389]:


df=f.diff(x)
ddf=f.diff(x,2)
ode=(1-x^2)*ddf-2*x*df+p*(p+1)*f
ode.full_simplify().show()


# In[390]:


Flist=[]
for i in range (0,7):
    p=i
    new=(-1/3251404800*p^14 - 1/464486400*p^13 + 1/17203200*p^12 + 1/2654208*p^11 - 2207/464486400*p^10 - 
         67/2457600*p^9 + 709207/3251404800*p^8 + 1/3251404800*(p^14 + 7*p^13 - 91*p^12 - 637*p^11 + 3003*p^10 + 21021*p^9
        - 44473*p^8 - 311311*p^7 + 296296*p^6 + 2074072*p^5 - 773136*p^4 - 5411952*p^3 + 518400*p^2 + 3628800*p)*x^7 + 
         96641/92897280*p^7 - 1/464486400*(p^14 + 7*p^13 - 105*p^12 - 721*p^11 + 3773*p^10 + 25641*p^9 - 58795*p^8 - 
        397243*p^7 + 403326*p^6 + 2716252*p^5 - 1068200*p^4 - 7182336*p^3 + 720000*p^2 + 4838400*p)*x^6 - 
         457157/77414400*p^6 + 1/154828800*(p^14 + 7*p^13 - 119*p^12 - 805*p^11 + 4879*p^10 + 31941*p^9 - 83197*p^8 - 
        533575*p^7 + 602084*p^6 + 3817072*p^5 - 1638784*p^4 - 10330320*p^3 + 1115136*p^2 + 7015680*p)*x^5 - 
         356257/16588800*p^5 - 1/92897280*(p^14 + 7*p^13 - 133*p^12 - 889*p^11 + 6321*p^10 + 39921*p^9 - 127759*p^8 - 
        760627*p^7 + 1033690*p^6 + 5941012*p^5 - 2978808*p^4 - 16831584*p^3 + 2066688*p^2 + 11612160*p)*x^4 + 
         5099509/58060800*p^4 + 1/92897280*(p^14 + 7*p^13 - 147*p^12 - 973*p^11 + 8099*p^10 + 49581*p^9 - 202561*p^8 - 
        1118719*p^7 + 2161824*p^6 + 10620232*p^5 - 7194992*p^4 - 33500208*p^3 + 5227776*p^2 + 23950080*p)*x^3 + 
         29407/138240*p^3 - 1/154828800*(p^14 + 7*p^13 - 161*p^12 - 1057*p^11 + 10213*p^10 + 60921*p^9 - 317683*p^8 - 
        1648171*p^7 + 4772726*p^6 + 20354572*p^5 - 27683656*p^4 - 91342272*p^3 + 23218560*p^2 + 72576000*p)*x^2 - 
         136919/235200*p^2 + 1/464486400*(p^14 + 7*p^13 - 175*p^12 - 1141*p^11 + 12663*p^10 + 73941*p^9 - 483205*p^8 - 
            2389303*p^7 + 9975196*p^6 + 38611552*p^5 - 98807520*p^4 - 264909456*p^3 + 321546240*p^2 + 460857600*p)*x - 
         1163/1680*p + 1)
    Flist.append(new)
    print('f(',i,')=',new)
    
show(Flist)


# These equations above are very similar to the Legendre Polynomials that we found earlier. The only issue that we can see is that the scalar in front of the equations from problem 2. Below I will show this! The actual Legendre polynomials are: $$P_0=1$$
# $$P_1=x$$ $$P_2=\frac{1}{2}(3x^2-1)$$ $$P_3=\frac{1}{2}(5x^3-3x)$$ $$P_4=\frac{1}{8}(35x^4-30x^2+3)$$ $$P_5=\frac{1}{8}(63x^5-70x^3+15x)$$ $$P_6=\frac{1}{16}(231x^6-315x^4+105x^2-5)$$

# Using the information collected in this problem and in problem 2 we see that:
# 
# $$P_0=1=f_0$$
# $$P_1=x=f_1$$ 
# $$P_2=\frac{1}{2}(3x^2-1)=\frac{3}{2}x^2-\frac{1}{2}=f_2$$
# $$P_3=\frac{1}{2}(5x^3-3x)=\frac{5}{2}x^3-\frac{3}{2}x=f_3$$
# $$P_4=\frac{1}{8}(35x^4-30x^2+3)=\frac{35}{8}*x^4-\frac{15}{4}x^2+\frac{3}{8}=f_4$$
# $$P_5=\frac{1}{8}(63x^5-70x^3+15x)=\frac{63}{8}x^5-\frac{35}{4}x^3+\frac{15}{8}x=f_5$$
# $$P_6=\frac{1}{16}(231x^6-315x^4+105x^2-5)=\frac{231}{16}x^6-\frac{315}{16}x^4+\frac{105}{16}x^2-\frac{5}{16}=f_6$$

# 10. Curiosity Element

# For my curiosity element I wanted to look at how calculus is used in Chemistry. As a chemistry major this is pretty important to me and I wanted to look at applications for this calculus that would be useful to me.
# 
# Calculus is used with a lot of the data analysis within chemistry but first I want to look at physical chemistry. 
# 
# Physical Chemistry uses linear algebra, multivariable, calculus, and limits. In particular when working with thermodynamics and kinetics (a class I will be taking next semester). Because most of the math that must be done for the semester is in partial derivatives, calculus knowledge is very important. 
# 
# Next, Analytical chemistry uses calculus to evaluate and predict pH concentrations. As well as using the Nernst-Planck equation. 
# 
# I'm going to look at the Nernst-Planck equation here and you'll see that derivatives are used a lot. 

# This equation is $$\frac{\partial c}{\partial t}= -\nabla *J$$
# 
# With $$J=-\left[D\nabla c-uc+\frac{Dze}{k_BT}c(\nabla\phi+\frac{\partial A}{\partial t})\right]$$
# With this being the change in the concentration with respect to time with J being the diffusion flux density. This equation describes the flux of ions under the influcene of the concentration gradient and the electric field.
# 
# This equation is the derivative to Fick's Law of diffusion and so this is a great example of how derivates are used in chemistry.

# Differential equations are also used when describing and predicting how chemical reactions proceed. There are many first order differential equations and variations of them. Physical chemistry uses second order differential equations as well. 
# 
# One thing that I have used in my chemistry laboratory work is regression analysis. Regression analysis is based in calculus and is something that I plan on using going forward. Linear regression uses first order differential equations and power series to solve machine learning problems. 
# 
# These are all ways that I can use calculus in chemistry.

# A bad calculus joke for you!!
# 
# A guy gets on a bus and starts threatening everybody: "I'll integrate you! I'll differentiate you!" So everybody gets scared and runs away. Only one person stays. The guy comes up to him and says: "Aren't you scared, I'll integrate you, I'll differentiate you!" And the other guy says: "No, I am not scared, I am $e^x$."
# 

# In[ ]:





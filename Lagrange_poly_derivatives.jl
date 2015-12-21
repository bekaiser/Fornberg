# Derivative tests for Fornberg (1998) Lagrange polynomial derivatives
# Bryan Kaiser
# 12/20/15

using DataArrays
using PyPlot
using PyCall
@pyimport numpy as np
@pyimport pylab as py


## ============================================================================
# function declaration

function weights(z,x,n,m)
# From Bengt Fornbergs (1998) SIAM Review paper.
#  	Input Parameters
#	z location where approximations are to be accurate,
#	x(0:nd) grid point locations, found in x(0:n)
#	n one less than total number of grid points; n must
#	not exceed the parameter nd below,
#	nd dimension of x- and c-arrays in calling program
#	x(0:nd) and c(0:nd,0:m), respectively,
#	m highest derivative for which weights are sought,
#	Output Parameter
#	c(0:nd,0:m) weights at grid locations x(0:n) for derivatives
#	of order 0:m, found in c(0:n,0:m)
#      	dimension x(0:nd),c(0:nd,0:m)

	c = zeros(n+1,m+1);
	c1 = 1.0;
	c4 = x[1+0]-z;
	for k=0:m
  		for j=0:n
    			c[1+j,1+k] = 0.0;
  		end
	end
	c[1+0,1+0] = 1.0;
	for  i=1:n
  		mn = min(i,m);
  		c2 = 1.0;
  		c5 = c4;
  		c4 = x[1+i]-z;
  		for j=0:i-1
    			c3 = x[1+i]-x[1+j];
    			c2 = c2*c3;
    			if (j == i-1)
      				for k=mn:-1:1
        				c[1+i,1+k] = c1*(k*c[1+i-1,1+k-1]-c5*c[1+i-1,1+k])/c2;
			      	end
      			c[1+i,1+0] = -c1*c5*c[1+i-1,1+0]/c2;
    			end
    			for k=mn:-1:1
      				c[1+j,1+k] = (c4*c[1+j,1+k]-k*c[1+j,1+k-1])/c3;
    			end
    			c[1+j,1+0] = c4*c[1+j,1+0]/c3;
  		end
  		c1 = c2;
	end
   	return c
end # weights function


## ============================================================================
# domain parameters

const Lx = 3000.0 # km, domain size
const Ly = Lx # km
const Lxcenter = 0.0 # x value @ the center of the grid
const Lycenter = 0.0 # y value @ the center of the grid
const N = 2^20 # series length (must be at least even)
const dx = Lx/float(N) # km, uniform longitudinal grid spacing
const dy = Ly/float(N) # km, uniform longitudinal grid spacing
x = collect(0.5*dx:dx:dx*N)-(Lx/2.0-Lxcenter) # km, centered uniform grid 


## ============================================================================
# signal

k = 2.0*pi/(Lx/5.0)
s = sin(x.*k)
dsdx = cos(x.*k).*k
d2sdx2 = -sin(x.*k).*k^2

figure(1)
plot(x,s,"b",x,dsdx,"r",x,d2sdx2,"g")
xlabel("x")
title("signal")
show()


##=============================================================================
# compute the uniform grid weights

c = weights(x[3],x[1:5],4,2)'; # 5 point stencil weights
dsF = zeros(1,N); # derivative of the input function
m = 1 # order of the derivative
   
dsF = zeros(Float64,N,1)
for j = 3:N-2      
	jm2 = j-2;
        jp2 = j+2;
        dsF[j,1] = vecdot(c[m+1,:],s[jm2:jp2]);
end
error = abs(dsdx-dsF)

@show(N)
figure(2)
semilogy(x[3:N-2],error[3:N-2],"k")
xlabel("x")
title("Lagrange polynomial 1st derivative error")
show()

const c = 299792458.
const mu0 = 4*pi*1e-7
const eps0 = 1/(mu0*c^2)

function Hertz_dipole (r, p, R, phi, f, t=0, epsr=1.)
    #=
    Calculate E field and B field strength of hertzian dipole(s)
    p: array of dipole moment
    R: array of dipole position
    r: observation point
    f: frequency
    t: time
    phi: dipole phase angles (0..2pi)
    return: field strength at observation point r at time t (3-tuple: Ex, Ey, Ez)
    =#
    rprime=r-R  # r'=r-R
    magrprime=sqrt(sum(rprime.^2)) # |r-R|
    w=2*pi*f  # omega
    k=w/c     # wave number
    krp=k*magrprime  # k*|r-R|
    rprime_cross_p = cross(rprime, p) # (r-R) x p
    rp_c_p_c_rp = cross(rprime_cross_p, rprime) # ((r-R) x p) x (r-R)
    expfac=exp(1im*(-w*t+krp-phi))/(4*pi*eps0*epsr)
    E=expfac*(w^2/(c^2*magrprime^3)*rp_c_p_c_rp+(1/magrprime^3-w*im/(c*magrprime^2))*(3*rprime*dot(rprime,p)/magrprime^2-p))
    B=expfac/(magrprime*c^3)*(w^2*rprime_cross_p)/magrprime*(1-c/(im*w*magrprime))
    return E,B
  end


function Hertz_dipole_ff (r, p, R, phi, f, t=0, epsr=1.)
  #=
  Calculate E field and B field strength in far field of hertzian dipole(s)
  p: array of dipole moment
  R: array of dipole position
  r: observation point
  f: frequency
  t: time
  phi: dipole phase angles (0..2pi)
  return: field strength at observation point r at time t (3-tuple: Ex, Ey, Ez)
  =#
  rprime=r-R  # r'=r-R
  magrprime=sqrt(sum(rprime.^2)) # |r-R|
  w=2*pi*f  # omega
  k=w/c     # wave number
  krp=k*magrprime  # k*|r-R|
  rprime_cross_p = cross(rprime, p) # (r-R) x p
  rp_c_p_c_rp = cross(rprime_cross_p, rprime) # ((r-R) x p) x (r-R)
  expfac=exp(1im*(-w*t+krp-phi))/(4*pi*eps0*epsr)
  E=expfac*(w^2/(c^2*magrprime^3)*rp_c_p_c_rp)
  B=expfac/(magrprime^2*c^3)*(w^2*rprime_cross_p)
  return E,B
end


function Hertz_dipole_nf (r, p, R, phi, f, t=0, epsr=1.)
    #=
  Calculate E field and B field strength in near field of hertzian dipole(s)
  p: array of dipole moment
  R: array of dipole position
  r: observation point
  f: frequency
  t: time
  phi: dipole phase angle (0..2pi)
  return: field strength at observation point r at time t (3-tuple: Ex, Ey, Ez)
  =#
  rprime=r-R  # r'=r-R
  magrprime=sqrt(sum(rprime.^2)) # |r-R|
  w=2*pi*f  # omega
  k=w/c     # wave number
  krp=k*magrprime  # k*|r-R|
  rprime_cross_p = cross(rprime, p) # (r-R) x p
  expfac=exp(1im*(-w*t+krp-phi))/(4*pi*eps0*epsr)
  E=expfac*((1/magrprime^3-w*im/(c*magrprime^2))*(3*rprime*dot(rprime,p)/magrprime^2-p))
  B=expfac/(magrprime*c^2)*(w^2*rprime_cross_p)/magrprime*(-1/(im*w*magrprime))
  return E,B
end

#Main Program

#observation points
const  nx=401
const  xmax=2
const  nz=201
const  zmax=1
const  x=linspace(-xmax,xmax,nx)
const  y=0
const  z=linspace(-zmax,zmax,nz)

#dipole
const  freq=1000e6
#dipole moment
#total time averaged radiated power P= 1 W dipole moment => |p|=sqrt(12πcP/µOω⁴)
const Pow = 1
const norm_p = sqrt(12*pi*c*Pow/(mu0*(2*pi*freq)^4))

const  p = [0,0,norm_p]
const  R = [0,0,0]
#dipole phase
const  phases_dip = 0

const  t0 = 1/freq/10
const  t1 = 5/freq
const  nt = int(t1/t0)
const  t = linspace(t0,t1,nt)

println("Computing the radiation...")
using PyPlot
pygui(false)


figure(num=1,figsize=(10,5), dpi=300)
for k=1:nt
  P=zeros(nx,nz)
  for i=1:nx
    for j=1:nz
      r = [x[i],y,z[j]]
      E,B = Hertz_dipole (r, p, R, phases_dip, freq, t[k])
      S = real(E).^2#0.5*numpy.cross(E.T,conjugate(B.T))
      P[i,j] = sum(S)
    end
  end
  percent = k/nt*100
  print("$percent/100")  #Radiation diagram
  pcolor(x,z,P[:,:]',cmap="hot")
  fname = "img_$k.png"
  clim(0,1000)
  axis("scaled")
  xlim(-xmax,xmax)
  ylim(-zmax,zmax)
  xlabel("x/m")
  ylabel("z/m")
  timestep=round(t[k]*1e10)/10
  title("t=$timestep ns")
  println ("Saving frame $fname")
  savefig("img_$k.png",bbox="tight")
  clf()
end

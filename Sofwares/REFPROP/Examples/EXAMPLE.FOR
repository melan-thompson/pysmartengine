c     Example program using routines in REFPROP
c       Units are K, kPa, mol/dm^3, mole fraction, J/mol, J/mol-K, m/s
c       uPa-s, W/m-K, and N/m

C   *** See the file MANUAL.TXT in the FORTRAN directory for more routines ***

      program EXAMPLE
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      parameter (ncmax=20)   !max number of components in mixture
      dimension x(ncmax),xliq(ncmax),xvap(ncmax),f(ncmax)
      character hrf*3, herr*255
      character*255 hf(ncmax),hfmix
      dimension xdew(ncmax),ybub(ncmax)  !Only used for TQFL2 below

c     character*3 htype,hmix,hcomp(1:ncmax)

c...Print out error messages when iprnterr=1
c...Print out warning and error messages when iprnterr=-1
c...Pause when equal to 3, -3, and print out error messages similar to 1, -1
c     common /prnterr/ iprnterr
c     iprnterr=-3

c...Use call to SETMOD to change the equation of state for any of the
c.....pure components from the default (recommended) values.
c     i=3
c     htype='EOS'
c     hmix='HMX'
c     hcomp(1)='FES'       !Nitrogen
c     hcomp(2)='FES'       !Oxygen
c     hcomp(3)='FES'       !Argon
c     call SETMOD (i,htype,hmix,hcomp,ierr,herr)

c...If the fluid files are located in a directory that the code cannot
c.....find, make a call to SETPATH
c     call SETPATH('C:\Program Files\Refprop')

c...Call SETUP to initialize the program and set the pure fluid component name
      nc=1                      !Number of components
      hf(1)='nitrogen.fld'      !Fluid name
      hfmix='hmx.bnc'           !Mixture file name
      hrf='DEF'                 !Reference state (DEF means default)
      call SETUP (nc,hf,hfmix,hrf,ierr,herr)
      if (ierr.ne.0) write (*,*) herr

c...Example code for calling a pseudo-pure fluid is given in EX-PPF.FOR

c...For a mixture, use the following setup instead of the lines above.
c     nc=3
c     hf(1)='nitrogen.fld'
c     hf(2)='argon.fld'
c     hf(3)='oxygen.fld'
c     hfmix='hmx.bnc'
c     hrf='DEF'
c     call SETUP (nc,hf,hfmix,hrf,ierr,herr)
c     if (ierr.ne.0) write (*,*) herr
c     x(1)=.7812d0     !Air composition
c     x(2)=.0092d0
c     x(3)=.2096d0

c...Call SETREF to set the reference state.  Valid inputs are:
c           'NBP':  h,s=0 at normal boiling point(s)
c           'ASH':  h,s=0 for sat liquid at -40 C (ASHRAE convention)
c           'IIR':  h=200 kJ/kg, s=1.0 kJ/kg-K for sat liq at 0 C (IIR convention)
c           'DEF':  default reference state as specified in fluid file
c           'OTH':  other, as specified by h0, s0, t0, p0
c...Only need to call SETREF if 'DEF' state is not used.
c...This call sets h=2000 J/mol and s=100 J/mol-K at 298.15 K and 101.325 MPa
c     call SETREF ('OTH',1,x,2000.d0,100.d0,298.15d0,101.325d0,ierr,herr)
c...If the reference state for a mixture is to be set to 'Currently defined
c...mixture' (as done in the graphical interface), then call SETREF with this
c...(be sure that the composition array x is initialized):
c     call SETREF (hrf,2,x,0.d0,0.d0,0.d0,0.d0,ierr,herr)


c...Get molecular weight (wm), triple point temperature, normal boiling point
c.....temperature, critical temperature, pressure, density, and compressibility
c.....factor, dipole moment, and universal gas constant for component #1
c.....(First input specifies component number)
      call INFO (1,wm,ttp,tnbp,tc,pc,dc,zc,acf,dip,rgas)
      write (*,1000) 'WM,ACF,DIP,TTP,TNBP   ',wm,acf,dip,ttp,tnbp
      write (*,1000) 'TC,PC,DC,RGAS         ',tc,pc,dc,rgas

c...Get approximate critical point for the mixture, only works for type 1
c.....mixtures.
      call CRITP (x,tcrit,pcrit,Dcrit,ierr,herr)
      write (*,1000) 'TCmix,PCmix,Dcmix     ',tcrit,pcrit,Dcrit,x(1)

c...Calculate molecular weight of a mixture
c     wm=WMOL(x)

c...Get saturation properties given t,x; for i=1: x is liquid phase
c.....                                   for i=2: x is vapor phase
      i=1
      t=100
      call SATT (t,x,i,p,dl,dv,xliq,xvap,ierr,herr)
      write (*,1000) 'P,Dl,Dv,xl(1),xv(1)   ',p,dl,dv,xliq(1),xvap(1)
      i=2
      call SATT (t,x,i,p,dl,dv,xliq,xvap,ierr,herr)
      write (*,1000) 'P,Dl,Dv,xl(1),xv(1)   ',p,dl,dv,xliq(1),xvap(1)

c...Calculate saturation properties at a given p. i is same as SATT
      i=2
      call SATP (p,x,i,t,dl,dv,xliq,xvap,ierr,herr)
      write (*,1000) 'T,Dl,Dv,xl(1),xv(1)   ',t,dl,dv,xliq(1),xvap(1)

c...Other saturation routines are given in SAT_SUB.FOR

      t=300
      p=20000

c...Calculate d from t,p,x
c...If phase is known: (j=1: Liquid, j=2: Vapor)
      j=1
      call TPRHO (t,p,x,j,0,d,ierr,herr)
      write (*,1000) 'T,P,D                 ',t,p,d

c...If phase is not known, call TPFLSH
c...Calls to TPFLSH are much slower than TPRHO since SATT must be called first.
c.....(If two phase, quality is returned as q)
      call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
      write (*,1000) 'T,P,D,H,CP            ',t,p,d,h,cp

c...Calculate pressure (p), internal energy (e), enthalpy (h), entropy (s),
c.....isochoric (cv) and isobaric (cp) heat capacities, speed of sound (w),
c.....and Joule-Thomson coefficient (hjt) from t,d,x
c.....(subroutines THERM2 and THERM3 contain more properties, see PROP_SUB.FOR)
      call THERM (t,d,x,p,e,h,s,cv,cp,w,hjt)

c...Calculate pressure
      call PRESS (t,d,x,p)

c...Calculate fugacity
      call FGCTY (t,d,x,f)

c...Calculate second and third virial coefficients
      call VIRB (t,x,b)
      call VIRC (t,x,c)
      write (*,1000) 'F,B,C                 ',f(1),b,c

c...Calculate the derivatives: dP/dD, d^2P/dD^2, dP/dT  (D indicates density)
C...(dD/dP, dD/dT, and dB/dT are also available, see PROP_SUB.FOR)
      call DPDD (t,d,x,dpdrho)
      call DPDD2 (t,d,x,d2pdd2)
      call DPDT (t,d,x,dpt)
      write (*,1000) 'dP/dD,d2P/dD2,dP/dT   ',dpdrho,d2pdd2,dpt


c...Calculate derivatives of enthalpy with respect to T, P, and D
      call DHD1(t,d,x,dhdt_d,dhdt_p,dhdd_t,dhdd_p,dhdp_t,dhdp_d)
      write (*,1000) 'Enthalpy derivatives  ',
     &       dhdt_d,dhdt_p,dhdd_t,dhdd_p/1000.d0,dhdp_t

c...Calculate surface tension
      tt=100
      call SURFT (tt,dl,x,sigma,ierr,herr)
      write (*,1000) 'T,SURF. TN.           ',tt,sigma

c...Calculate viscosity (eta) and thermal conductivity (tcx)
      call TRNPRP (t,d,x,eta,tcx,ierr,herr)
      write (*,1000) 'VIS.,TH.CND.          ',eta,tcx*1000.D0

c...General property calculation with inputs of t,d,x
      call TDFLSH (t,d,x,pp,dl,dv,xliq,xvap,q,e,h1,s,cv,cp,w,ierr,herr)
      write (*,1000) 'T, D, P from TDFLSH   ',t,d,pp/1000.d0

c...General property calculation with inputs of p,d,x
      call PDFLSH (p,d,x,tt,dl,dv,xliq,xvap,q,e,h1,s,cv,cp,w,ierr,herr)
      write (*,1000) 'T, D, P from PDFLSH   ',tt,d,p/1000.d0

c...General property calculation with inputs of p,h,x
      call PHFLSH (p,h,x,tt,dd,dl,dv,xliq,xvap,q,e,s,cv,cp,w,ierr,herr)
      write (*,1000) 'T, D, P from PHFLSH   ',tt,dd,p/1000.d0

c...General property calculation with inputs of p,s,x
      call PSFLSH (p,s,x,tt,dd,dl,dv,xliq,xvap,q,e,h1,cv,cp,w,ierr,herr)
      write (*,1000) 'T, D, P from PSFLSH   ',tt,dd,p/1000.d0

c...General property calculation with inputs of d,h,x
      call DHFLSH (d,h,x,tt,pp,dl,dv,xliq,xvap,q,e,s,cv,cp,w,ierr,herr)
      write (*,1000) 'T, D, P from DHFLSH   ',tt,d,pp/1000.d0

c...General property calculation with inputs of t,h,x
c     kr--flag specifying desired root for multi-valued inputs:
c         1=return lower density root
c         2=return higher density root
      kr=1
      call THFLSH (t,h,x,
     &             kr,pp,dd,dl,dv,xliq,xvap,q,e,s,cv,cp,w,ierr,herr)
      write (*,1000) 'T, D, P from THFLSH   ',t,dd,pp/1000.d0

c...Other general property calculation routines are given in FLSH_SUB.FOR
c...and FLASH2.FOR

c...Calculate melting pressure
      t=100.d0
      call MELTT (t,x,p,ierr,herr)
      write (*,1000) 'Melting pressure(MPa) ',p/1000.d0,t

c...Calculate melting temperature
c     call MELTP (p,x,tt,ierr,herr)
c     write (*,1000) 'Melting temperature(K)',tt,p/1000.d0

c...Calculate sublimation pressure
      t=200
      call SUBLT (t,x,p,ierr,herr)
      write (*,1000) 'Sublimation pr.(kPa)  ',p,t

c...Calculate sublimation temperature
c     call SUBLP (p,x,tt,ierr,herr)
c     write (*,1000) 'Sublimation temp.(K)  ',tt,p

c...Get limits of the equations and check if t,d,p is a valid point
c...Equation of state
c     call LIMITK ('EOS',1,t,d,p,tmin,tmax,Dmax,pmax,ierr,herr)
c...Viscosity equation
c     call LIMITK ('ETA',1,t,d,p,tmin,tmax,Dmax,pmax,ierr,herr)
c...Thermal conductivity equation
c     call LIMITK ('TCX',1,t,d,p,tmin,tmax,Dmax,pmax,ierr,herr)

c...Calculate pressure given a temperature and quality.  In this case,
c...q is set to 0.2, or 20% of the substance is vapor and 80% is liquid,
c...on a molar basis (to use a mass based q, set kq=2).
c...The substance is changed here to a mixture to illustrate the effect.
c...Note that xdew and ybub must be dimensioned to 20, but are not used (as
c...well as pbub, pdew, Dlbub, and Dvdew) if ksat=0.  See manual.txt for more.
      nc=2
      hf(1)='methane'
      hf(2)='ethane'
      hfmix='hmx.bnc'
      hrf='DEF'
      call SETUP (nc,hf,hfmix,hrf,ierr,herr)
      if (ierr.ne.0) write (*,*) herr
      x(1)=.7d0
      x(2)=.3d0
      t=100
      q=0.2
      kq=1
      ksat=0
      call TQFL2 (t,q,x,kq,ksat,pbub,pdew,Dlbub,Dvdew,ybub,xdew,
     &            p,Dl,Dv,xliq,xvap,ierr,herr)
      if (ierr.ne.0) write (*,*) herr
      write (*,1000) 'Sat. pressure (MPa)   ',p/1000.d0,t

C...Other routines are given in UTILITY.FOR

      stop
 1000 format (1x,a22,5f11.4)
      end

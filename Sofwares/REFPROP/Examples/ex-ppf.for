c     Example program using routines in REFPROP to call a pseudo-pure fluid
c       Units are K, kPa, mol/dm^3, mole fraction, J/mol, J/mol-K, m/s
c       uPa-s, W/m-K, N/m

      program EXAMPLE
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      parameter (ncmax=20)   !max number of components in mixture
      dimension x(ncmax),xliq(ncmax),xvap(ncmax),f(ncmax)
      character hrf*3, herr*255
      character*255 hf(ncmax),hfmix

      write (*,*)
      write (*,*) 'Example output for a pseudo-pure fluid'
c...Call SETUP to initialize the program for a pseudo-pure fluid
c.....A pseudo-pure fluid representes the properties of a mixture at a
c.....specific composition using a single equation of state similar to those
c.....used for other pure fluids.  All REFPROP calls are the same as for a
c.....regular pure fluid, except those to SATT and SATP, which are explained
c.....below.  The composition array x is not used and can be set to anything.
      i=1
      hf(1)='R410A.PPF'
      hfmix='hmx.bnc'
      hrf='DEF'
      call SETUP (i,hf,hfmix,hrf,ierr,herr)
      if (ierr.ne.0) write (*,*) herr

      t=300
      p=20000
      call TPFLSH (t,p,x,d,dl,dv,xliq,xvap,q,e,h,s,cv,cp,w,ierr,herr)
      write (*,1000) 'T,P,D,H,CP            ',t,p/1000.d0,d,h,cp

c...For saturation states, SATT or SATP must be called once for the liquid
c.....phase and once for the vapor phase since composition is not part of the
c.....equation of state
c...Call SATT to get liquid phase pressure and density at 300 K:
      t=300.d0
      i=1
      call SATT (t,x,i,p,dl,dv,xliq,xvap,ierr,herr)
      write (*,1000) 'P,Dl,Dv,xl(1),xv(1)   ',p,dl,dv,xliq(1),xvap(1)
c...Call SATT to get vapor phase pressure and density at 300 K:
      i=2
      call SATT (t,x,i,p,dl,dv,xliq,xvap,ierr,herr)
      write (*,1000) 'P,Dl,Dv,xl(1),xv(1)   ',p,dl,dv,xliq(1),xvap(1)

c...Call SATP to get liquid phase temperature and density at 1 MPa:
      p=1000.d0
      i=1
      call SATP (p,x,i,t,dl,dv,xliq,xvap,ierr,herr)
      write (*,1000) 'T,Dl,Dv,xl(1),xv(1)   ',t,dl,dv,xliq(1),xvap(1)
c...Call SATP to get vapor phase temperature and density at 1 MPa:
      i=2
      call SATP (p,x,i,t,dl,dv,xliq,xvap,ierr,herr)
      write (*,1000) 'T,Dl,Dv,xl(1),xv(1)   ',t,dl,dv,xliq(1),xvap(1)

      stop
 1000 format (1x,a22,5f11.4)
      end

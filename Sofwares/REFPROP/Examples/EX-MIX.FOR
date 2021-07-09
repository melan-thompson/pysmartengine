c  Example program showing how to setup a mixture using the information
c  contained in a .mix file (located in the ...\Refprop\Mixtures directory).

      program EXAMPLE
      implicit double precision (a-h,o-z)
      implicit integer (i-k,m,n)
      parameter (ncmax=20)        !max number of components in mixture
      dimension x(ncmax),xliq(ncmax),xvap(ncmax)
      character hrf*3
      character*255 hf(ncmax),hfmix,hflnme,herr

c...If the fluid files are located in a directory that the code cannot
c.....find, make a call to SETPATH
      call SETPATH('C:\Program Files\Refprop')

      hflnme='R410A.mix'
      hfmix='hmx.bnc'
      hrf='DEF'

c...The SETMIX subroutine returns the number of fluids [ncc], fluid
c...names [hf()], and the composition of the mixture [x()].
      call SETMIX (hflnme,hfmix,hrf,ncc,hf,x,ierr,herr)
      if (ierr.ne.0) write (*,*) herr
      if (ierr.gt.0) stop

      t=300

c...Mixture properties use the exact same routines as those described
c...in EXAMPLE.FOR

c...Calculate saturated liquid properties
      call SATT (t,x,1,p,dl,dv,xliq,xvap,ierr,herr)
      if (ierr.ne.0) write (*,*) herr
      write (*,1000) t,p,dl,dv,xliq(1),xvap(1)

c...Calculate saturated vapor properties
      call SATT (t,x,2,p,dl,dv,xliq,xvap,ierr,herr)
      if (ierr.ne.0) write (*,*) herr
      write (*,1000) t,p,dl,dv,xliq(1),xvap(1)

 1000 format (15f11.4)
      end

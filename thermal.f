
c.   .....   .c

      subroutine KawasakiSwap(acc)

      implicit none
      include 'parameter.inc'
      include 'system.inc'
      integer ix, iy
      integer ixp, iyp
      integer ip, ipp, imol
      integer ActSp
      double precision ran2, rn
      double precision dEn, De_bond
      integer s1, s2, s3, a1, a2, a3, k
      integer o1, o2, o3
      integer acc, acsp
      integer d_length, d_covered_length

      acc=0
      imol=0
      k=1
      d_length=0
      d_covered_length=0

      do while (imol .ne. 1) 

         ix=int(ran2(Iseed)*dble(Lx))
         iy=int(ran2(Iseed)*dble(Ly)) 
         ixp=ix
         iyp=iy

         rn=ran2(Iseed)

         if(rn .lt. 1.d0/3.d0) then
            ixp=Mod(ix+1,Lx)
         end if 

         if((rn .ge. 1.d0/3.d0).and.(rn .lt. 2.d0/3.d0)) then
            iyp=Mod(iy+1,Ly)
         end if 

         if(rn .ge. 2.d0/3.d0) then 
            if(Mod(iy,2) .eq. 0) ixp=Mod(ix-1+Lx,Lx)
            if(Mod(iy,2) .eq. 1) ixp=Mod(ix+1,Lx)
            iyp=Mod(iy+1,Ly)
         end if 

         imol=1

c new boundary conditions 
         if(fixed_boundary .eqv. .true.) then
            if((iy .eq. 0).or.(iy .eq. Ly-1)) imol=0
            if((iyp .eq. 0).or.(iyp .eq. Ly-1)) imol=0
         end if

      end do

      s1=sp(ix,iy)
      a1=anch(ix,iy)
      o1=occ(ix,iy)
      s2=sp(ixp,iyp)
      a2=anch(ixp,iyp)
      o2=occ(ixp,iyp)
      dEn=0.d0

      if(s1*s2 .lt. 0) then 

c. neighbours of ix iy (s1) 
         s3=sp(mod(ix-1+Lx,Lx),iy)
         a3=anch(mod(ix-1+Lx,Lx),iy)
         o3=occ(mod(ix-1+Lx,Lx),iy)
         dEn=dEn+De_bond(s1,a1,s3,a3,jnn,ja,jb,k)
         if(s1*s3 .lt. 0) then
            d_length=d_length-1
            if(o1.ne.-1) then
               if(o1.eq.o3) then
                  d_covered_length=d_covered_length-1
               end if
            end if
         else
            d_length=d_length+1
             if(o1.ne.-1) then
               if(o1.eq.o3) then
                  d_covered_length=d_covered_length+1
               end if
            end if
         end if
         
         s3=sp(mod(ix+1,Lx),iy)
         a3=anch(mod(ix+1,Lx),iy)
         o3=occ(mod(ix+1,Lx),iy)
         dEn=dEn+De_bond(s1,a1,s3,a3,jnn,ja,jb,k)
         if(s1*s3 .lt. 0) then
            d_length=d_length-1
            if(o1.ne.-1) then
               if(o1.eq.o3) then
                  d_covered_length=d_covered_length-1
               end if
            end if
         else
            d_length=d_length+1
            if(o1.ne.-1) then
               if(o1.eq.o3) then
                  d_covered_length=d_covered_length+1
               end if
            end if
         end if
         
         s3=sp(ix,mod(iy-1+Ly,Ly))
         a3=anch(ix,mod(iy-1+Ly,Ly))
         o3=occ(ix,mod(iy-1+Ly,Ly))
         dEn=dEn+De_bond(s1,a1,s3,a3,jnn,ja,jb,k)
         if(s1*s3 .lt. 0) then
            d_length=d_length-1
            if(o1.ne.-1) then
               if(o1.eq.o3) then
                  d_covered_length=d_covered_length-1
               end if
            end if
         else
            d_length=d_length+1
            if(o1.ne.-1) then
               if(o1.eq.o3) then
                  d_covered_length=d_covered_length+1
               end if
            end if
         end if
         
         s3=sp(ix,mod(iy+1,Ly))
         a3=anch(ix,mod(iy+1,Ly))
         o3=occ(ix,mod(iy+1,Ly))
         dEn=dEn+De_bond(s1,a1,s3,a3,jnn,ja,jb,k)
         if(s1*s3 .lt. 0) then
            d_length=d_length-1
            if(o1.ne.-1) then
               if(o1.eq.o3) then
                  d_covered_length=d_covered_length-1
               end if
            end if
         else
            d_length=d_length+1
            if(o1.ne.-1) then
               if(o1.eq.o3) then
                  d_covered_length=d_covered_length+1
               end if
            end if
         end if
         
         if(Mod(iy,2).eq.0) then 
            s3=sp(mod(ix-1+Lx,Lx),Mod(iy+1,Ly))
            a3=anch(mod(ix-1+Lx,Lx),Mod(iy+1,Ly))
            o3=occ(mod(ix-1+Lx,Lx),Mod(iy+1,Ly))
         else
            s3=sp(mod(ix+1,Lx),Mod(iy+1,Ly))
            a3=anch(mod(ix+1,Lx),Mod(iy+1,Ly))
            o3=occ(mod(ix+1,Lx),Mod(iy+1,Ly))
         end if 
         dEn=dEn+De_bond(s1,a1,s3,a3,jnn,ja,jb,k)
         if(s1*s3 .lt. 0) then
            d_length=d_length-1
            if(o1.ne.-1) then
               if(o1.eq.o3) then
                  d_covered_length=d_covered_length-1
               end if
            end if
         else
            d_length=d_length+1
            if(o1.ne.-1) then
               if(o1.eq.o3) then
                  d_covered_length=d_covered_length+1
               end if
            end if
         end if
         
         if(Mod(iy,2).eq.0) then 
            s3=sp(mod(ix-1+Lx,Lx),Mod(iy-1+Ly,Ly))
            a3=anch(mod(ix-1+Lx,Lx),Mod(iy-1+Ly,Ly))
            o3=occ(mod(ix-1+Lx,Lx),Mod(iy-1+Ly,Ly))
         else
            s3=sp(mod(ix+1,Lx),Mod(iy-1+Ly,Ly))
            a3=anch(mod(ix+1,Lx),Mod(iy-1+Ly,Ly))
            o3=occ(mod(ix+1,Lx),Mod(iy-1+Ly,Ly))
         end if         
         dEn=dEn+De_bond(s1,a1,s3,a3,jnn,ja,jb,k)
         if(s1*s3 .lt. 0) then
            d_length=d_length-1
            if(o1.ne.-1) then
               if(o1.eq.o3) then
                  d_covered_length=d_covered_length-1
               end if
            end if
         else
            d_length=d_length+1
            if(o1.ne.-1) then
               if(o1.eq.o3) then
                  d_covered_length=d_covered_length+1
               end if
            end if
         end if
         
c. we now update ip to ipp 
         sp(ix,iy)=s2
         
c. neighbours of ixp iyp
         s3=sp(mod(ixp-1+Lx,Lx),iyp)
         a3=anch(mod(ixp-1+Lx,Lx),iyp)
         o3=occ(mod(ixp-1+Lx,Lx),iyp)
         dEn=dEn+De_bond(s2,a2,s3,a3,jnn,ja,jb,k)
         if(s2*s3 .lt. 0) then
            d_length=d_length-1
            if(o2.ne.-1) then
               if(o2.eq.o3) then
                  d_covered_length=d_covered_length-1
               end if
            end if
         else
            d_length=d_length+1
            if(o2.ne.-1) then
               if(o2.eq.o3) then
                  d_covered_length=d_covered_length+1
               end if
            end if
         end if
         
         s3=sp(mod(ixp+1,Lx),iyp)
         a3=anch(mod(ixp+1,Lx),iyp)
         o3=occ(mod(ixp+1,Lx),iyp)
         dEn=dEn+De_bond(s2,a2,s3,a3,jnn,ja,jb,k)
         if(s2*s3 .lt. 0) then
            d_length=d_length-1
            if(o2.ne.-1) then
               if(o2.eq.o3) then
                  d_covered_length=d_covered_length-1
               end if
            end if
         else
            d_length=d_length+1
            if(o2.ne.-1) then
               if(o2.eq.o3) then
                  d_covered_length=d_covered_length+1
               end if
            end if
         end if
         
         s3=sp(ixp,mod(iyp-1+Ly,Ly))
         a3=anch(ixp,mod(iyp-1+Ly,Ly))
         o3=occ(ixp,mod(iyp-1+Ly,Ly))
         dEn=dEn+De_bond(s2,a2,s3,a3,jnn,ja,jb,k)
         if(s2*s3 .lt. 0) then
            d_length=d_length-1
            if(o2.ne.-1) then
               if(o2.eq.o3) then
                  d_covered_length=d_covered_length-1
               end if
            end if
         else
            d_length=d_length+1
            if(o2.ne.-1) then
               if(o2.eq.o3) then
                  d_covered_length=d_covered_length+1
               end if
            end if
         end if
         
         s3=sp(ixp,mod(iyp+1,Ly))
         a3=anch(ixp,mod(iyp+1,Ly))
         o3=occ(ixp,mod(iyp+1,Ly))
         dEn=dEn+De_bond(s2,a2,s3,a3,jnn,ja,jb,k)
         if(s2*s3 .lt. 0) then
            d_length=d_length-1
            if(o2.ne.-1) then
               if(o2.eq.o3) then
                  d_covered_length=d_covered_length-1
               end if
            end if
         else
            d_length=d_length+1
            if(o2.ne.-1) then
               if(o2.eq.o3) then
                  d_covered_length=d_covered_length+1
               end if
            end if
         end if
         
         if(Mod(iyp,2).eq.0) then 
            s3=sp(mod(ixp-1+Lx,Lx),Mod(iyp+1,Ly))
            a3=anch(mod(ixp-1+Lx,Lx),Mod(iyp+1,Ly))
            o3=occ(mod(ixp-1+Lx,Lx),Mod(iyp+1,Ly))
         else
            s3=sp(mod(ixp+1,Lx),Mod(iyp+1,Ly))
            a3=anch(mod(ixp+1,Lx),Mod(iyp+1,Ly))
            o3=occ(mod(ixp+1,Lx),Mod(iyp+1,Ly))
         end if
         dEn=dEn+De_bond(s2,a2,s3,a3,jnn,ja,jb,k)
         if(s2*s3 .lt. 0) then
            d_length=d_length-1
            if(o2.ne.-1) then
               if(o2.eq.o3) then
                  d_covered_length=d_covered_length-1
               end if
            end if
         else
            d_length=d_length+1
            if(o2.ne.-1) then
               if(o2.eq.o3) then
                  d_covered_length=d_covered_length+1
               end if
            end if
         end if
         
         if(Mod(iyp,2).eq.0) then 
            s3=sp(mod(ixp-1+Lx,Lx),Mod(iyp-1+Ly,Ly))
            a3=anch(mod(ixp-1+Lx,Lx),Mod(iyp-1+Ly,Ly))
            o3=occ(mod(ixp-1+Lx,Lx),Mod(iyp-1+Ly,Ly))
         else
            s3=sp(mod(ixp+1,Lx),Mod(iyp-1+Ly,Ly))
            a3=anch(mod(ixp+1,Lx),Mod(iyp-1+Ly,Ly))
            o3=occ(mod(ixp+1,Lx),Mod(iyp-1+Ly,Ly))
         end if 
         dEn=dEn+De_bond(s2,a2,s3,a3,jnn,ja,jb,k)
         if(s2*s3 .lt. 0) then
            d_length=d_length-1
            if(o2.ne.-1) then
               if(o2.eq.o3) then
                  d_covered_length=d_covered_length-1
               end if
            end if
         else
            d_length=d_length+1
            if(o2.ne.-1) then
               if(o2.eq.o3) then
                  d_covered_length=d_covered_length+1
               end if
            end if
         end if
      
         if(ran2(Iseed) .lt. dexp(-dEn) ) then

            acc = 1 
            sp(ixp,iyp)=s1
            EnCk = EnCk + dEn
            length = length + d_length
            length_covered = length_covered+ d_covered_length
         else
         
            sp(ix,iy)=s1
         
         end if

      end if
      

      end

      
c. ...... .c

c..........!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!..........c
c.                                                   .c
c. single update has not been adapted to the new BCs .c
c.                                                   .c
c.....................................................c
      
      subroutine single_update(acc)

      implicit none
      include 'parameter.inc'
      include 'system.inc'
      integer ix, iy, iym
      integer ixp, iyp, ix2, ixm
      integer ip, ipp, ingh
      integer ActSp
      double precision ran2, rn
      double precision dEn, De_bond 
      integer acc, s1, s2, a1, a2
      integer k
      integer d_length

      acc = 0
      k=1
 
      iy=int(ran2(Iseed)*dble(Ly))
      if(fixed_boundary .eqv. .true.) then
         ix=int(ran2(Iseed)*dble(Lx-2))+1
      else
         ix=int(ran2(Iseed)*dble(Lx))
      end if

      ixm=mod(ix-1+Lx,Lx)
      ixp=mod(ix+1,Lx)
      iym=mod(iy-1+Ly,Ly)
      iyp=mod(iy+1,Ly)
      ix2=ixp
      if(mod(iy,2).eq.0) ix2=ixm
      s1=sp(ix,iy)
      a1=anch(ix,iy)

      s2=sp(ixp,iy)
      a2=anch(ixp,iy)
      dEn=De_bond(s1,a1,s2,a2,jnn,ja,jb,k)
      if(s1*s2 .lt. 0) then
         d_length=-1
      else
         d_length=1
      end if
      
      s2=sp(ixm,iy)
      a2=anch(ixm,iy)
      dEn=dEn+De_bond(s1,a1,s2,a2,jnn,ja,jb,k)
      if(s1*s2 .lt. 0) then
         d_length=d_length-1
      else
         d_length=d_length+1
      end if
      
      s2=sp(ix,iyp)
      a2=anch(ix,iyp)
      dEn=dEn+De_bond(s1,a1,s2,a2,jnn,ja,jb,k)
      if(s1*s2 .lt. 0) then
         d_length=d_length-1
      else
         d_length=d_length+1
      end if
      
      s2=sp(ix2,iyp)
      a2=anch(ix2,iyp)
      dEn=dEn+De_bond(s1,a1,s2,a2,jnn,ja,jb,k)
      if(s1*s2 .lt. 0) then
         d_length=d_length-1
      else
         d_length=d_length+1
      end if
      
      s2=sp(ix,iym)
      a2=anch(ix,iym)
      dEn=dEn+De_bond(s1,a1,s2,a2,jnn,ja,jb,k)
      if(s1*s2 .lt. 0) then
         d_length=d_length-1
      else
         d_length=d_length+1
      end if
      
      s2=sp(ix2,iym)
      a2=anch(ix2,iym)
      dEn=dEn+De_bond(s1,a1,s2,a2,jnn,ja,jb,k)
      if(s1*s2 .lt. 0) then
         d_length=d_length-1
      else
         d_length=d_length+1
      end if

      if( ran2(Iseed) .lt. dexp(-dEn) ) then
          acc = 1 
          sp(ix,iy) = -sp(ix,iy)
          EnCk = EnCk + dEn
          length = length + d_length
       end if

       end subroutine

c     . to support energy calculation
c     anchor-anchor interaction=0

      function e_bond(sp, ap, spp, app, jnn, ja, jb, k)

      implicit none
      double precision e_bond
      double precision jnn, ja, jb
      integer sp, spp
      integer ap, app

      integer k

c. sp <=> spp     
      e_bond=-jnn*sp*spp

c. (ap <=> spp) + (ap <=> sp)*k/6     
      if(ap .eq. 2) then
         e_bond=e_bond-ja*ap*spp/2
         e_bond=e_bond-ja*ap*sp*k/12
      end if
      if(ap .eq. 3) then
         e_bond=e_bond-jb*ap*spp/3
         e_bond=e_bond-jb*ap*sp*k/18
      end if

c. app <=> sp      
      if(app .eq. 2) then
         e_bond=e_bond-ja*app*sp/2
      end if
      if(app .eq. 3) then
         e_bond=e_bond-jb*app*sp/3
      end if

      end function 

c. when swapping sp => sp=+-1
      function De_bond(sp, ap, spp, app, jnn, ja, jb, k)

      implicit none
      double precision De_bond
      double precision e_bond
      double precision jnn, ja, jb
      integer sp, spp
      integer ap, app
      integer k

      De_bond=e_bond(-sp, ap, spp, app, jnn, ja, jb, k)

      De_bond=De_bond
     &     -e_bond(sp, ap, spp, app, jnn, ja, jb, k)

      if (abs(sp) .ne. 1) then
         write(*,*) 'problem in De_bond'
      end if
      
c      if(abs(spp) .eq. 1) then
c         De_bond=2*jnn*sp*spp
c      end if

c      if(abs(sp*spp) .eq. 2) then
c         De_bond=2*ja*sp*spp/2
c      end if

c      if(abs(sp*spp) .eq. 3) then
c         De_bond=2*jb*sp*spp/3
c      end if

      end function 

c. statistics of the number of tiles at the CL .c
      subroutine update_n_fav_cont()

      implicit none
      include 'parameter.inc'
      include 'system.inc'
      integer it, ian 
      integer iDx, iDy, ix, iy, s1, s2
      integer a1, a2
      integer ixp, iyp
      integer ix_CM, iy_CM, ithe
      integer ifav
      double precision en, janch 

      do it=0,N_tiles-1
         ix_CM=ixiy_tiles(0,it)
         iy_CM=ixiy_tiles(1,it)
         ithe=theta_tiles(it)
         ifav=0
         do ian=0,N_modules-1
            iDx=mask_tile(0,list_modules(ian))
            iDy=mask_tile(1,list_modules(ian))
            call calc_ix_iy(ix_CM,iy_CM,
     &           iDx,iDy,ithe,ix,iy)
            a1=anch(ix,iy)
            if(a1 .eq. 0) write(*,*) 'inconsistency in
     &         update_n_fav_cont'
            if(a1 .eq. 2) janch = ja
            if(a1 .eq. 3) janch = jb
            
            s1=sp(ix,iy)
            en=-janch*s1
            if(en.lt.0) ifav=ifav+1

            ixp=Mod(ix+1,Lx)
            s2=sp(ixp,iy)
            en=-janch*s2
            if(en.lt.0) ifav=ifav+1

            ixp=Mod(ix-1+Lx,Lx)
            s2=sp(ixp,iy)
            en=-janch*s2
            if(en.lt.0) ifav=ifav+1

            iyp=Mod(iy+1,Ly)
            s2=sp(ix,iyp)
            en=-janch*s2
            if(en.lt.0) ifav=ifav+1

            iyp=Mod(iy-1+Ly,Ly)
            s2=sp(ix,iyp)
            en=-janch*s2
            if(en.lt.0) ifav=ifav+1

            if(Mod(iy,2) .eq. 0) ixp=Mod(ix-1+Lx,Lx)
            if(Mod(iy,2) .eq. 1) ixp=Mod(ix+1,Lx)
            iyp=Mod(iy+1,Ly)
            s2=sp(ixp,iyp)
            en=-janch*s2
            if(en.lt.0) ifav=ifav+1

            iyp=Mod(iy-1+Ly,Ly)
            s2=sp(ixp,iyp)
            en=-janch*s2
            if(en.lt.0) ifav=ifav+1

         end do
         
         n_favorable_contacts(ifav)=n_favorable_contacts(ifav)+1.0
         N_samples=N_samples+1
         
      end do
           
      end subroutine 
      
      function ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      DOUBLE PRECISION ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1.d0/IM1,IMM1=IM1-1, 
     & IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, 
     & NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2d-7,RNMX=1.d0-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      end
! (C) Copr. 1986-92 Numerical Recipes Software 1(-V%'2150)-3. ccccccccc



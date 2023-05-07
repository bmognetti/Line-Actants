
      program line_actants

      implicit none
      include 'parameter.inc'
      include 'system.inc'

c. ... ... ... ... ... ... ... ... ... ... .c

      double precision RANnyu
      double precision rnd, ran2
      integer i1, i2, i3
      integer ix, iy, iz
      integer iout, iout1
      integer isp
      integer ixM, iyM
      integer Dx, Dy

      double precision En

      integer NMC, NSamp
      integer acc, Nacc_s, Ntot_s
      integer Nacc_d, Ntot_d      
      integer n_tiles_x, n_tiles_y 

      integer length_ck, length_covered_ck
      integer length_max

      double precision p_spin

      double precision p_phi

      logical ovl, overlap, local 
      
      open(5,file = "config.dat",action="write",position="append")
      open(17,file="config_spin.dat",action="write",position="append")
      open(10,file="config_tiles.dat",action="write",
     &     position="append")
      open(11,file="config_tiles_2.dat",action="write",
     &     position="append")
      open(25,file="system.dat")
      open(35,file="MCdetails.dat",action="write",
     &     position="append")
      open(46,file="tiles_at_CL.dat",action="write",
     &     position="append")
      open(47,file="length.dat",action="write",
     &     position="append")
      open(51,file="length_covered.dat",action="write",
     &     position="append")
      open(48,file="spins.dat")
      open(49,file="tiles.dat")

c. ... ... ... ... .c

      read(25,*) Lx, Ly, fixed_boundary
      read(25,*) N_tiles, n_tiles_x, n_tiles_y 
      read(25,*) jnn
      read(25,*) ja, jb
      read(25,*) NMC, iout
      read(25,*) Iseed, p_phi

      iout1=iout

      if(mod(Ly,2).eq.1) then
         write(35,*) 'Ly must be even'
         write(*,*) 'Ly must be even'
         return
      end if

c      if(mod(n_tiles_x,2).eq.0) then
c         write(35,*) 'n_tiles_x must be odd'
c         write(*,*) 'n_tiles_x must be odd'
c         return
c      end if

c      if(mod(n_tiles_y,2).eq.0) then
c         write(35,*) 'n_tiles_y must be odd'
c         write(*,*) 'n_tiles_y must be odd'
c         return
c      end if

      rnd=ran2(-Iseed)
      do i1=1,5
         write(35,*) ran2(Iseed)
      end do

c. initialise the mask with the tile's shape                       .c
c.   variable initialised= N_tiles, N_sites, N_modules             .c
c.                      mask_tile,list_modules,type_modules        .c

      call load_mask(n_tiles_x,n_tiles_y)

c.      do i1=0,N_sites-1
c.         write(*,*) i1, mask_tile(0,i1), mask_tile(1,i1)
c.      end do

cccccccccccccccccccccccccccccccc
c                              cc
c initialisation of the system ccc
c                              cc
cccccccccccccccccccccccccccccccc

      if(iout.eq.0) then 

cccccccccccccccccccccccccccccccccc
c. random initial configuration .c
cccccccccccccccccccccccccccccccccc

c      do i2=0,Lx-1
c         do i3=0,Ly-1

c            if(ran2(Iseed) .lt. p_phi) then
c               sp(i2,i3)=-1
c            else
c               sp(i2,i3)=1
c            end if
c            anch(i2,i3)=0
c            occ(i2,i3)=-1
cc. set up BC .c
c            if(fixed_boundary .eqv. .true.) then

cc     . new BCs here  i2 => i3 in if

c               if(i3 .eq. 0) sp(i2,i3)=-1 
c               if(i3 .eq. Ly-1) sp(i2,i3)=1

c            end if
c         end do
c      end do     

cccccccccccccccccccccccccccccccccccccccccccc
c.  striped BC, consistent with fixed BCs  c
cccccccccccccccccccccccccccccccccccccccccccc

      do i2=0,Lx-1
         do i3=0,Ly-1
            if(i3 .lt. Ly/2) then
               sp(i2,i3)=-1
            else
               sp(i2,i3)=1
            end if
            anch(i2,i3)=0
            occ(i2,i3)=-1
         end do
      end do     

c. initialise the position and the orientation of the N_tiles .c

      i1=0
      do while (i1 .le. N_tiles-1)
         ix=int(ran2(Iseed)*dble(Lx))
         iy=int(ran2(Iseed)*dble(Ly))
         i3=int(ran2(Iseed)*6)
cc.      debug
c         ix=Lx/2
c         iy=Ly/2
c         i3=0
cc.      end debug 
         ovl=overlap(ix,iy,i3,i1)
         if(ovl.eqv. .false.) then
            ixiy_tiles(0,i1)=ix
            ixiy_tiles(1,i1)=iy
            theta_tiles(i1)=i3
c.      occupy update occ() and anch()
            call occupy(ix,iy,i3,i1)
            i1=i1+1
c. the bloc below is now (sp+anchor version) done in occupy, sp always +-1 
c     we now insert the anchors/modules
c            do i2=0,N_modules-1
c               Dx=mask_tile(0,list_modules(i2))
c               Dy=mask_tile(1,list_modules(i2))
c               call calc_ix_iy(ix,iy,Dx,Dy,i3,ixM,iyM)
c               if(type_modules(i2).eq.2) sp(ixM,iyM)=2
c               if(type_modules(i2).eq.3) sp(ixM,iyM)=3
c            end do
            write (*,*) 'tiles inserted=', i1
         end if
      end do

c. ... termalisation ... ... ... ... ... ... ... .c
c      do i1=1,NMC/10
c        do i2=1,Lx*Ly
c          call single_update(acc)
c        end do
c     end do
      iout=0
      iout1=0
      call PrintOut(iout)
      call PrintOut_tiles(iout1)
      iout=iout+1
      iout1=iout1+1

cccccccccccccccccccccccccccccccccccccccccccccc
c                                            cc
c reading starting configurations from files ccc
c                                            cc
cccccccccccccccccccccccccccccccccccccccccccccc      
      else

c        reading the spin configuration 
         do i2=0,Lx-1
            do i3=0,Ly-1
               read(48,*) isp
               sp(i2,i3) = isp
               anch(i2,i3)=0
               occ(i2,i3)=-1
            end do
         end do

c        loading the tile configuration
         i1=0
         do while (i1 .le. N_tiles-1)
            read(49,*) isp, ix, iy, i3

            if(isp .ne. i1) then
               write(*,*) "inconsistency when loading tiles, e1"
               return
            end if
            ovl=overlap(ix,iy,i3,i1)
            if(ovl .eqv. .true.) then
               write(*,*) "inconsistency when loading tiles, e2"
               return
            end if
            ixiy_tiles(0,i1)=ix
            ixiy_tiles(1,i1)=iy
            theta_tiles(i1)=i3
            call occupy(ix,iy,i3,i1)
            i1=i1+1
         end do
         iout=iout+1
         iout1=iout

      end if

c.  ....  ....  ....  ....  ....  ....  ....  ....  ....  .c

c     . initialisation of the observables
      do i1=0,N_sites_max
         n_favorable_contacts(i1)=0
      end do
      length_max=Lx*Ly
      do i1=0,length_max
         length_distribution(i1)=0
         length_covered_distribution(i1)=0
      end do
      N_samples=0
      N_samples_length=0
      N_samples_covered_length=0
      Nacc_s=0
      Ntot_s=0
      Nacc_d=0
      Ntot_d=0

      call TotEn(En,length_ck,length_covered_ck)
      EnCk=En
      length=length_ck
      length_covered=length_covered_ck

      i1=0
      write(15,*) i1, En, length_ck,length_covered_ck

c      write(*,*) i1, En, length_ck,length_covered_ck

      
      p_spin=0.9
      if(N_tiles .eq. 0) p_spin=1.1

ccccccccccccccccccccc
cc                 cc
cc   main MC loop  cc 
cc                 cc
ccccccccccccccccccccc

      local = .true.
c      local = .false.

      do i1=1,NMC+1

         do i2=1,Lx*Ly+2*N_tiles

            if(ran2(Iseed).lt.p_spin) then
               call KawasakiSwap(acc)
               Nacc_s=Nacc_s+acc
               Ntot_s=Ntot_s+1
            else
c               if(ran2(Iseed).lt.0.98) then
c                  local=.true.
c               else
c                  local=.false.
c               end if
               call displace_tiles(local,acc)
               Nacc_d=Nacc_d+acc
               Ntot_d=Ntot_d+1
            end if

         end do

c        do i2=1,Lx*Ly
c          call single_update(acc) 
c        end do

c        do i2=0,N_tiles-1
c           call displace_tiles(local,acc)
c        end do

c     . debug
c        call TotEn(En)
c        write(*,*) En, EnCk
c        call PrintOut(1)
c        return
c     . end debug

         if(mod(i1*50,NMC).eq.0) then

            call PrintOut(iout)
            call PrintOut_tiles_CL(iout)
            call PrintOut_length_distribution(iout,length_max)
            iout=iout+1
            call TotEn(En,length_ck,length_covered_ck)
            write(35,*) En, EnCk, length, length_ck,
     &           length_covered, length_covered_ck
            write(*,*) i1
c     ja=ja+0.1
c           jb=jb-0.1
         end if

         if(mod(i1*50,NMC).eq.0) then
            call PrintOut_tiles(iout1)
            iout1=iout1+1
         end if
      
         call update_n_fav_cont()
        if (length .lt. length_max) then
           length_distribution(length)=
     &          length_distribution(length)+1
           N_samples_length=N_samples_length+1
        end if
        if (length_covered .lt. length_max) then
           length_covered_distribution(length_covered)=
     &     length_covered_distribution(length_covered)+1
           N_samples_covered_length=N_samples_covered_length+1
        end if
        
      end do

      write(35,*) 'acceptances'
      write(35,*) 'spin updates (Nmov, accepted, %)'
      write(35,*) Ntot_s, Nacc_s, 1.*Nacc_s/Ntot_s
      write(35,*) 'tile updates (Nmov, accepted, %)'
      write(35,*) Ntot_d, Nacc_d, 1.*Nacc_d/Ntot_d

      close(5)
      close(17)
      close(10)
      close(11)
      close(25)
      close(35)
      close(46)
      close(47)
      close(48)
      close(49)

      end

cccccccccccccccccccccccccccc
c  end of the main program c
cccccccccccccccccccccccccccc

c      subroutine SaveConf()
c      implicit none
c      include 'parameter.inc'
c      include 'system.inc'
c      integer ix, iy, iz
c      integer i1
c      open(1005,file='startconf.dat')
c      do i1=1,Np
c        write(1005,*) x(i1), y(i1), z(i1), the(i1)
c      end do
c      close(1005)
c      end

c      subroutine LoadConf()
c      implicit none 
c      include 'parameter.inc'
c      include 'system.inc'
c      integer ix, iy, iz
c      integer ixp, iyp, izp
c      integer ip

ccccccccccccccccccccccccccccc
c total energy              c
ccccccccccccccccccccccccccccc

      subroutine TotEn(En,ln,ln_cv)

      implicit none 
      include 'parameter.inc'
      include 'system.inc'

      double precision En, e_bond
      integer  ix, iy, ixmin, ixmax, iymin, iymax
      integer ixp,iyp
      integer s1,s2,a1,a2,o1,o2
c. k/6 is the number of times the sp-ap interaction is counted
      integer k
      integer ln, ln_cv

      ixmin=0
      ixmax=Lx-1
      iymin=0
      iymax=Ly-1    
      k=2

      if(fixed_boundary .eqv. .true.) then
c. new BC
         iymin=0
         iymax=Ly-2
      end if

      En = 0.d0
      ln = 0
      ln_cv = 0
      
      do ix=ixmin,ixmax
c     . new BCs
         do iy=iymin,iymax
c     do iy=0,Ly-1
            s1=sp(ix,iy)
            a1=anch(ix,iy)
            o1=occ(ix,iy) 
            
            ixp=Mod(ix+1,Lx)
c     En=En-jnn*sp(ix,iy)*sp(ixp,iy)
            s2=sp(ixp,iy)
            o2=occ(ixp,iy) 
            a2=anch(ixp,iy)
            En=En+e_bond(s1,a1,s2,a2,jnn,ja,jb,k)
            if(s1*s2.lt.0) then
               ln=ln+1
               if(o1.eq.o2) then
                  if(o1 .ne. -1) then
                     ln_cv=ln_cv+1
                  end if
               end if
            end if
            iyp=Mod(iy+1,Ly)
c     En=En-jnn*sp(ix,iy)*sp(ix,iyp)
            s2=sp(ix,iyp)
            a2=anch(ix,iyp)
            o2=occ(ix,iyp)
            En=En+e_bond(s1,a1,s2,a2,jnn,ja,jb,k)
            if(s1*s2.lt.0) then
               ln=ln+1
               if(o1.eq.o2) then
                  if(o1 .ne. -1) then
                     ln_cv=ln_cv+1
                  end if
               end if
            end if
            
            if(Mod(iy,2) .eq. 0) ixp=Mod(ix-1+Lx,Lx)
            if(Mod(iy,2) .eq. 1) ixp=Mod(ix+1,Lx)
c     En=En-jnn*sp(ix,iy)*sp(ixp,iyp)
            s2=sp(ixp,iyp)
            a2=anch(ixp,iyp)
            o2=occ(ixp,iyp)
            En=En+e_bond(s1,a1,s2,a2,jnn,ja,jb,k)
            if(s1*s2.lt.0) then
               ln=ln+1
               if(o1.eq.o2) then
                  if(o1 .ne. -1) then
                     ln_cv=ln_cv+1
                  end if
               end if
            end if           
         end do
      end do

c      if(fixed_boundary .eqv. .true.) then
c         ix=1
c         ixp=0
c         k=0
c         do iy=0,Ly-1
c            s1=sp(ix,iy)
c            a1=anch(ix,iy)
c     En=En-jnn*sp(ix,iy)*sp(ixp,iy)
c            s2=sp(ixp,iy)
c            a2=anch(ixp,iy)
c            En=En+e_bond(s1,a1,s2,a2,jnn,ja,jb,k)
c            if(s1*s2.lt.0) ln=ln+1
            
c            if(Mod(iy,2) .eq. 0) then     
c               iyp=Mod(iy-1+Ly,Ly)
c     En=En-jnn*sp(ix,iy)*sp(ixp,iyp)
c               s2=sp(ixp,iyp)
c               a2=anch(ixp,iyp)
c               En=En+e_bond(s1,a1,s2,a2,jnn,ja,jb,k)
c               if(s1*s2.lt.0) ln=ln+1
c            end if 
c         end do
c         ix=Lx-2
c         ixp=Lx-1
c         do iy=0,Ly-1
c            s1=sp(ix,iy)
c            a1=anch(ix,iy)
c            if(Mod(iy,2) .eq. 1) then     
c               iyp=Mod(iy-1+Ly,Ly)
cc     En=En-jnn*sp(ix,iy)*sp(ixp,iyp)
c               s2=sp(ixp,iyp)
c               a2=anch(ixp,iyp)
c               En=En+e_bond(s1,a1,s2,a2,jnn,ja,jb,k)
c               if(s1*s2.lt.0) ln=ln+1
c            end if
c         end do
c      end if

      end

c.                              .c 
c.   print out configurations   .c cc c
c.                              .c

      subroutine PrintOut(n)

      implicit none 
      include 'parameter.inc'
      include 'system.inc'
      integer ix, iy, n, iprint, itile
      integer ixcm, iycm, ior, ii, iii
      double precision xcm, ycm
      double precision Vec_rot(0:1,0:3), vx, vy

      do ix=0,Lx-1
         do iy=0,Ly-1
            iprint= sp(ix,iy)
            write(17,*) n, ix+0.5*Mod(iy,2), 
     &                iy*0.8660254038,       
     &                iprint
            if(anch(ix,iy) .ne. 0) iprint=anch(ix,iy)
            write(5,*) n, ix+0.5*Mod(iy,2), 
     &                iy*0.8660254038,       
     &                iprint
        end do
      end do

      do ix=0,Lx-1
         do iy=0,Ly-1
c            write(*,*) occ(ix,iy)
            if(occ(ix,iy) .ne. -1) then               
               write(10,*) n, ix+0.5*Mod(iy,2), 
     &                iy*0.8660254038
            end if   
         end do
      end do

c      do itile=0, N_tiles-1
         
c         ixcm=ixiy_tiles(0,itile)
c         iycm=ixiy_tiles(1,itile)
c         xcm=ixcm+0.5*Mod(iycm,2)
c         ycm=iycm*0.8660254038
         
c         ior=theta_tiles(itile)

c         Vec_rot(0,0)=Vec_corner(0,0)
c         Vec_rot(1,0)=Vec_corner(1,0)
c         Vec_rot(0,1)=Vec_corner(0,1)
c         Vec_rot(1,1)=Vec_corner(1,1)
c         Vec_rot(0,2)=Vec_corner(0,2)
c         Vec_rot(1,2)=Vec_corner(1,2)
c         Vec_rot(0,3)=Vec_corner(0,3)
c         Vec_rot(1,3)=Vec_corner(1,3)

c         ii=1
c         do while(ii.le.ior)
c            do iii=0,3
c               vx=Vec_rot(0,iii)*0.5
c     &              -Vec_rot(1,iii)*0.8660254038
c               vy=Vec_rot(0,iii)*0.8660254038
c     &              +Vec_rot(1,iii)*0.5
c               Vec_rot(0,iii)=vx 
c               Vec_rot(1,iii)=vy 
c            end do
c            ii=ii+1
c         end do
c         do iii=0,3
c               Vec_rot(0,iii)=Vec_rot(0,iii)+xcm
c               Vec_rot(1,iii)=Vec_rot(1,iii)+ycm
c            end do
c           
c         write(11,*) n,itile,ixcm,iycm,ior,Vec_rot(0,0),Vec_rot(1,0),  
c     &           Vec_rot(0,1), Vec_rot(1,1), Vec_rot(0,2), Vec_rot(1,2),
c     &           Vec_rot(0,3), Vec_rot(1,3)
c      end do 

      end


      subroutine PrintOut_tiles(n)

      implicit none 
      include 'parameter.inc'
      include 'system.inc'
      integer ix, iy, n, iprint, itile
      integer ixcm, iycm, ior, ii, iii
      double precision xcm, ycm
      double precision Vec_rot(0:1,0:3), vx, vy
      

      do itile=0, N_tiles-1
         
         ixcm=ixiy_tiles(0,itile)
         iycm=ixiy_tiles(1,itile)
         xcm=ixcm+0.5*Mod(iycm,2)
         ycm=iycm*0.8660254038
         
         ior=theta_tiles(itile)

         Vec_rot(0,0)=Vec_corner(0,0)
         Vec_rot(1,0)=Vec_corner(1,0)
         Vec_rot(0,1)=Vec_corner(0,1)
         Vec_rot(1,1)=Vec_corner(1,1)
         Vec_rot(0,2)=Vec_corner(0,2)
         Vec_rot(1,2)=Vec_corner(1,2)
         Vec_rot(0,3)=Vec_corner(0,3)
         Vec_rot(1,3)=Vec_corner(1,3)

         ii=1
         do while(ii.le.ior)
            do iii=0,3
               vx=Vec_rot(0,iii)*0.5
     &              -Vec_rot(1,iii)*0.8660254038
               vy=Vec_rot(0,iii)*0.8660254038
     &              +Vec_rot(1,iii)*0.5
               Vec_rot(0,iii)=vx 
               Vec_rot(1,iii)=vy 
            end do
            ii=ii+1
         end do
         do iii=0,3
               Vec_rot(0,iii)=Vec_rot(0,iii)+xcm
               Vec_rot(1,iii)=Vec_rot(1,iii)+ycm
            end do
           
         write(11,*) n,itile,ixcm,iycm,ior,Vec_rot(0,0),Vec_rot(1,0),  
     &           Vec_rot(0,1), Vec_rot(1,1), Vec_rot(0,2), Vec_rot(1,2),
     &           Vec_rot(0,3), Vec_rot(1,3)
      end do
      
      end

      
      subroutine PrintOut_tiles_CL(n)

      implicit none 
      include 'parameter.inc'
      include 'system.inc'
      integer i1, n

      do i1=0,N_sites_max
         if (i1 .le. 7*N_modules) then
            write(46,*) n,i1,
     &           1.*n_favorable_contacts(i1)/N_samples
            n_favorable_contacts(i1)=0
         else
            if(n_favorable_contacts(i1).gt.0) then
               write(35,*) "inconsistency in Print_Out_t_CL"
               write(*,*) "inconsistency in Print_Out_t_CL"
            end if
         end if   
      end do
      N_samples=0
      end

      subroutine PrintOut_length_distribution(n,length_max)

      implicit none 

      include 'parameter.inc'
      include 'system.inc'
      integer i1, n, length_max

      do i1=0,length_max
         write(47,*) n,i1,
     &   1.*length_distribution(i1)/N_samples_length,
     &   1.*length_covered_distribution(i1)/N_samples_covered_length     
         length_distribution(i1)=0
         length_covered_distribution(i1)=0
      end do
      N_samples_length=0
      N_samples_covered_length=0

      end
      
c ********************************************************************
c     RANDOM NUMBER GENERATOR RANNYU
c     SET SEED BY RANSET, RECOVER BY RANGET
c ********************************************************************

c      function RANnyu()
c      double precision  twom12,rannyu
c      parameter (twom12 = 1/4096.d0)
c      common /rnyucm/ m1,m2,m3,m4,l1,l2,l3,l4

cc     This is rannyu as modified by A. Sokal 9/26/85.
cc        It is linear congruential with modulus m = 2**48, increment c = 1,
cc        and multiplier a = (2**36)*m1 + (2**24)*m2 + (2**12)*m3 + m4.
cc        The multiplier is stored in common (see subroutine setrn)
cc      and is set to a = 31167285 (recommended by Knuth, vol. 2,
cc      2nd ed., p. 102).

c      i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1
c      i2 = l2*m4 + l3*m3 + l4*m2
c      i3 = l3*m4 + l4*m3
c      i4 = l4*m4  +  1
c      l4 = mod(i4, 4096)
c      i3 = i3 + i4/4096
c      l3 = mod(i3, 4096)
c      i2 = i2 + i3/4096
c      l2 = mod(i2, 4096)
c      l1 = mod(i1 + i2/4096, 4096)
c      rannyu = twom12*(dble(l1) +
c     +       twom12*(dble(l2) +
c     +       twom12*(dble(l3) +
c     +       twom12*(dble(l4)))))
c      return
c      end

c **************************************************************

c      subroutine setrn(iseed)
c      common /rnyucm/ m(4),l(4)
c      integer iseed(4)
cc
cc     Multiplier is 31167285 = (2**24) + 3513*(2**12) + 821.
cc        Recommended by Knuth, vol. 2, 2nd ed., p. 102.
cc     (Generator is linear congruential with odd increment
cc        and maximal period, so seed is unrestricted: it can be
cc        either even or odd.)
cc
cc     data m /   0,   1,3513, 821/
cc     data l /   0,   0,   0,   1/
cc
c      do 10 i = 1, 4
c         l(i) = iseed(i)
c10    continue
c      m(1) = 0
c      m(2) = 1
c      m(3) = 3513
c      m(4) = 821
c      return
c      end
c
c **********************************************************
c
c      subroutine savern(iseed)
c      common /rnyucm/ m(4),l(4)
c      integer iseed(4)
cc
c      do 10 i = 1, 4
c         iseed(i) = l(i)
c10    continue
c      return
c      end

c     
c **********************************************************
c

     

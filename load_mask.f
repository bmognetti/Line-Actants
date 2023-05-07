
c.   initialise the mask with the tile's shape
c.   variable initialised= N_tiles, N_sites, N_modules
c.   mask_tile,list_modules,type_modules       

      subroutine load_mask(n_tiles_x,n_tiles_y)

      implicit none
      include 'parameter.inc'
      include 'system.inc'
      integer i1, it, i, ii, i_site
      integer ix, iy, Dix, Diy, ixCM, iyCM
      integer n_tiles_x, n_tiles_y
      integer nym, nyp, nxm, nxp
      double precision xx, yy, xCM, yCM

      N_sites=n_tiles_x*n_tiles_y

      N_sites_ev=-1
c      list_border(0:N_sites_Max-1)

      if(mod(n_tiles_x,2).eq.0) then
         nxm=-n_tiles_x/2+1
         nxp=n_tiles_x/2
      else
         nxm=-(n_tiles_x-1)/2
         nxp=(n_tiles_x-1)/2
      end if

      if(mod(n_tiles_y,2).eq.0) then
         nym=-n_tiles_y/2+1
         nyp=n_tiles_y/2
      else
         nym=-(n_tiles_y-1)/2
         nyp=(n_tiles_y-1)/2
      end if
      
      i_site=0
      do iy=nym, nyp
         ii=IABS(iy)
         if(iy .lt. 0) then
            ii=-ii
            ii=ii/2
         end if
         if(iy .gt. 0) then
            ii=(ii+1)/2
         end if
         do ix=nxm, nxp
            mask_tile(0,i_site) = ix-ii
            mask_tile(1,i_site) = iy
            if( (iy .eq. nym) .or.
     &          (iy .eq. nyp) .or.
     &          (ix .eq. nxm) .or.
     &          (ix .eq. nxp) ) then
               N_sites_ev=N_sites_ev+1
               list_border(N_sites_ev)=i_site
            end if
            i_site=i_site+1
c            write(*,*) ix-ii, iy
         end do
      end do

c debug
c      write(*,*) N_sites_ev
c      do i1=0,N_sites_ev
c         write(*,*) list_border(i1)
c      end do
      
      do i1=0,3
         Vec_corner(0,i1)=0.0
         Vec_corner(1,i1)=0.0
      end do

      xCM=Lx/2.+0.5*Mod(Ly/2,2)
      yCM=Ly/2.*0.8660254038

      do i1=0,N_sites-1

         Dix=mask_tile(0,i1)
         Diy=mask_tile(1,i1)
         call calc_ix_iy(Lx/2,Ly/2,Dix,Diy,0,ix,iy)
         xx=ix+0.5*Mod(iy,2)-xCM
         yy=iy*0.8660254038-yCM

         if(xx.lt.Vec_corner(0,1)) then
            Vec_corner(0,1)=xx
            Vec_corner(0,2)=xx
         end if
         if(xx.gt.Vec_corner(0,0)) then
            Vec_corner(0,0)=xx
            Vec_corner(0,3)=xx
         end if
         if(yy.lt.Vec_corner(1,2)) then
            Vec_corner(1,2)=yy
            Vec_corner(1,3)=yy
         end if
         if(yy.gt.Vec_corner(1,0)) then
            Vec_corner(1,0)=yy
            Vec_corner(1,1)=yy
         end if
      end do

      Vec_corner(0,0)=Vec_corner(0,0)+0.5
      Vec_corner(1,0)=Vec_corner(1,0)+0.8660254038/2
      
      Vec_corner(0,1)=Vec_corner(0,1)-0.5
      Vec_corner(1,1)=Vec_corner(1,1)+0.8660254038/2

      Vec_corner(0,2)=Vec_corner(0,2)-0.5
      Vec_corner(1,2)=Vec_corner(1,2)-0.8660254038/2

      Vec_corner(0,3)=Vec_corner(0,3)+0.5
      Vec_corner(1,3)=Vec_corner(1,3)-0.8660254038/2
  

c      do ii=0,3
c         write(*,*) Vec_corner(0,ii), Vec_corner(1,ii)
c      end do

      open(45,file="mask.dat")

      read(45,*) N_modules

      do i1=0,N_modules-1

         read(45,*) i_site, it
         list_modules(i1) = i_site
         type_modules(i1) = it

      end do

      close(45)

      end

c. this function check overlap between i_tile and the rest of the tiles and boundary conditions

      function overlap(iCM_new_x,iCM_new_y,the_new,i_tile)

      implicit none
      include 'parameter.inc'
      include 'system.inc'
      integer i1, it, i, ii, i_site
      integer ix, iy, Dix, Diy
      integer iCM_new_x,iCM_new_y,the_new,i_tile
      logical overlap

      overlap=.false.

c      do i1=0,N_sites-1
c         Dix=mask_tile(0,i1)
c         Diy=mask_tile(1,i1)
cc.        write (*,*) i1
c         call calc_ix_iy(iCM_new_x,iCM_new_y,Dix,Diy,the_new,ix,iy)
cc.         write (*,*) ix, iy
c         if (fixed_boundary .eqv. .true.) then
cc. new BCs
c            if ((iy .le. 0) .or. (iy .ge. Ly-1)) then
c               overlap=.true.
c               return
c            end if
c         end if
c         if (occ(ix,iy) .ne. -1) then
c            if (occ(ix,iy) .ne. i_tile) then
c               overlap=.true.
c               return
c            end if
c         end if
c      end do

      do ii=0,N_sites_ev
         i1=list_border(ii)
         Dix=mask_tile(0,i1)
         Diy=mask_tile(1,i1)
c.        write (*,*) i1
         call calc_ix_iy(iCM_new_x,iCM_new_y,Dix,Diy,the_new,ix,iy)
c.         write (*,*) ix, iy
         if (fixed_boundary .eqv. .true.) then
c. new BCs
            if ((iy .le. 0) .or. (iy .ge. Ly-1)) then
               overlap=.true.
               return
            end if
         end if
         if (occ(ix,iy) .ne. -1) then
            if (occ(ix,iy) .ne. i_tile) then
               overlap=.true.
               return
            end if
         end if
      end do

      end function

c. ... ... ... ... ... .c

      subroutine occupy(iCMx,iCMy,ithe,itile)
      implicit none
      include 'parameter.inc'
      include 'system.inc'
      integer iCMx,iCMy,itile,ithe
      integer Dix, Diy
      integer isite, ix, iy, i1, ln_cv
      integer ixp,iyp,o1,o2,s1,s2

c this could be greatly optimised in local moves      
      do isite=0,N_sites-1
         Dix=mask_tile(0,isite)
         Diy=mask_tile(1,isite)
         call calc_ix_iy(iCMx,iCMy,Dix,Diy,ithe,ix,iy)
         occ(ix,iy)=itile
      end do

      ln_cv=0
      do isite=0,N_sites-1
         Dix=mask_tile(0,isite)
         Diy=mask_tile(1,isite)
         call calc_ix_iy(iCMx,iCMy,Dix,Diy,ithe,ix,iy)
         s1=sp(ix,iy)
         o1=occ(ix,iy) 
            
         ixp=Mod(ix+1,Lx)
         s2=sp(ixp,iy)
         o2=occ(ixp,iy) 
         if(s1*s2.lt.0) then
            if(o1 .ne. -1) then
               if(o1.eq.o2) then
                  ln_cv=ln_cv+1
               end if
            end if
         end if
         
         iyp=Mod(iy+1,Ly)
         s2=sp(ix,iyp)
         o2=occ(ix,iyp)
         if(s1*s2.lt.0) then
            if(o1.eq.o2) then
               if(o1 .ne. -1) then
                  ln_cv=ln_cv+1
               end if
            end if
         end if
            
         if(Mod(iy,2) .eq. 0) ixp=Mod(ix-1+Lx,Lx)
         if(Mod(iy,2) .eq. 1) ixp=Mod(ix+1,Lx)
         s2=sp(ixp,iyp)
         o2=occ(ixp,iyp)
         if(s1*s2.lt.0) then
            if(o1.eq.o2) then
               if(o1 .ne. -1) then
                  ln_cv=ln_cv+1
               end if
            end if
         end if           
      end do
      length_covered=length_covered+ln_cv
         

c      do i1=0,N_sites_ev
c         isite=list_border(i1)
c         Dix=mask_tile(0,isite)
c         Diy=mask_tile(1,isite)
c         call calc_ix_iy(iCMx,iCMy,Dix,Diy,ithe,ix,iy)
c         occ(ix,iy)=itile
c      end do

      do i1=0,N_modules-1
         isite=list_modules(i1)
         Dix=mask_tile(0,isite)
         Diy=mask_tile(1,isite)
         call calc_ix_iy(iCMx,iCMy,Dix,Diy,ithe,ix,iy)
         anch(ix,iy)=type_modules(i1)
      end do

      end subroutine


      subroutine unoccupy(iCMx,iCMy,ithe,itile)
      implicit none
      include 'parameter.inc'
      include 'system.inc'
      integer iCMx,iCMy,itile,ithe
      integer Dix, Diy
      integer isite, ix, iy, i1, ln_cv
      integer ixp,iyp,o1,o2,s1,s2


      ln_cv=0
      do isite=0,N_sites-1
         Dix=mask_tile(0,isite)
         Diy=mask_tile(1,isite)
         call calc_ix_iy(iCMx,iCMy,Dix,Diy,ithe,ix,iy)
         s1=sp(ix,iy)
         o1=occ(ix,iy) 
            
         ixp=Mod(ix+1,Lx)
         s2=sp(ixp,iy)
         o2=occ(ixp,iy) 
         if(s1*s2.lt.0) then
            if(o1 .ne. -1) then
               if(o1.eq.o2) then
                  ln_cv=ln_cv+1
               end if
            end if
         end if
         
         iyp=Mod(iy+1,Ly)
         s2=sp(ix,iyp)
         o2=occ(ix,iyp)
         if(s1*s2.lt.0) then
            if(o1.eq.o2) then
               if(o1 .ne. -1) then
                  ln_cv=ln_cv+1
               end if
            end if
         end if
            
         if(Mod(iy,2) .eq. 0) ixp=Mod(ix-1+Lx,Lx)
         if(Mod(iy,2) .eq. 1) ixp=Mod(ix+1,Lx)
         s2=sp(ixp,iyp)
         o2=occ(ixp,iyp)
         if(s1*s2.lt.0) then
            if(o1.eq.o2) then
               if(o1 .ne. -1) then
                  ln_cv=ln_cv+1
               end if
            end if
         end if           
      end do
      length_covered=length_covered-ln_cv
      
      do isite=0,N_sites-1
         Dix=mask_tile(0,isite)
         Diy=mask_tile(1,isite)
         call calc_ix_iy(iCMx,iCMy,Dix,Diy,ithe,ix,iy)
         if(occ(ix,iy) .ne. itile) then
            write(*,*) "a-inconsistency in 'unoccupy'"
         end if
         occ(ix,iy)=-1
      end do

c      do i1=0,N_sites_ev
c         isite=list_border(i1)
c         Dix=mask_tile(0,isite)
c         Diy=mask_tile(1,isite)
c         call calc_ix_iy(iCMx,iCMy,Dix,Diy,ithe,ix,iy)
c         if(occ(ix,iy) .ne. itile) then
c            write(*,*) "a-inconsistency in 'unoccupy'"
c         end if
c         occ(ix,iy)=-1
c      end do
      
      do i1=0,N_modules-1
         isite=list_modules(i1)
         Dix=mask_tile(0,isite)
         Diy=mask_tile(1,isite)
         call calc_ix_iy(iCMx,iCMy,Dix,Diy,ithe,ix,iy)
         if(anch(ix,iy) .ne. type_modules(i1)) then
            write(*,*) "b-inconsistency in 'unoccupy'"
         end if
         anch(ix,iy)=0
      end do

      end subroutine

c     . this function provide the coordinates ix iy of the site specified by the tile's center of mass and orientation and Dix, Diy as provided by mask

      subroutine calc_ix_iy(iCM_x,iCM_y,Dix,Diy,the,ix,iy)

      implicit none
      include 'parameter.inc'
      include 'system.inc'
      integer ix, iy
      integer iCM_x,iCM_y,Dix,Diy,the
      integer Dix_rot, Diy_rot

c     theta=0 (Dix,Diy) => (Dix,Diy)
      if(the .eq. 0) then
         Dix_rot=Dix
         Diy_rot=Diy
      end if

c     theta=1 (Dix,Diy) => (0,Dix)+(-Diy,Diy)=(-Diy,Dix+Diy)
      if(the .eq. 1) then
         Dix_rot=-Diy
         Diy_rot=Dix+Diy
      end if

c     theta=2 (Dix,Diy) => (-Dix,Dix)+(-Diy,0)=(-Dix-Diy,Dix)
      if(the .eq. 2) then
         Dix_rot=-Dix-Diy
         Diy_rot=Dix
      end if

c     theta=3 (Dix,Diy) => (-Dix,-Diy)
      if(the .eq. 3) then
         Dix_rot=-Dix
         Diy_rot=-Diy
      end if

c     theta=4 (Dix,Diy) => (0,-Dix)+(Diy,-Diy)=(Diy,-Diy-Dix)
      if(the .eq. 4) then
         Dix_rot=Diy
         Diy_rot=-Diy-Dix
      end if

c     theta=5 (Dix,Diy) => (Dix,-Dix)+(Diy,0)=(Dix+Diy,-Dix)
      if(the .eq. 5) then
         Dix_rot=Dix+Diy
         Diy_rot=-Dix
      end if

      ix=iCM_x+Dix_rot
      iy=iCM_y+Diy_rot

      if(Diy_rot>0) then
         if(Mod(iCM_y,2).eq.0) then
            ix=ix+Diy_rot/2
         else
            ix=ix+(Diy_rot+1)/2
         end if
      end if

      if(Diy_rot<0) then
         if(Mod(iCM_y,2).eq.0) then
            ix=ix+(Diy_rot-1)/2
         else
            ix=ix+Diy_rot/2
         end if
      end if

c new BCs
      
      if(ix .lt. 0) ix=ix+Lx
      if(ix .ge. Lx) ix=ix-Lx
      
      if(fixed_boundary .eqv. .false.) then
         if(iy .lt. 0) iy=iy+Ly
         if(iy .ge. Ly) iy=iy-Ly
      end if

      end subroutine 

c.   if local we displace the tile on neighbouring sites
c     .   if local we translate or rotate with equal probability

      subroutine displace_tiles(local,acc)

      implicit none

      include 'parameter.inc'
      include 'system.inc'
      integer it
      integer ix_CM_old,iy_CM_old,ithe_old
      integer ix_CM_new,iy_CM_new,ithe_new
      integer Dx, Dy
      integer ix, iy
      integer i1, i2, i3, i4
      integer acc
      double precision ran2, ir, ir2, ir3
      double precision E_anc_o, E_anc_n
c      double precision E_sp_o, E_sp_n
c      double precision Wr_o, Wr_n
      double precision z_acc

      logical new

c     anc_n/o are the position of the anchors in the new/old position
c     vac_n/o are the holes to be filled. vac_n is a subset of anc_o
c     there are N_s vacancies
c      integer N_s
      integer anc_o(0:1,0:N_modules-1)
      integer anc_n(0:1,0:N_modules-1)
c      integer vac_n(0:1,0:N_modules-1)
c      integer vac_o(0:1,0:N_modules-1)
c      integer spin_o(0:N_modules-1)
c      integer spin_n(0:N_modules-1)

      logical overlap, ovl
      logical local

      acc=0

c     . choose a tile
      it=int(ran2(Iseed)*N_tiles)
      ix_CM_old=ixiy_tiles(0,it)
      iy_CM_old=ixiy_tiles(1,it)
      ithe_old=theta_tiles(it)

      if(local .eqv. .true.) then
         ir3=ran2(Iseed)
         if(ir3 .lt. 0.5) then
            ir=int(ran2(Iseed)*6)
            if(ir .eq. 0) then
               ix_CM_new= ix_CM_old+1
               iy_CM_new= iy_CM_old
            end if
            if(ir .eq. 1) then
               if(mod(iy_CM_old,2).eq.0) then
                  ix_CM_new= ix_CM_old
               else
                  ix_CM_new= ix_CM_old+1
               end if
               iy_CM_new= iy_CM_old+1
            end if
            if(ir .eq. 2) then
               if(mod(iy_CM_old,2).eq.0) then
                  ix_CM_new= ix_CM_old-1
               else
                  ix_CM_new= ix_CM_old
               end if
               iy_CM_new= iy_CM_old+1
            end if
            if(ir .eq. 3) then
               ix_CM_new= ix_CM_old-1
               iy_CM_new= iy_CM_old
            end if
            if(ir .eq. 4) then
               if(mod(iy_CM_old,2).eq.0) then
                  ix_CM_new= ix_CM_old-1
               else
                  ix_CM_new= ix_CM_old
               end if
               iy_CM_new= iy_CM_old-1
            end if
            if(ir .eq. 5) then
               if(mod(iy_CM_old,2).eq.0) then
                  ix_CM_new= ix_CM_old
               else
                  ix_CM_new= ix_CM_old+1
               end if
               iy_CM_new= iy_CM_old-1
            end if

            ix_CM_new=mod(ix_CM_new+Lx,Lx)
            iy_CM_new=mod(iy_CM_new+Ly,Ly)
            ithe_new=ithe_old
         else
            ir2=ran2(Iseed)
            if(ir2 .lt. 0.5) then
               ithe_new=mod(ithe_old-1+6,6)
            else
               ithe_new=mod(ithe_old+1,6)
            end if
            ix_CM_new=ix_CM_old
            iy_CM_new=iy_CM_old
         end if
      else

c. before, below was .true.
         if(fixed_boundary .eqv. .false.) then
            iy_CM_new=int(ran2(Iseed)*Ly)
         else
            iy_CM_new=int(ran2(Iseed)*(Ly-1))+1
         end if

         ix_CM_new=int(ran2(Iseed)*Lx)
         ithe_new=int(ran2(Iseed)*6)

      end if

c     . debug
c      iy_CM_new=iy_CM_old-4
c      ix_CM_new=ix_CM_old
c      ithe_new=ithe_old
c     . end debug

      ovl=overlap(ix_CM_new,iy_CM_new,
     &     ithe_new,it)

      if(ovl .eqv. .false.) then

c initialise anc_o anc_n
c we do not need this block of code
         call  init_list(anc_o,anc_n,
     &     ix_CM_old,iy_CM_old,ithe_old,
     &     ix_CM_new,iy_CM_new,ithe_new)

c     . debug ----------------
c         write(*,*) anc_o
c         write(*,*) anc_n
c         write(*,*) vac_o
c         write(*,*) vac_n
c         write(*,*) spin_o
c         write(*,*) N_s

c     . old configuration -----

c we do not need this block of code          
c. we set spin(vac_o) to 0
c         do i1=0,N_s-1
c            ix=vac_o(0,i1)
c            iy=vac_o(1,i1)
c            sp(ix,iy)=0
c         end do
c     . calculate the anchor energy
         call Eanchors(anc_o,E_anc_o)
c     . resample the old spins
c         new=.false. 
c         call sample_vacancies(vac_o,N_s,Wr_o,E_sp_o,spin_o,new)

cc     . new configuration -----
cc     . update spins and anchors
c         do i1=0,N_modules-1
c            ix=anc_n(0,i1)
c            iy=anc_n(1,i1) 
c            if(type_modules(i1).eq.2) sp(ix,iy)=2
c            if(type_modules(i1).eq.3) sp(ix,iy)=3
c         end do
c         do i1=0,N_s-1
c            ix=vac_n(0,i1)
c            iy=vac_n(1,i1)
c            sp(ix,iy)=0
c            spin_n(i1)=0
c         end do
c     . calculate the anchor energy         
         call Eanchors(anc_n,E_anc_n)
c     . sample the new spin and calculate Wn
c         new=.true.
c         call sample_vacancies(vac_n,N_s,Wr_n,E_sp_n,spin_n,new)

c     . acceptance
         z_acc=dexp(-E_anc_n+E_anc_o)
c         z_acc=z_acc*Wr_n/Wr_o

c. debug
c         z_acc=0.
c
         
         if(ran2(Iseed).lt.z_acc) then 
c. in this case, the spin configuration is up to date 
            call unoccupy(ix_CM_old,iy_CM_old,ithe_old,it)
            ixiy_tiles(0,it)=ix_CM_new
            ixiy_tiles(1,it)=iy_CM_new
            theta_tiles(it)=ithe_new
            call occupy(ix_CM_new,iy_CM_new,ithe_new,it)
c. update of the configurational energy 
            EnCk = EnCk + E_anc_n - E_anc_o
            acc=1
c         else
c     . rejection, we reverse the spin configuration
c     .    set anchors back
c            do i1=0,N_modules-1
c               ix=anc_o(0,i1)
c               iy=anc_o(1,i1) 
c               if(type_modules(i1).eq.2) sp(ix,iy)=2
c               if(type_modules(i1).eq.3) sp(ix,iy)=3
c            end do
c            do i1=0,N_s-1
c               ix=vac_o(0,i1)
c               iy=vac_o(1,i1)
c               sp(ix,iy)=spin_o(i1)
c            end do          
         end if
      end if

      end subroutine 

c ---------------------------------------------------------
c     vac_n (vac_o) is a subset of anc_o (anc_n)
c     vac_n is selected from anc_n and randomly permutated
      subroutine init_list(anc_o,anc_n,
     &     ix_CM_old,iy_CM_old,ithe_old,
     &     ix_CM_new,iy_CM_new,ithe_new)

      implicit none
      include 'parameter.inc'
      include 'system.inc'

      integer anc_o(0:1,0:N_modules-1)
      integer anc_n(0:1,0:N_modules-1)

      integer ix_CM_old,iy_CM_old,ithe_old
      integer ix_CM_new,iy_CM_new,ithe_new

      integer ix,iy, iDx, iDy
      integer i1,i2,i3,i4,i5


      do i1=0,N_modules-1
         
c. ... ... .preparing anc_o. ... .c
         iDx=mask_tile(0,list_modules(i1))
         iDy=mask_tile(1,list_modules(i1))
         call calc_ix_iy(ix_CM_old,iy_CM_old,
     &        iDx,iDy,ithe_old,ix,iy)
         anc_o(0,i1)=ix
         anc_o(1,i1)=iy

c. filling anc_n
         call calc_ix_iy(ix_CM_new,iy_CM_new,
     &        iDx,iDy,ithe_new,ix,iy)
         anc_n(0,i1)=ix
         anc_n(1,i1)=iy
         
      end do

      end subroutine 

c. ... ... ... ... .c
      
      subroutine Eanchors(anc,En)
      
      implicit none
      include 'parameter.inc'
      include 'system.inc'

      integer i1,i2,i3,i4
      integer s1, s2
      integer ix, iy
      integer ixm, iym, ixp, iyp, ix2
      integer anc(0:1,0:N_modules-1)
      double precision janch
      double precision En
      
      En=0
      do i1=0,N_modules-1

         ix=anc(0,i1)
         iy=anc(1,i1)
         
         if(type_modules(i1).eq.2) then
            janch=ja
         end if
         
         if(type_modules(i1).eq.3) then
            janch=jb
         end if
         
c. same site interaction
         s1=sp(ix,iy)
         En=En-janch*s1
         
         ixm=mod(ix-1+Lx,Lx)
         ixp=mod(ix+1,Lx)
         iym=mod(iy-1+Ly,Ly)
         iyp=mod(iy+1,Ly)
         ix2=ixp
         if(mod(iy,2).eq.0) ix2=ixm

         s2=sp(ixp,iy)
         En=En-janch*s2

         s2=sp(ixm,iy)
         En=En-janch*s2

         s2=sp(ix,iyp)
         En=En-janch*s2

         s2=sp(ix2,iyp)
         En=En-janch*s2

         s2=sp(ix,iym)
         En=En-janch*s2

         s2=sp(ix2,iym)
         En=En-janch*s2
         
      end do

      end subroutine 
      
c. ... ... ... ... .c



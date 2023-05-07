
c. the cluster iclb is merged into the cluster icla (if iclb>icla)

      subroutine mergeclusters(icla,iclb)

      implicit none
      include 'parameter.inc'
      include 'system.inc'
      integer icla, iclb
      integer iclm, iclp
      integer icl, ipar

      iclm=icla
      iclp=iclb
      if(iclb .lt. icla) then
        iclm=iclb
        iclp=icla
      end if

      NClust = NClust-1
c        write(*,*) icla,iclb
      NpIn(iclm)=NpIn(iclm)+NpIn(iclp)
      do icl=iclp,Nclust+1
        NpIn(icl)=NpIn(icl+1)
      end do
      do ipar=1,Np
        if(cl(ipar).eq.iclp) cl(ipar)=iclm
        if(cl(ipar).gt.iclp) cl(ipar)=cl(ipar)-1
      end do

      end

c.

      subroutine clusterstatistics(icycl)

      implicit none
      include 'parameter.inc'
      include 'system.inc'
      integer icla, iclb
      integer icl, ipar
      integer ix, iy, iz
      integer i, ingh
      integer icycl

      do ipar=1,Np
        cl(ipar) = 0
        NpIn(ipar) = 0
      end do
      NClust=0

      do ipar=1,Np
        cl(ipar)=NClust+1
        NpIn(NClust+1)=1
        NClust=NClust+1
        ix=x(ipar)
        iy=y(ipar)
        iz=z(ipar)
        do i=1,6
          call neigh(ingh,i,ix,iy,iz)
c         write(*,*) i, Nclust, ingh
          if(ingh.ne.0) then
            if((cl(ingh).ne.cl(ipar)).and.(cl(ingh).gt.0)) then
              call mergeclusters(cl(ipar),cl(ingh))
            end if
          end if
        end do
      end do

      write(101,*) icycl, Nclust
      do icl=1,NClust
        write(102,*) icl, NpIn(icl)
      end do

      end

c.

      subroutine neigh(ingh,i,ix,iy,iz)

      implicit none
      include 'parameter.inc'
      include 'system.inc'
      integer ingh, i , ix, iy, iz

      if(i.eq.1) then
        ingh=sp(mod(ix-1+Lx,Lx),iy,iz)
      end if
      if(i.eq.2) then
        ingh=sp(mod(ix+1,Lx),iy,iz)
      end if
      if(i.eq.3) then
        ingh=sp(ix,mod(iy-1+Ly,Ly),iz)
      end if
      if(i.eq.4) then
        ingh=sp(ix,mod(iy+1,Ly),iz)
      end if
      if(i.eq.5) then
        if(Mod(iy,2).eq.0) then
          ingh=sp(mod(ix-1+Lx,Lx),Mod(iy+1,Ly),iz)
        else
          ingh=sp(mod(ix+1,Lx),Mod(iy+1,Ly),iz)
        end if
      end if
      if(i.eq.6) then
        if(Mod(iy,2).eq.0) then
          ingh=sp(mod(ix-1+Lx,Lx),Mod(iy-1+Ly,Ly),iz)
        else
          ingh=sp(mod(ix+1,Lx),Mod(iy-1+Ly,Ly),iz)
        end if
      end if

      end



      

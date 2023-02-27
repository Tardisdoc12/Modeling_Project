      program gcmc_tbim_bulk 
c                                                                      
c       *********************************************************       
c       *                                                       *       
c       *       Semi-Grand Canonical Monte Carlo simulation     *
c       *                within the TBIM Model                  * 
c       *             for bulk metallic alloy                   * 
c       *                                                       * 
c       *                  Metropolis Algorithm                 * 
c       *********************************************************
c 
c                                                                       
c       ******************************************************          
c       *                                                    *          
c       *                  PARAMETERS                        *
c       *                                                    *
c       *   maximal number of atoms          : nmax          *
c       *   max. number of neighbours        : nvmax         *
c       *   time steps number                : npas          * 
c       *   periodicity for writing          : npw           *         
c       *   temperature  (Kelvin)            : temp          * 
c       *   composition                      : cmet          *
c       *   atoms number of metal 1          : nat1          *
c       *   chemical potential difference    : dmu           *
c       *     dmu= mu_B - mu_A if we change A in B           *
c       *                                                    * 
c       ******************************************************                
c 
c
      include 'tbim.common'
      real rran
      character*2 leg(10)
      character*12 filename(0:1000)
c
      pi=3.141592653589793d0                                            
c     boltzmann constant in eV/K
      cbol=8.62d-05
c
      rac2=dsqrt(2.d0)                                                  
      rac3=dsqrt(3.d0)
      rac5=dsqrt(5.d0)                                                  
      rac8=dsqrt(8.d0)
c
      open (11,file='gcmc_tbim_bulk.in',status='old')
      open (13,file='pr.xyz',status='unknown')
      open (14,file='concentration.dat',status='unknown')
      open (18,file='pr.in',status='unknown')
      open (19,file='pinit.in',status='unknown')
      open (21,file='energy.dat',status='unknown')
c
      read(11,*) npas,npw
      read(11,*) temp,dmu,idmumax,ddmu
      read(11,*) imax,jmax,kmax
      read(11,*) iconf 
      read(11,*) irand
      npeq=20*npw
c      npeq=npw
c
      print*,'starting random number:',irand
c
c     INIZIALIZZO RANMAR
c
      ij=int(rran(irand)*10000)
      kl=int(rran(irand)*10000)
      call rmarin(ij,kl)
c
      print*,'Nombre de macropas :',npas
      print*,' writing =',npw,' starting averages =',npeq
c
c     parametres du TBIM 
      call bimeta(leg)
c
c     le reseau est construit sur le parametre de maille arete=4.d0
      call reseau(imax,jmax,kmax,leg,iconf)
c
      print*,'Nombre de micropas = ',npas*natot
c
      print*,' '
      print*,'bulk of ',(leg(i),i=1,3)
      print*,' '
      print*,'Temperature = ',temp,'K'
c
c     -----------------------------------------------------------------
c     initialisation en energie
c
      call voisin
c
      call energy(ener)
      print*,'energie initiale =',ener
c
c     ----------------------------------------------------------------
c     loop on chemical potential 
c
      do 40 idmu=0,idmumax
      dmu=dmu+ddmu
c
      print*,' '
      print*,'Delta mu = ',dmu
      print*,'  ipas, <c2>, c2, <enrj>, <rej>'
c
      esum=0.d0
      c2sum=0.d0
      irej=0
c
      call subnomi(idmu,filename)
      open(42,file=filename(idmu),status='unknown')
c
c     -----------------------------------------------------------------
c     Monte Carlo loop
c
      do ipas=1,npas
c
          call mc_exchange(temp,dmu)
c
c     moyennes sur l energie,  etc...
      if(ipas.gt.npeq)  then
c
         call energy(ener)
         if (mod(ipas,npw).eq.0) then
            etot=0.d0
            do i=1,natot
               etot=etot+eneri(i)
            enddo
            drift=abs(ener-etot)
            if(drift.gt.0.001d0) print*,'drift in energy:',drift
         endif
         esum=esum+ener
         c2sum=c2sum+1.d0-dfloat(nat1)/dfloat(natot)
c
c     writing the outputs
c
         if(mod(ipas,npw).eq.0) then
c
         rewind(13)
         rewind(42)
         a=apar(1)
         write(13,*) natot
         write(13,*) bx,by,bz
         write(42,*) natot
         write(42,*) bx,by,bz 
         do i=1,natot
            write(13,500) leg(itype(i)),x(i)*a,y(i)*a,z(i)*a,itype(i)
            write(42,500) leg(itype(i)),x(i)*a,y(i)*a,z(i)*a,itype(i)
         enddo
         rej=dfloat(irej)/dfloat((ipas-npeq)*natot)*100.d0
         enrj=esum/(ipas-npeq)
         c2mean=c2sum/(ipas-npeq)
         cmet2=1.d0-dfloat(nat1)/dfloat(natot)
         print 200,dfloat(ipas)/npw,c2mean,cmet2,enrj,rej
c
         endif
      endif
c
      enddo 
c
      write(14,400) dmu,c2mean,enrj
      write(21,400) c2mean,enrj
c
      print*,' dmu, <c2>, <enrj> '
      print 400, dmu,c2mean,enrj
   40 continue
c
c     ################################################################
c
c  
  100 format(1x,2a2,4f10.5)    
  500 format(1x,a2,3f16.5,i4)
  200 format(1x,3f6.2,3f8.2)
  300 format(1x,5f8.4)
  400 format(1x,3f14.6,i5)
  600 format(1x,30f6.2)
c
      stop 
      end
c     -----------------------------------------------------------------
      subroutine bimeta(leg)
c
      include 'tbim.common'
c
      character*2 leg(10)
c
c ****************************************************************
c
c     leg(i)   : nature du metal i
c     apar(i) : parametre de maille du metal i (in angstrom)
c     ecoh(i) : energie de cohesion du metal i (in eV)
c
*******************************************************************
c
c         Cu
      leg(1)='Cu'
      apar(1)=3.62d0
      ecoh(1)=-3.50d0
c         Ag
      leg(2)='Ag'
      apar(2)=4.09d0  
      ecoh(2)=-2.95d0  
c         Ni
      leg(3)='Ni'
      apar(3)=3.52d0
      ecoh(3)=-4.44d0
c      
c
c ************************************************
c
c     systeme Cu-Ag  
      v1(4)=-0.0135d0
      tau(4)=(ecoh(2)-ecoh(1))/12.d0
      tc(4)=2086.d0    
c     systeme Ni-Cu
      v1(5)=-0.0015d0
      tau(5)=(ecoh(1)-ecoh(3))/12.d0
      tc(5)=487.d0
c     systeme Ni-Ag
      v1(6)=-0.1000d0
      tau(6)=(ecoh(2)-ecoh(3))/12.d0    
      tc(6)=4172.d0    
c
      print*,(leg(i),i=1,3)
      print*,' V1=',(v1(3+i),i=1,3)
      print*,' tau=',(tau(3+i),i=1,3)
c
c     pair interactions:
      do i=1,3
      epsi(i)=ecoh(i)/6.d0
      enddo
      do i=4,6
      epsi(i)=cbol*tc(i)/6.d0
      enddo 
c
      return
      end
c     ----------------------------------------------------------------- 
      subroutine energy(ener)
c
      include 'tbim.common'
c        
      ener=0.d0
c     cas d'un binaire: on décrit l'énergie en fonction de la
c     concentration en Ag (itype(i)=2) / TBIM
      do i=1,natot
         eneri(i)=0.d0
         if(itype(i).eq.2) then
           do j=1,nvois1(i)
              k=ivois1(j,i)
              if(itype(k).eq.2) eneri(i)=eneri(i)+v1(4)
           enddo
          eneri(i)=nvois1(i)*(tau(4)-v1(4))+eneri(i)
         endif
         ener=ener+eneri(i)
      enddo
c
c     Hamiltonien en interactions de paires pour le ternaire
c
      do i=1,natot
         eneri(i)=0.d0
         do j=1,nvois1(i)
            k=ivois1(j,i)
            if((itype(i).eq.1).and.(itype(k).eq.1)) then
               itypik=1
            else if((itype(i).eq.2).and.(itype(k).eq.2)) then
               itypik=2
            else if((itype(i).eq.3).and.(itype(k).eq.3)) then
               itypik=3
            else if(((itype(i).eq.1).and.(itype(k).eq.2)).or.
     ^             ((itype(i).eq.2).and.(itype(k).eq.1))) then
               itypik=4
            else if(((itype(i).eq.1).and.(itype(k).eq.3)).or.
     ^            ((itype(i).eq.3).and.(itype(k).eq.1))) then
               itypik=5
            else
               itypik=6
            endif
            eneri(i)=eneri(i)+epsi(itypik)
         enddo
         ener=ener+eneri(i)
      enddo
c        
      return
      end
c     ----------------------------------------------------------------- 
      subroutine modifener(deipc,ipc1)
c
      include 'tbim.common'
c        
      deipc=0.d0
      i=ipc1
      eneri(i)=0.d0
      if(itype(i).eq.2) then
        do j=1,nvois1(i)
           k=ivois1(j,i)
           if(itype(k).eq.2) eneri(i)=eneri(i)+v1(4)
        enddo
        eneri(i)=nvois1(i)*(tau(4)-v1(4))+eneri(i)
      endif
      deipc=deipc+eneri(i)-enerio(i)
c
      do n=1,nvois1(ipc1)
         i=ivois1(n,ipc1)
         eneri(i)=0.d0
         if(itype(i).eq.2) then
           do j=1,nvois1(i)
              k=ivois1(j,i)
              if(itype(k).eq.2) eneri(i)=eneri(i)+v1(4)
           enddo
         eneri(i)=nvois1(i)*(tau(4)-v1(4))+eneri(i)
         endif
         deipc=deipc+eneri(i)-enerio(i)
      enddo
c
      return
      end
c     ------------------------------------------------------------
      subroutine reseau(imax,jmax,kmax,leg,iconf)
c
c     construction d un slab fcc(100)
c
      include 'tbim.common'
c
      character*2 symb,leg(10)
c
      nats=(imax+1)*(jmax+1)*2
      natot=(imax+1)*(jmax+1)*(kmax+1)*4
c
      if (natot.gt.nmax0) then
         write (*,505)
         stop
      endif
c
c     ------------------------------
c     construction du slab (100) pur
c     ------------------------------
      if(iconf.eq.0) then
        l=0
c
        do k=0,2*kmax+1
        do j=0,jmax
        do i=0,imax
           l=l+1
           pk=dfloat(mod(k,2))/2.d0
           x(l)=dfloat(i)+pk
           y(l)=dfloat(j)
           z(l)=dfloat(k)/2.d0
           x(l+1)=dfloat(i)+0.5d0-pk
           y(l+1)=dfloat(j)+0.5d0
           z(l+1)=dfloat(k)/2.0d0
           l=l+1
        enddo
        enddo
        enddo
c
        if(l.ne.natot) then
          write (*,510) l,natot
          stop
        endif
c
c     dimension de la boite : bx,by,bz
c
        bx=dfloat(imax+1)
        by=dfloat(jmax+1)
        bz=dfloat(kmax+1)
c
        do i=1,natot
           itype(i)=1
        enddo
      else
c        _________________________________
c        Lecture de configuration initiale
c        _________________________________
         print*,'Lecture de configuration initiale'
         open(2,file='p.in',status='old')
         read(2,*) natot
         read(2,*) bx,by,bz 
         do i=1,natot
            read(2,500) symb,x(i),y(i),z(i)
            if(symb.eq.leg(1)) itype(i)=1
            if(symb.eq.leg(2)) itype(i)=2
         enddo
      endif
c
      nat1=0
      nat2=0
c
      do i=1,natot
         if(itype(i).eq.1) nat1=nat1+1
         if(itype(i).eq.2) nat2=nat2+1
      enddo
      print*,nat1,nat2
      cmet2=dfloat(nat2)/natot
      print*,cmet2
c
      if((nat1+nat2).ne.natot) then
        print*,'Pb de composition chimique'
        stop 12
      endif
c
c     writing the initial positions:
         a=1.d0
         if(iconf.eq.0) a=apar(2)
         rewind 18
         write(18,*) natot
         write(18,*) bx,by,bz
        rewind 19
         write(19,*) natot
         write(19,*) bx,by,bz         
         do i=1,natot
c            itype(i)=1
            write(18,500) leg(itype(i)),x(i)*a,y(i)*a,z(i)*a,itype(i)
            write(19,500) leg(itype(i)),x(i),y(i),z(i)
         enddo
c
      print*,' '
      print*,natot,' atoms : ',nat1,leg(1),' and',nat2,leg(2)
      print*,' '
      print*,'concentration en ',leg(2),': ',cmet2
      print*,' '
c
  500 format(1x,a2,3f16.5,i4)
  505 format(1x,' attention, trop d atomes : natot > nmax')
  510 format(2i5,1x,' attention, le cristal est mal construit ')
      return
      end
c
c     ------------------------------------------------------------------
      subroutine voisin 
c                                                     
      include 'tbim.common'
c
      dcut1=1.d0/rac2+0.01d0
      dcut2=1.d0+0.01d0
      dcut3=rac3/rac2+0.01d0
       print*,'dcut',dcut1,dcut2,dcut3
c
      do i=1,natot
         nvois1(i)=0
      enddo
c
      do i=1,natot-1  
         do j=i+1,natot
c
            xij=x(j)-x(i)  
            yij=y(j)-y(i)  
            zij=z(j)-z(i)
c
            if (dabs(xij+bx).lt.dabs(xij)) xij=xij+bx
            if (dabs(xij-bx).lt.dabs(xij)) xij=xij-bx
            if (dabs(yij+by).lt.dabs(yij)) yij=yij+by
            if (dabs(yij-by).lt.dabs(yij)) yij=yij-by
            if (dabs(zij+bz).lt.dabs(zij)) zij=zij+bz
            if (dabs(zij-bz).lt.dabs(zij)) zij=zij-bz
c
            dij=dsqrt(xij*xij+yij*yij+zij*zij)
c
            if (dij.lt.dcut1) then
               nvois1(i)=nvois1(i)+1
               nvois1(j)=nvois1(j)+1
               ivois1(nvois1(i),i)=j
               ivois1(nvois1(j),j)=i
            endif
c
          enddo
          if(i.eq.1) print*,i,nvois1(i)
      enddo
c
      return                                                            
      end     
c     ________________________________________________________________
      subroutine mc_exchange(temp,dmu)
c
      include 'tbim.common'
c
      do iex=1,natot
c
c     tirage au hazard de l atome a transmuter
c
      call ranmar(aaa)
      ipc1 = int(natot * aaa + 1)
c
c        on conserve les donnees avant echange 
         do i=1,natot
            enerio(i)=eneri(i)
         enddo
c
c        echange
c
         if(ibin.eq.4) then
           if(itype(ipc).eq.3) goto 500
           it12=itype(ipc)
           itype(ipc)=mod(it12,2)+1
           dn1=it12-itype(ipc)
           nat(1)=nat(1)+dn1
           nat(2)=nat(2)-dn1
         endif
         if(ibin.eq.5) then
           if(itype(ipc).eq.2) goto 500
           it13=itype(ipc)
           if(it13.eq.1) itype(ipc)=3
           if(it13.eq.3) itype(ipc)=1
           dn1=(it13-itype(ipc))/2.d0
           nat(1)=nat(1)+dn1
           nat(3)=nat(3)-dn1
         endif
         if(ibin.eq.6) then
           if(itype(ipc).eq.1) goto 500
           it23=itype(ipc)
           itype(ipc)=mod(it23,3)+2
           dn1=it23-itype(ipc)
           nat(2)=nat(2)+dn1
           nat(3)=nat(3)-dn1
         endif
c
         call modifener(deipc,ipc1)
c         call energy(ener)
c
c         print*,ener,enero,dn1 
c
c         de = ener - enero + dn1 * dmu
         de=deipc + dn1 * dmu
         if(de.le.0.d0)  goto 500
         ex = exp(-de/(cbol*temp))
         call ranmar(aaa)
         if(aaa.lt.ex) goto 500
c
c        la nouvelle configuration est refusee, retour a l ancienne
c
         if(ipas.gt.npeq) irej=irej+1
         do i=1,natot
            eneri(i)=enerio(i)
         enddo
        if(ibin.eq.4) then
           itype(ipc)=it12
           nat(1)=nat(1)-dn1
           nat(2)=nat(2)+dn1
         endif
         if(ibin.eq.5) then
           itype(ipc)=it13
           nat(1)=nat(1)-dn1
           nat(3)=nat(3)+dn1
         endif
         if(ibin.eq.6) then
           itype(ipc)=it23
           nat(2)=nat(2)-dn1
           nat(3)=nat(3)+dn1
         endif
c
c         print*,'nat1',nat1,de,'refus'
c
  500 continue
c      print*,'nat1',nat1,de 
      enddo
      return
      end
c     -----------------------------------------------------------------
      FUNCTION RRAN(ISEED)
      PARAMETER(IA=7141,IC=54773,IM=259200)
      SAVE
      ISEED=MOD(ISEED*IA+IC,IM)
      RRAN=FLOAT(ISEED)/FLOAT(IM)
      RETURN
      END
c     ---------------------------------------------------------------
      SUBROUTINE RMARIN(IJ,KL)
      implicit integer (i-n),double precision (a-h,o-z)
      COMMON/RASET1/U(97),C,CD,CM,I97,J97
      I=MOD(IJ/177,177)+2
      J=MOD(IJ,177)+2
      K=MOD(KL/169,178)+1
      L=MOD(KL,169)
      DO II=1,97
      X=0.
      Y=.5
      DO JJ=1,24
      M=MOD(MOD(I*J,179)*K,179)
      I=J
      J=K
      K=M
      L=MOD(53*L+1,169)
      IF (MOD(L*M,64).GE.32) X=X+Y
      ENDDO
      Y=.5*Y
      ENDDO
      U(II)=X
      C=362436./16777216.
      CD=7654321./16777216.
      CM=16777213./16777216.
      I97=97
      J97=33
      RETURN
      END
!
      SUBROUTINE RANMAR(RVEC)
      implicit integer (i-n),double precision (a-h,o-z)
      double precision rvec
      COMMON/RASET1/U(97),C,CD,CM,I97,J97
!
      UNI=U(I97)-U(J97)
      IF (UNI.LT.0.) UNI=UNI+1.
      U(I97)=UNI
      I97=I97-1
      IF (I97.EQ.0) I97=97
      J97=J97-1
      IF (J97.EQ.0) J97=97
      C=C-CD
      IF (C.LT.0.) C=C+CM
      UNI=UNI-C
      IF (UNI.LT.0.) UNI=UNI+1.
      RVEC=UNI
!
      RETURN
      END
c     --------------------------------------------------------------
      subroutine subnomi(mnomi,filename)
c
c     probleme quand on part de 190 atomes...
c
      character*12 filename(0:1000)
c
c     mnomi est le numero de l atome depose
c      do 10 i=1,mnomi
c
        i=mnomi
        if(i.lt.10) then
        filename(i)='p0.xyz'
        write(filename(i)(2:2),'(i1)') i
         elseif((i.ge.10).and.(i.lt.100)) then
        filename(i)='p00.xyz'
        write(filename(i)(2:3),'(i2)') i
         elseif((i.lt.1000).and.(i.ge.100)) then
         filename(i)='p000.xyz'
        write(filename(i)(2:4),'(i3)') i
         elseif((i.lt.10000).and.(i.ge.1000)) then
         filename(i)='p0000.xyz'
        write(filename(i)(2:5),'(i3)') i
         endif
c       enddo

        return
        end

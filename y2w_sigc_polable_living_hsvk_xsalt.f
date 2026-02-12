      program platem
      implicit double precision (a-h,o-z)
      include 'PBlive.inc'
      dimension cchat(0:maxel),achat(0:maxel),
     *ttfdmon(0:maxel),ttfdnmon(0:maxel),fp(0:maxel)
c     *ttcchat(0:maxel),ttachat(0:maxel)
c     ALTERNATING CHARGES!
      do i = 1,maxel
      cchat(i) = 1.d0
      achat(i) = 1.d0      
      enddo
      ifc = 38
      ife = 36      
      ifnc = 34      
      ins = 49
      issolv = 56
      isd = 58
      idpot = 63
      ksurf = 65
      ikh = 67      
      pi = acos(-1.d0)
      onethird = 1.d0/3.d0
      twothirds = 2.d0/3.d0
      bk = 1.38066D-23
      avno = 6.02214D23
      epszero = 8.854D-12
      elch = 1.602D-19
      T = 298.d0
      epsilon = 78.3
      fourpi = 4.d0*pi
      twopi = 2.d0*pi
      fourpi = 4.d0*pi
      rtwelve = 1.d0/12.d0
      rthree = 1.d0/3.d0
      volfact = fourpi/3.d0
      rnine = 1.d0/9.d0
      rthree = 1.d0/3.d0
      rphi = fourpi*0.2d0 
      aphi = fourpi*0.5d0
      pis = pi/6.d0
      pit = pi/3.d0
      pif = pi/4.d0     
      es22 = -32.d0*pi/9.d0
      a1 = 1.d0
      a2 = 2.45696d0
      b1 = 1.d0
      b2 = 4.10386d0
      c1 = -1.d0
      c2 = -3.75503
      AA1 = 2.d0*c1-2.d0*a1-4.d0
      AA2 = 2.d0*c2-2.d0*a2-4.d0
      BB1 = 3.d0-b1+a1-3.d0*c1
      BB2 = 3.d0-b2+a2-3.d0*c2
      ddtol = 0.0000001d0
      delecok = 1.D-11
      smet = 1.d-10
      unorm = elch*elch/(4.d0*pi*epszero*epsilon*smet)
      bkT = bk*T
      T = bkT/unorm
      rrT = 1.d0/T
      bjerrum = rrT
      rbjerrum = rrT
      write(*,*) 'unorm,bkT /J = ',unorm,bkT
      write(*,*) 'T,bjerrum = ',T,bjerrum
c     CLOSE TO THE WALLS, THE DENSITY IS ASSUMED TO BE ZERO
      open (ifc,file='fcdfil',form='formatted')
      open (ifnc,file='fncdfil',form='formatted')
      open (ife,file='fesfil',form='formatted')               
      open (ins,file='dpcinp',form='formatted')
      open (idpot,file='donnan',form='formatted')
      open (ksurf,file='surfdens',form='formatted')      
      open (ikh,file='separation',form='formatted')                  
      rewind idpot      
      rewind ifc
      rewind ife
      rewind ifnc      
      rewind ins
      rewind ksurf
      rewind ikh
      read(ins,*) 
      read(ins,*) bdm,bdes
      read(ins,*) 
      read(ins,*) nval,nsval
      nval = 1
      nsval = 1
      read(ins,*) 
      read(ins,*) rkpp,rkpn      
c      if (rkpn.gt.rkpp) then
c      write(*,*) 'rkpp must be gt. rkpn'
c      stop
c      endif
      rlambda = (rkpp+rkpn)/2.d0
      rDelta = (rkpp-rkpn)/2.d0
      write(*,*) 'rkpp,rkpn = ',rkpp,rkpn      
      write(*,*) 'rlambda,rDelta = ',rlambda,rDelta
      write(*,*)
      write(*,*) 
c     bulk values:
      ekpn = dexp(-rkpn)
      ekpp = dexp(-rkpp)

      rK = 1.d0-ekpp-ekpn
      write(*,*) 'rK = ',rK

c     avr = ekpn/(1.d0-ekpp-ekpn)+1.d0
      avr = 1.d0/rK      
      write(*,*) 'avr = ',avr

      bdpol = 2.d0*bdm/avr
      write(*,*) 'bdpol = ',bdpol
      ronedens = bdpol*rK
      write(*,*) 'ronedens = ',ronedens
      
      abulkchat = dsqrt(bdm/(0.5d0*rK*bdpol))
c      bulkchat = 1.d0/(1.d0-ekpn-ekpp)
      write(*,*) 'abulkchat = ',abulkchat      
c      write(*,*) 'bulkchat = ',bulkchat

c      bahat = bulkchat
c      bchat = bulkchat
c      chat = 1.d0+ekpn*bahat+ekpp*bchat
c      ahat = 1.d0+ekpn*bchat+ekpp*bahat      
c      write(*,*) 'ahat,chat = ',ahat,chat
c      bc1sum = bulkchat
c      bc2sum = bulkchat
c      write(*,*) 'bc1sum,bc2sum = ',bc1sum,bc2sum      
c      bcch = (1.d0+ekpn*bc1sum+ekpp*bc2sum)
c      bach = (1.d0+ekpp*bc1sum+ekpn*bc2sum)      
c      write(*,*) 'bcch,bach = ',bcch,bach

c      cavr = (1.-ekpn)*(1.-ekpp)/(1.-ekpn-ekpp)**2
c      write(*,*) 'cavr = ',cavr,rK/(1.-ekpn-ekpp)**2
c      write(*,*) '0.5*bdpol*cavr,bdm = ',0.5*bdpol*cavr,bdm
c      write(*,*) '0.5d0*rK*bdpol*bchat**2 = ',0.5d0*rK*bdpol*bchat**2
c      write(*,*) '0.5*rK*bdpol*bulkchat**2 = ',0.5*rK*bdpol*bulkchat**2
      write(*,*) '0.5*rK*bdpol*abulkchat**2=',0.5*rK*bdpol*abulkchat**2             

c      stop
      
      read(idpot,*) donn
      frdonn = 0.001d0
      read(ksurf,*) surfdens      
      read(ins,*) 
      read(ins,*) htrams
      read(ikh,*) h      
      read(ins,*)
      read(ins,*) dz
      read(ins,*)      
      read(ins,*) dmm
      read(ins,*)      
      read(ins,*) dms
      read(ins,*)      
      read(ins,*) kread
      read(ins,*)      
      read(ins,*) ioimaxm
      read(ins,*)
      read(ins,*) trams
      read(ins,*)      
c     read(ins,*) dm1,ds1,bl
      read(ins,*) dm1,dmw,bl,vk
      ds1 = dm1      
      write(*,*) 'dm1,dmw = ',dm1,dmw
      dhs = dm1
      write(*,*) 'vk,bl = ',vk,bl
      rbl = 1.d0/bl
      bl2 = bl*bl
      bl3 = bl2*bl
      rbl3 = 1.d0/bl3
      dhs2 = dhs*dhs
      dhs3 = dhs2*dhs
      rdhs3 = 1.d0/dhs3         

      rnval = real(nval)
      rrnval = 1.d0/rnval
      rnsval = real(nsval)
      rrnsval = 1.d0/rnsval
      halfh = 0.5d0*h
      nhalfh = int(halfh/dz+0.01)
      cbdm = 0.d0
      bdnm = bdm
      bds = bdes
      write(*,*) 'surdfdens = ',surfdens 
      write(*,*) 'bdm,bdes = ',bdm,bdes
      write(*,*) 'bdnm,bds = ',bdnm,bds
      tdmm = 1.d0-dmm
      tdms = 1.d0-dms
c      rrT = 1.d0/T
      rdz = 1.d0/dz
c     twopidz = twopi*dz
      twopidz = 0.5d0*dz*rbl      
      irdz = int(rdz+0.001d0)
      ibl = int(bl/dz+0.01d0)      

      closew = 0.5d0*dmw
      write(*,*) 'closew = 0.5*dmw = ',closew
      nfack = int((h+dmw)/dz+0.01d0)
      imitt = int((h+dmw)/dz+0.01d0)/2                  
      sclosew = closew      
      istart = int(closew/dz+0.01d0)
      iblnw = ibl+int(closew/dz+0.01d0)
      csurf = -2.d0*pi*surfdens*rrT*(h+dmw)                
      write(*,*) 'csurf = ',csurf            

      nfin = 2*imitt
      istp1 = istart+1 
      islut = int((h+0.5d0*dmw)/dz+0.01d0)
      istp1s = istp1
      isluts = islut
      inw = istp1+ibl-1      
      write(*,*) 'h = ',h
      write(*,*) 'closew,sclosew = ',closew,sclosew
      write(*,*) 'ism,ibl = ',ism,ibl
      write(*,*) 'nfack,nfin = ',nfack,nfin
      write(*,*) 'istp1,islut = ',istp1,islut
      write(*,*) 'imitt = ',imitt

      dm2 = dm1*dm1
      dm3 = dm2*dm1
      rdm3 = 1.d0/dm3
      ds2 = ds1*ds1
      ds3 = ds2*ds1
      rds3 = 1.d0/ds3      
      phzmm = -(6.d0*pi*dm2/5.d0)
      rphi = rphi*dm2
      aphi = aphi*dm2
      phzss = phzmm
      phzms = phzmm
      
      es22ss = 0.d0
      es22mm = es22ss
      es22ms = es22ss
      Pblj = 0.5d0*es22mm*bdm**2

      q1 = ds1
      p1 = dm1
      q2 = q1*q1
      q3 = q2*q1
      p2 = p1*p1
      p3 = p2*p1
      dm2 = dm1*dm1
      dm3 = dm2*dm1
      rp3 = 1.d0/p3
      rq3 = 1.d0/q3
      rdm3 = 1.d0/dm3
      write(*,*) 'dm1,q1,p1 = ',dm1,q1,p1            
      ism = int(dm1/dz+0.01d0)
      isms = int(ds1/dz+0.01d0)
      ismms = int(p1/dz+0.01d0)
      write(*,*) 'ism,isms,ismms = ',ism,isms,ismms

      bdt = 2.d0*(vk*bdm+bdes)*dhs3
      xsib = 1.d0-pis*bdt
      rxsib = 1.d0/xsib
      rxsibsq = rxsib*rxsib
      aex1= -(c1+1.d0)*dlog(xsib)-
     *0.5d0*(AA1*pis*bdt+BB1*(pis*bdt)**2)*rxsibsq
      Vdae1dV = -pis*bdt*rxsib*(c1+1.d0-
     *0.5d0*(AA1+2.d0*BB1*pis*bdt)*rxsib-
     *pis*bdt*rxsibsq*(AA1+BB1*pis*bdt))
      rNpda1 = pis*bdm*rxsib*((c1+1.d0)-
     *0.5d0*AA1*rxsib*(1.d0+2.d0*pis*bdt*rxsib)-
     *BB1*pis*rxsib*bdt*(1.d0+pis*bdt*rxsib))
      rNsda1 = pis*bds*rxsib*((c1+1.d0)-
     *0.5d0*AA1*rxsib*(1.d0+2.d0*pis*bdt*rxsib)-
     *BB1*pis*rxsib*bdt*(1.d0+pis*bdt*rxsib))      
      baex1 = aex1
      exchpp = vk*(aex1+2.d0*rNpda1*dhs3*vk+2.d0*rNsda1*dhs3)
      exchs = aex1+2.d0*rNpda1*dhs3*vk+2.d0*rNsda1*dhs3
      exPp = -2.d0*bdm*Vdae1dV*vk
      exPs = -2.d0*bdes*Vdae1dV      
      exP = exPp+exPs
c     Pbhs = 2.d0*bdpol+exPp
c     Pbhs = exPp
      Pbhs = bdpol+bds+bdes+exP     
      Pb = Pbhs
      scalem = exchpp/2.d0
      emscale = 2.d0*scalem
      chemps = dlog(bds)+exchs
      chempes = dlog(bdes)+exchs      
      scales = chemps
      scalees = chempes            
      write(*,*) 'bdm,bdes = ',bdm,bdes
      write(*,*) 'max no. of iterations = ',ioimaxm
      write(*,*) 'bjerrum,rbjerrum = ',bjerrum, rbjerrum
      write(*,*) 'h,dz = ',h,dz
      write(*,*) 'exchpp = ',exchpp
      write(*,*) 'coion chemical potential = ',chemps
      write(*,*) 'counterion chemical potential = ',chempes            
      write(*,*) 'total bulk pressure = ',Pb      
      write(*,*) 'rrT = ',rrT
      write(*,*) 'nval,nsval = ',nval,nsval
      write(*,*) 'surfdens = ',surfdens
      write(*,*) 'istp1,islut,inw = ',istp1,islut,inw
      write(*,*) 'ism,ibl = ',ism,ibl
      write(*,*) 'dmm,dms (density mixing param. mon.,solv.) = ',dmm,dms
      write(*,*) 'bond length (bl): ',bl
      write(*,*) 'monomer AND solvent hs diameter (dhs): ',dhs
      write(*,*) 'NOTE: MONOMERS/SOLVENT PARTICLES HAVE THE SAME SIZE!'       

      z = 0.5d0*dz
      iz = 1
      zp = z-dz
      do jz = 1,maxel
      zp = zp+dz
      diffz = abs(z-zp)
      iii = iabs(jz-iz)
      phiz = -2.d0*pi*diffz*rrT*dz
      Phimm(iii) = phiz
      enddo      
      if (kread.eq.0) then
      donn = 0.d0
      do iz = istp1,imitt       
      fdmon(iz) = bdm
      fdnmon(iz) = bdm
      cchat(iz) = dsqrt(bdm/(0.5d0*rK*bdpol))
      achat(iz) = dsqrt(bdm/(0.5d0*rK*bdpol))
      fdes(iz) = bdes
      fdsol(iz) = bds            
      enddo
      else
      do iz = 1,imitt                
      read(ifc,*) trams,fdmon(iz),cchat(iz)
      read(ifnc,*) trams,fdnmon(iz),achat(iz)
      read(ife,*) trams,fdsol(iz),fdes(iz)                 
      enddo
      endif

      sum = 0.d0
      do i = istp1,imitt
      sum = sum+fdmon(i)-fdnmon(i)-fdsol(i)+fdes(i)
      enddo
      sum = sum*dz
      elec = surfdens+sum
      write(*,*) 'surfdens,elec (st) = ',surfdens,elec          

      donnB = donn
      
      elecA = elec      
      write(*,*) 'elecA,donn = ',elecA,donn
      if (dabs(elecA).lt.delecok) goto 985
      donnA = donn
      if (dabs(donn).lt.0.001d0) then
      ddonn = 0.05d0
      else
      ddonn = frdonn*donn
      endif
      if (elecA.lt.0.d0) then
      donnB = donnA+ddonn
      else
      donnB = donnA-ddonn
      endif
      ddd = donnB
 4949 continue
      donnB = ddd
      write(*,*) 'donnB,ddd = ',donnB,ddd
      sumfdm = 0.d0
      sumfdnm = 0.d0      
      sumfds = 0.d0
      sumfdes = 0.d0
      do i = istp1s,imitt
      sumfdm = fdmon(i)*dexp(-(donnB-donn))+sumfdm
      sumfdnm = fdnmon(i)*dexp((donnB-donn))+sumfdnm      
      sumfdes = fdes(i)*dexp(-(donnB-donn))+sumfdes      
      sumfds = fdsol(i)*dexp((donnB-donn))+sumfds
      enddo
      elecB = (sumfdm+sumfdes-sumfds-sumfdnm)*dz+surfdens
      ddd = donnA-elecA*(donnB-donnA)/(elecB-elecA)      
      write(*,*) donn,donnB,elecB
      if (dabs(elecB).lt.delecok) goto 985
c      write(*,*) 'elecA,elecB = ',elecA,elecB
c      write(*,*) donnA,donnB,ddd
c      write(*,*) 'donnA,donnB = ',donnA,donnB,ddonn
c      write(*,*) (donnB-donnA)/(elecB-elecA),ddd
c      write(*,*) 
c      if (iop.gt.12) stop      
      donnA = donnB
      elecA = elecB
      donnB = ddd
      goto 4949
 985  continue
      write(*,*) 
      write(*,*) 'donn,donnB,elecB = ',donn,donnB,elecB

      sum = 0.d0
      do i = istp1,imitt
      fdmon(i) = fdmon(i)*dexp(-(donnB-donn))
      fdnmon(i) = fdnmon(i)*dexp((donnB-donn))
      fdes(i) = fdes(i)*dexp(-(donnB-donn)) 
      fdsol(i) = fdsol(i)*dexp((donnB-donn))      
      sum = sum+fdmon(i)-fdnmon(i)-fdsol(i)+fdes(i)      
      enddo
      celec = surfdens+sum*dz
      write(*,*) 'surfdens,celec (st) = ',surfdens,celec
      donn = donnB
      
      do i = 1,istp1-1
      fdmon(i) = 0.d0
      fdnmon(i) = 0.d0
      cchat(i) = 0.d0
      achat(i) = 0.d0
      ehbclam(i) = 0.d0
      ehbclanm(i) = 0.d0
      fdsol(i) = 0.d0
      fdes(i) = 0.d0      
      enddo
      jz = imitt+1
      do iz = imitt+1,nfin      
      jz = jz-1
      fdmon(iz) = fdmon(jz)      
      fdnmon(iz) = fdnmon(jz)
      fdsol(iz) = fdsol(jz)      
      fdes(iz) = fdes(jz)      
      cchat(iz) = cchat(jz)
      achat(iz) = achat(jz)
      enddo
      
      ddmax = -10000.
      niter = 0
 100  continue
      niter = niter+1
c      write(*,*) niter
      if (mod(niter,100).eq.0) then
      sum = 0.d0
      do i = istp1,imitt         
      sum = sum+fdmon(i)-fdnmon(i)-fdsol(i)+fdes(i)      
      enddo
      celec = surfdens+sum*dz      
      write(*,*) 'niter,ddmax,celec =',niter,ddmax,celec
      endif
      if (niter.gt.ioimaxm) goto 200      
      CALL CALCCD
      CALL AVEC
      CALL CALCELAM
      ddmax = -10000.
      
      sumfdm = 0.d0
      sumfdnm = 0.d0
      sumfds = 0.d0
      sumfdes = 0.d0
      do iz = istp1,imitt      
      kz = iz-ibl
      factkz = 0.5d0
      if (kz.lt.istp1) then
      kz = istp1
      factkz = 1.d0
      endif
      c1sum = ehbclanm(kz)*achat(kz)*factkz
      c2sum = ehbclam(kz)*cchat(kz)*factkz
c      a1sum = ehbclam(kz)*cchat(kz)*factkz
c      a2sum = ehbclanm(kz)*achat(kz)*factkz
c     c1sum = a2sum
c     c2sum = a1sum
      do jz = kz+1,iz+ibl-1
      c1sum = ehbclanm(jz)*achat(jz)+c1sum
      c2sum = ehbclam(jz)*cchat(jz)+c2sum
c      a1sum = ehbclam(jz)*cchat(jz)+a1sum
c      a2sum = ehbclanm(jz)*achat(jz)+a2sum            
      enddo
      jz = iz+ibl
      c1sum = (0.5d0*ehbclanm(jz)*achat(jz)+c1sum)*twopidz
      c2sum = (0.5d0*ehbclam(jz)*cchat(jz)+c2sum)*twopidz
c      a1sum = (0.5d0*ehbclam(jz)*cchat(jz)+a1sum)*twopidz
c      a2sum = (0.5d0*ehbclanm(jz)*achat(jz)+a2sum)*twopidz
c      a1sum = c2sum
c      a2sum = c1sum
      tchat = ehbclam(iz)*(1.d0+ekpn*c1sum+ekpp*c2sum)
      tahat = ehbclanm(iz)*(1.d0+ekpn*c2sum+ekpp*c1sum)

      ttfdmon(iz) = 0.5d0*rK*bdpol*tchat**2
      ttfdnmon(iz) = 0.5d0*rK*bdpol*tahat**2
      
      sumfdm = ttfdmon(iz)+sumfdm
      sumfdnm = ttfdnmon(iz)+sumfdnm
cc      fp(iz) = 0.5d0*rK*bdpol*ekappa*tcchat*ehbclam(iz)
cc      fp(iz) = fp(iz)+0.5d0*rK*bdpol*ekappa*tachat*ehbclanm(iz)
      fp(iz) = 0.5d0*rK*bdpol*tchat*ehbclam(iz)
      fp(iz) = fp(iz)+0.5d0*rK*bdpol*tahat*ehbclanm(iz)      
c      ttcchat(iz) = tcchat
c     ttachat(iz) = tachat
      sumfds = eblam(iz)+sumfds
      sumfdes = eeblam(iz)+sumfdes      
      enddo
      stfdm = sumfdm*dz
      stfdnm = sumfdnm*dz      
      stfdes = sumfdes*dz      
      stfds = sumfds*dz
      elecA = stfdm+stfdes-stfds-stfdnm+surfdens

c      write(*,*) sumfdm,sumfds,ttfdmon(imitt)
      
      ddd = donn
      if (dabs(elecA).lt.delecok) goto 5959
      
      donnA = donn
      if (dabs(donn).lt.0.001d0) then
      ddonn = 0.05d0
      else
      ddonn = frdonn*donn
      endif
      if (elecA.lt.0.d0) then
      donnB = donnA+ddonn
      else
      donnB = donnA-ddonn
      endif
      ddd = donnB
 5959 continue
      donnB = ddd
c      write(*,*) 'donnB,ddd = ',donnB,ddd
      sumfdm = 0.d0
      sumfdnm = 0.d0      
      sumfds = 0.d0
      sumfdes = 0.d0
      do i = istp1s,imitt
      sumfdm = ttfdmon(i)*dexp(-(donnB-donn))+sumfdm
      sumfdnm = ttfdnmon(i)*dexp(donnB-donn)+sumfdnm      
      sumfdes = eeblam(i)*dexp(-(donnB-donn))+sumfdes      
      sumfds = eblam(i)*dexp(donnB-donn)+sumfds
      enddo
      elecB = (sumfdm+sumfdes-sumfds-sumfdnm)*dz+surfdens
      ddd = donnA-elecA*(donnB-donnA)/(elecB-elecA)      
c      write(*,*) donn,donnB,elecB
      if (dabs(elecB).lt.delecok) goto 585
c      write(*,*) 'elecA,elecB = ',elecA,elecB
c      write(*,*) donnA,donnB,ddd
c      iop = iop+1
c      write(*,*) 'donnA,donnB = ',donnA,donnB,ddonn
c      write(*,*) (donnB-donnA)/(elecB-elecA),ddd
c      write(*,*) 
c      if (iop.gt.12) stop      
      donnA = donnB
      elecA = elecB
      donnB = ddd
      goto 5959
 585  continue
      
      ddonn = donnB-donn
      ddmax = 0.d0      
      sumfdm = 0.d0
      sumfdnm = 0.d0
      sumfds = 0.d0
      sumfdes = 0.d0      
      do i = istp1s,imitt
      tfdm = ttfdmon(i)*dexp(-ddonn)
      ddiff = abs(tfdm-fdmon(i))/tfdm
      if (ddiff.gt.ddmax) ddmax = ddiff
      fdmon(i) = fdmon(i)*dmm+tdmm*tfdm      
      sumfdm = tfdm+sumfdm
      tfdm = ttfdnmon(i)*dexp(ddonn)
      ddiff = abs(tfdm-fdnmon(i))/tfdm
      if (ddiff.gt.ddmax) ddmax = ddiff
      fdnmon(i) = fdnmon(i)*dmm+tdmm*tfdm      
      sumfdnm = tfdm+sumfdnm
cc      cchat(i) = cchat(i)*dmm+tdmm*ttcchat(i)
cc      achat(i) = achat(i)*dmm+tdmm*ttachat(i)
      cchat(i) = dsqrt(fdmon(i)/(0.5d0*rK*bdpol))
      achat(i) = dsqrt(fdnmon(i)/(0.5d0*rK*bdpol))
      tfds = eblam(i)*dexp(ddonn)
      ddiff = abs(tfds-fdsol(i))/tfds
      if (ddiff.gt.ddmax) ddmax = ddiff
      fdsol(i) = fdsol(i)*dms+tdms*tfds
      sumfds = tfds+sumfds            
      tfdes = eeblam(i)*dexp(-ddonn)
      ddiff = abs(tfdes-fdes(i))/tfdes
      if (ddiff.gt.ddmax) ddmax = ddiff
      fdes(i) = fdes(i)*dms+tdms*tfdes
      sumfdes = tfdes+sumfdes            
      enddo
      elec = (sumfdm+sumfdes-sumfds-sumfdnm)*dz+surfdens      
c      write(*,*) elec,cchat(imitt),tfdm
      
c      stop
      
      jz = imitt+1
c      sum = 0.d0
      do iz = imitt+1,islut      
      jz = jz-1
      fdmon(iz) = fdmon(jz)
      fdnmon(iz) = fdnmon(jz)
      cchat(iz) = cchat(jz)
      achat(iz) = achat(jz)
      fdsol(iz) = fdsol(jz)
      fdes(iz) = fdes(jz)      
c      sum = sum+fdmon(iz)-fdnmon(iz)
      enddo
c      stop
      celec = surfdens+sum*dz
c      write(*,*) celec
c      stop
      donn = donn+ddonn                  
      if (ddmax.lt.ddtol.and.dabs(celec).lt.delecok) goto 200
      goto 100
 200  continue

      sum = 0.d0
      do i = istp1,imitt         
      sum = sum+fdmon(i)-fdnmon(i)-fdsol(i)+fdes(i)
      enddo
      elec = surfdens+sum*dz            
      write(*,*) 'elec,surfdens = ',elec,surfdens
      rewind ifc
      rewind ifnc
      rewind ife
      appsdens = surfdens
      rewind 90
      z = -0.5d0*dz
      sumfp = 0.d0
      sumfm = 0.d0
      sums = 0.d0
      do i = 1,imitt                            
      z = z+dz
      write(ifc,*) z,fdmon(i),cchat(i)
      write(ifnc,*) z,fdnmon(i),achat(i)
      write(ife,*) z,fdsol(i),fdes(i)                  
      write(90,*) z,fp(i),rK*bdpol*ehbclam(i)**2,
     *rK*bdpol*ehbclanm(i)**2
      sumfp = sumfp+fp(i)
      sums = sums+
     *rK*bdpol*ehbclam(i)**2+rK*bdpol*ehbclanm(i)**2
      sumfm = fdmon(i)+fdnmon(i)+sumfm
      enddo
      write(*,*) 'sumfm/sumfp,avr = ',sumfm/sumfp,avr
      write(*,*) 'sums = ',sums
      do i = imitt+1,nfin
      z = z+dz
      write(ifc,*) z,fdmon(i),cchat(i)
      write(ifnc,*) z,fdnmon(i),achat(i)
      write(ife,*) z,fdsol(i),fdes(i)            
      enddo
      
      avfdm = 0.d0
      avfdnm = 0.d0
      avfds = 0.d0
      avfdes = 0.d0
      sumFid = 0.d0
      csumFid = 0.d0
      sumuu = 0.d0
      sumuw = 0.d0
      sumuljm = 0.d0
      altsumF = 0.d0
      sumexFreen = 0.d0
      csumexFreen = 0.d0      
      z = closew-0.5d0*dz
      do iz = istp1s,imitt      
      z = z+dz
      fdm = fdmon(iz)
      fdnm = fdnmon(iz)
      fds = fdsol(iz)
      tfdes = fdes(iz)      
      avfdm = fdm+avfdm
      avfdnm = fdnm+avfdnm
      avfds = fds+avfds
      avfdes = tfdes+avfdes            

      suu = 0.d0
      suw = 0.d0
      do jz = istp1s,isluts      
      iii = iabs(jz-iz)
      suu = (fdmon(jz)-fdnmon(jz)-fdsol(jz)+fdes(jz))*Phimm(iii)+suu
      enddo
      sumuu = 0.5d0*(fdm-fdnm-fds+tfdes)*suu+sumuu
      sumuljm = 0.d0
      sumuw = (fdm-fdnm-fds+tfdes)*csurf+sumuw               
      
      bclamb = 2.d0*(dlog(ehbclam(iz))-scalem)
      belamb = dlog(ebelam(iz))-emscale
      bclanmb = 2.d0*(dlog(ehbclanm(iz))-scalem)
      belanmb = dlog(ebelanm(iz))-emscale
      bslamb = dlog(eblam(iz))-scales
      beslamb = dlog(eeblam(iz))-scalees            

      aabclamb = 2.d0*(dlog(ehbclam(iz))-scalem)+donn
      aabelamb = dlog(ebelam(iz))-emscale+donn
      aabclanmb = 2.d0*(dlog(ehbclanm(iz))-scalem)-donn
      aabelanmb = dlog(ebelanm(iz))-emscale-donn
      aabslamb = dlog(eblam(iz))-scales-donn      
      aabeslamb = dlog(eeblam(iz))-scalees+donn      
      
      cdt = cdtot(iz)*dm3
      pcdt = pis*cdt
      xsi = (1.d0-pcdt)
      rxsi = 1.d0/xsi
      sqrxsi = rxsi*rxsi
      flog = dlog(xsi)
      aex1 = -(c1+1.d0)*flog-0.5d0*(AA1+BB1*pcdt)*pcdt*sqrxsi
      aex2 = -(c2+1.d0)*flog-0.5d0*(AA2+BB2*pcdt)*pcdt*sqrxsi
      exFreen = (fdm*vk+fdnm*vk+fds+tfdes)*aex1
      sumexFreen = exFreen+sumexFreen
      csumexFreen = (fdm*vk+fdnm*vk+fds+tfdes)*aex1+csumexFreen
      sumFid = sumFid+fdm*bclamb+fdnm*bclanmb-fp(iz)+
c      sumFid = sumFid+fdm*bclamb+fdnm*bclanmb-(fdm+fdnm)/avr+      
     *fds*(dlog(fds)-1.d0-chemps)+tfdes*(dlog(tfdes)-1.d0-chempes)      
c      csumFid = csumFid+fdm*bclamb+fdnm*bclanmb-fp(iz)+
      csumFid = csumFid+fdm*bclamb+fdnm*bclanmb-(fdm+fdnm)/avr+
c      csumFid = csumFid+fdm*bclamb+fdnm*bclanmb-fp(iz)+
c     *-(fdm+fdnm)*(1.-ekpp-ekpn)**2/rK+
     *fds*(dlog(fds)-1.d0-chemps)+tfdes*(dlog(tfdes)-1.d0-chempes)      
      altsumF = altsumF+fdm*bclamb+fdnm*bclanmb-fp(iz)+
c      altsumF = altsumF+fdm*bclamb+fdnm*bclanmb-(fdm+fdnm)/avr+  
     *fds*aabslamb-fds+tfdes*aabeslamb-tfdes      
      enddo
      write(*,*) 'z = ',z
      sumexFreen = 2.d0*sumexFreen*dz
      csumexFreen = 2.d0*csumexFreen*dz
      write(*,*) 'sumexFreen, csumexFreen = ',sumexFreen,csumexFreen
      sumF = 2.d0*sumFid*dz+sumexFreen
      csumF = 2.d0*csumFid*dz+sumexFreen
      altsumF = 2.d0*altsumF*dz      
      write(*,*) sumF,csumF
      sumuljm = 2.d0*sumuljm*dz      
      sumuu = 2.d0*sumuu*dz
      sumuw = 2.d0*sumuw*dz
      write(*,*) 'sumuu = ',sumuu
      write(*,*) 'sumuw = ',sumuw
      write(*,*) 'sumuljm = ',sumuljm
      sumF = sumF+sumuu+sumuw
      csumF = csumF+sumuu+sumuw
      altsumF = altsumF-sumuu
      write(*,*) sumF,csumF,altsumF      
      write(*,*) 'sumF = ',sumF
      write(*,*) 'csumF = ',csumF
      write(*,*) 'altsumF = ',altsumF
      omega = sumF
      comega = csumF
      aomega = altsumF      
      write(*,*) 'grand pot. (excl. surf-surf):  ',omega
      write(*,*) omega
      write(*,*) comega
      write(*,*) 'alt. grand pot. (excl. surf-surf):  ',aomega
      write(*,*) aomega      

      uss = -2.d0*pi*surfdens*surfdens*(h+dmw)*rrT                      
      write(*,*) 'u(surf-surf), uss  = ',uss
      totomega = omega+uss
      ctotomega = comega+uss
      atotomega = aomega+uss
      write(*,*)            
      write(*,*) 'total grand pot. (incl. surf-surf):  ',totomega
      write(*,*) totomega+Pbhs*(h+dmw)
      write(*,*) ctotomega+Pbhs*(h+dmw)
      write(*,*) atotomega-2.d0*surfdens*donn+Pbhs*(h+dmw)            
      write(*,*) donn
      write(*,*) surfdens
      write(*,*)      
      write(*,*) 'elec,donn (fin) = ',elec,donn
      write(*,*) 'surfdens = ',surfdens      
      rewind idpot
      write(idpot,*) donn

      fp1S = fdmon(istp1)
      fn1S = fdmon(istp1+2)
      c0Skv = fdmon(istp1+1)
      c1Skv = (fp1S-fn1S)*0.5d0
      c2Skv = (fp1S+fn1S-2.d0*c0Skv)*0.5d0
C      ** EXTRAPOLATION BY USE OF QUADRATIC EXPRESSION **
      fmm = (c0Skv+c1Skv*1.5d0+c2Skv*2.25d0)
      write(*,*) 'cationic monomer contact density - quad. extr.: ',fmm
      fp1S = fdnmon(istp1)
      fn1S = fdnmon(istp1+2)
      c0Skv = fdnmon(istp1+1)
      c1Skv = (fp1S-fn1S)*0.5d0
      c2Skv = (fp1S+fn1S-2.d0*c0Skv)*0.5d0
C      ** EXTRAPOLATION BY USE OF QUADRATIC EXPRESSION **
      fnm = (c0Skv+c1Skv*1.5d0+c2Skv*2.25d0)
      write(*,*) 'anionic monomer contact density - quad. extr.: ',fnm

      fp1S = fdsol(istp1)
      fn1S = fdsol(istp1+2)
      c0Skv = fdsol(istp1+1)
      c1Skv = (fp1S-fn1S)*0.5d0
      c2Skv = (fp1S+fn1S-2.d0*c0Skv)*0.5d0
C      ** EXTRAPOLATION BY USE OF QUADRATIC EXPRESSION **
      fsm = (c0Skv+c1Skv*1.5d0+c2Skv*2.25d0)
      write(*,*) 'solvent contact density - quad. extr.: ',fsm
      fp1S = fdes(istp1)
      fn1S = fdes(istp1+2)
      c0Skv = fdes(istp1+1)
      c1Skv = (fp1S-fn1S)*0.5d0
      c2Skv = (fp1S+fn1S-2.d0*c0Skv)*0.5d0
C      ** EXTRAPOLATION BY USE OF QUADRATIC EXPRESSION **
      fesm = (c0Skv+c1Skv*1.5d0+c2Skv*2.25d0)
      write(*,*) 'co-solvent contact density - quad. extr.: ',fesm      
      Pidc = fmm+fnm+fsm+fesm
      write(*,*) 
      write(*,*) 'ideal part of contact pressure: '
      write(*,'(1e25.14)') Pidc
      write(*,*) 
      write(*,*) 'pressure contributions: '
      write(*,'(3e25.14)') Pidc,-Pbhs,-2.d0*pi*surfdens*surfdens*rrT
      write(*,*) Pidc-Pbhs,-Pbhs-2.d0*pi*surfdens*surfdens*rrT
      write(*,*) Pidc-Pbhs-2.d0*pi*surfdens*surfdens*rrT
      write(*,*) Pidc-Pbhs-2.d0*pi*surfdens*surfdens*rrT
      write(*,*) 
      
      avfdm = avfdm*dz/h
      avfdnm = avfdnm*dz/h
      avfds = avfds*dz/h
      avfdes = avfdes*dz/h            
      write(*,*) 'avfdm = ',avfdm
      write(*,*) 'avfdnm = ',avfdnm
      write(*,*) 'avfds = ',avfds
      write(*,*) 'avfdes = ',avfdes            
      write(*,*) 'donn = ',donn
      write(*,*) donn

      fdonn = donn
      write(*,*) 'fdonn = donn = ',fdonn
      write(*,*) fdonn
       
c     z = zmitt
      z = (h+dm1)/2.d0      
      du = 0.d0
      zp = closew-0.5d0*dz
      write(*,*) 'zmitt,zp(init) = ',z,zp
      do j = istp1,islut      
      zp = zp+dz
      diffz = dabs(zp-z)
      phic = -2.d0*pi*diffz*rrT*dz      
      du = phic*(fdmon(j)-fdnmon(j)-fdsol(j)+fdes(j))+du      
      enddo
      write(*,*) 'zmitt,zp(final) = ',z,zp      
      pmid = du+fdonn+csurf
      write(*,*) 'potential at mid plane: ',pmid,du+donn,du
      write(*,*) pmid
      write(*,*) du+donn+csurf
      write(*,*) du+csurf

      z = 0.d0
      du = 0.d0
      zp = closew-0.5d0*dz      
      write(*,*) 'zlw,zp(init) = ',z,zp
      do j = istp1,islut            
      zp = zp+dz
      diffz = dabs(zp-z)
      phic = -2.d0*pi*diffz*rrT*dz
      du = phic*(fdmon(j)-fdnmon(j)-fdsol(j)+fdes(j))+du
      enddo
      write(*,*) 'zlw,zp(final) = ',z,zp            
      psurf = du+fdonn+csurf
      write(*,*) 'potential at surface: ',du+fdonn,du+donn,du
      write(*,*) psurf
      write(*,*) du+donn+csurf
      write(*,*) du+csurf

      write(*,*) 'net surface potential, psurf-pmid: '
      write(*,*) psurf-pmid

      write(*,*) 'grand potential-surfdens*donn:'
      write(*,*) totomega-surfdens*donn
c      write(*,*) atotomega+surfdens*donn            

      rewind 92
      rewind 94
      summon = 0.d0
      sumsalt = 0.d0
      z = -0.5d0*dz            
      write(*,*) 'z+0.5d0*dz = ',z+0.5d0*dz
      do i = 1,imitt
      z = z+dz
      zp = -0.5d0*dz
      du = 0.d0
      do j = 1,islut            
      zp = zp+dz
      diffz = dabs(zp-z)
      kkk = iabs(j-i)
      du = Phimm(kkk)*(fdmon(j)-fdnmon(j)-fdsol(j)+fdes(j))+du
      enddo
      write(92,*) z,du+csurf
      write(94,*) z,fdmon(i)-fdnmon(i)-fdsol(i)+fdes(i)
      summon = summon+fdmon(i)-fdnmon(i)
      sumsalt = sumsalt-fdsol(i)+fdes(i)
      enddo
      summon = summon*dz
      sumsalt = sumsalt*dz
      write(*,*)'summon, summon+surfdens = ',summon, summon+surfdens
      write(*,*)'sumsalt, sumsalt+surfdens = ',sumsalt, sumsalt+surfdens
      write(*,*) 'summon/sumsalt = ',summon/sumsalt
      write(*,*) 'ddmax,niter: ',ddmax,niter      
      STOP
      END

      subroutine CALCCD
      implicit double precision (a-h,o-z)
      include 'PBlive.inc'
      fp1S = vk*(fdmon(istp1)+fdnmon(istp1))+fdsol(istp1)+fdes(istp1)
      fn1S = vk*(fdmon(istp1+2)+fdnmon(istp1+2))+
     *fdsol(istp1+2)+fdes(istp1+2)
      c0Skv =
     *vk*(fdmon(istp1+1)+fdnmon(istp1+1))+fdsol(istp1+1)+fdes(istp1+1)
      c1Skv = (fp1S-fn1S)*0.5d0
      c2Skv = (fp1S+fn1S-2.d0*c0Skv)*0.5d0
C      ** EXTRAPOLATION BY USE OF QUADRATIC EXPRESSION **
      fwc = (c0Skv+c1Skv*1.5d0+c2Skv*2.25d0)
      if (fwc.lt.0.d0) fwc = 0.d0
      z = sclosew-0.5d0*dz
      iii = min(istp1+ism-1,imitt)      
      do iz = istp1s,imitt
      z = z+dz
      zs = closew
      z1 = zs
      z2 = zs+0.5d0*dz
      f1 = fwc
      f2 = vk*(fdmon(istp1)+fdnmon(istp1))+fdsol(istp1)+fdes(istp1)
      sancint = f1*rthree*(z1-z)**3
      sanciq = (f1+f2)*(z2-z1)
      sanci4 = 2.d0*(f2-f1)*((z2-z)**4-(z1-z)**4)
      zs = z2-dz
      do jz = istp1,iz+ism-2
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      f1 = f2
      f2 = vk*(fdmon(jz+1)+fdnmon(jz+1))+fdsol(jz+1)+fdes(jz+1)
      sanciq = (f1+f2)*(z2-z1)+sanciq
      sanci4 = (f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4
      enddo
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      f1 = f2
      f2 = vk*(fdmon(iz+ism)+fdnmon(iz+ism))+fdsol(iz+ism)+fdes(iz+ism)
      sancint = sancint-f2*rthree*(z2-z)**3+
     *(sanciq+(f1+f2)*(z2-z1))*0.5d0*dm2+
     *((f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4)*rdz*rtwelve
      cdtot(iz) = 0.75d0*sancint*rdm3
      enddo

      do iz = iii+1,imitt      
      z = z+dz
      zs = z-dm1
      z1 = zs
      z2 = zs+dz
      f1 = vk*(fdmon(iz-ism)+fdnmon(iz-ism))+fdsol(iz-ism)+fdes(iz-ism)
      f2 =
     *vk*(fdmon(iz-ism+1)+fdnmon(iz-ism+1))+
     *fdsol(iz-ism+1)+fdes(iz-ism+1)
      sancint = f1*rthree*(z1-z)**3
      sanciq = (f1+f2)*(z2-z1)
      sanci4 = (f2-f1)*((z2-z)**4-(z1-z)**4)
      do jz = iz-ism+1,iz+ism-2
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      f1 = f2
      f2 = vk*(fdmon(jz+1)+fdnmon(jz+1))+fdsol(jz+1)+fdes(jz+1)
      sanciq = (f1+f2)*(z2-z1)+sanciq
      sanci4 = (f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4
      enddo
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      f1 = f2
      f2 = vk*(fdmon(iz+ism)+fdnmon(iz+ism))+fdsol(iz+ism)+fdes(iz+ism)
      sancint = sancint-f2*rthree*(z2-z)**3+
     *(sanciq+(f1+f2)*(z2-z1))*0.5d0*dm2+
     *((f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4)*rdz*rtwelve
      cdtot(iz) = 0.75d0*sancint*rdm3
      enddo
      return
      end

      subroutine AVEC
      implicit double precision (a-h,o-z)
      include 'PBlive.inc'
      do iz = istp1s,imitt
      cdt = cdtot(iz)*dm3
      pcdt = pis*cdt
      xsi = (1.d0-pcdt)
      rxsi = 1.d0/xsi
      sqrxsi = rxsi*rxsi
      flog = dlog(xsi)            
      ae1(iz) = -(c1+1.d0)*flog-0.5d0*(AA1+BB1*pcdt)*pcdt*sqrxsi      
      daex1 = rxsi*(c1+1.d0-0.5d0*AA1*rxsi*(1.d0+2.d0*pcdt*rxsi)-
     *BB1*pcdt*rxsi*(1.d0+pcdt*rxsi))
      convp(iz) =
     *(vk*fdmon(iz)+vk*fdnmon(iz)+fdsol(iz)+fdes(iz))*daex1*pis*dhs3      
      enddo
      jz = imitt+1
      do iz = imitt+1,islut      
      jz = jz-1      
      convp(iz) = convp(jz)
      ae1(iz) = ae1(jz)
c      ae2(iz) = ae2(jz)
      enddo
      return
      end

      subroutine CALCELAM
      implicit double precision (a-h,o-z)
      include 'PBlive.inc'
c      dimension trams(0:maxel)
      fp1S = convp(istp1s)
      c0Skv = convp(istp1s+1)
      fn1S = convp(istp1s+2)
      c1Skv = (fp1S-fn1S)*0.5d0
      c2Skv = (fp1S+fn1S-2.d0*c0Skv)*0.5d0
C      ** EXTRAPOLATION BY USE OF QUADRATIC EXPRESSION **
      fwc = (c0Skv+c1Skv*1.5d0+c2Skv*2.25d0)
      z = closew-0.5d0*dz
      iii = min(istp1s+ism-1,imitt)
      do iz = istp1,iii
      z = z+dz
      zs = sclosew
      z1 = zs
      z2 = zs+0.5d0*dz
      f1 = fwc
      f2 = fp1S
      sancint = f1*rthree*(z1-z)**3
      sanciq = (f1+f2)*(z2-z1)
      sanci4 = 2.d0*(f2-f1)*((z2-z)**4-(z1-z)**4)
      zs = z2-dz
      jjj = iz+ism      
      do jz = istp1s,jjj-2
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      f1 = f2
      f2 = convp(jz+1)
      sanciq = (f1+f2)*(z2-z1)+sanciq
      sanci4 = (f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4
      enddo
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      f1 = f2
      f2 = convp(jjj)
      sancint = sancint-f2*rthree*(z2-z)**3+
     *(sanciq+(f1+f2)*(z2-z1))*0.5d0*dm2+
     *((f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4)*rdz*rtwelve
      trams = sancint*0.75d0*rdhs3
      
      du = 0.d0
      do jz = istp1s,isluts
      kkk = iabs(jz-iz)
      du = Phimm(kkk)*(fdmon(jz)-fdnmon(jz)-fdsol(jz)+fdes(jz))+du      
      enddo
      
      strams = trams+ae1(iz)      
      emtrams = strams*vk
      cmtrams = strams*vk
      eblam(iz) = dexp(-strams+scales+du+donn+csurf)
      eeblam(iz) = dexp(-strams+scalees-du-donn-csurf)                  
c      ebelam(iz) = dexp(-emtrams+emscale-(du+donn+csurf))
      ehbclam(iz) = dexp(-0.5d0*(cmtrams+du+donn+csurf)+scalem)
c      ebelanm(iz) = dexp(-emtrams+emscale+du+donn+csurf)
      ehbclanm(iz) = dexp(-0.5d0*(cmtrams-du-donn-csurf)+scalem)      
      enddo

      do iz = iii+1,imitt
      z = z+dz
      zs = z-dm1
      z1 = zs
      z2 = zs+dz
      f1 = convp(iz-ism)
      f2 = convp(iz-ism+1)
      sancint = f1*rthree*(z1-z)**3
      sanciq = (f1+f2)*(z2-z1)
      sanci4 = (f2-f1)*((z2-z)**4-(z1-z)**4)
      jjj = iz+ism      
      do jz = iz-ism+1,jjj-2
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      f1 = f2
      f2 = convp(jz+1)
      sanciq = (f1+f2)*(z2-z1)+sanciq
      sanci4 = (f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4
      enddo
      zs = zs+dz
      z1 = zs
      z2 = zs+dz
      f1 = f2
      f2 = convp(jjj)
      sancint = sancint-f2*rthree*(z2-z)**3+
     *(sanciq+(f1+f2)*(z2-z1))*0.5d0*dm2+
     *((f2-f1)*((z2-z)**4-(z1-z)**4)+sanci4)*rdz*rtwelve
      trams = sancint*0.75d0*rdhs3

      du = 0.d0
      do jz = istp1s,isluts
      kkk = iabs(jz-iz)
      du = Phimm(kkk)*(fdmon(jz)-fdnmon(jz)-fdsol(jz)+fdes(jz))+du      
      enddo
      
      strams = trams+ae1(iz)      
      emtrams = strams*vk
      cmtrams = strams*vk
      eblam(iz) = dexp(-strams+scales+du+donn+csurf)
      eeblam(iz) = dexp(-strams+scalees-du-donn-csurf)                  
c      ebelam(iz) = dexp(-emtrams+emscale-(du+donn+csurf))
      ehbclam(iz) = dexp(-0.5d0*(cmtrams+du+donn+csurf)+scalem)
c      ebelanm(iz) = dexp(-emtrams+emscale+du+donn+csurf)
      ehbclanm(iz) = dexp(-0.5d0*(cmtrams-du-donn-csurf)+scalem)            
      enddo
      jz = imitt+1
      do iz = imitt+1,islut
      jz = jz-1
c      ebelam(iz) = ebelam(jz)            
      ehbclam(iz) = ehbclam(jz)
c      ebelanm(iz) = ebelanm(jz)
      ehbclanm(iz) = ehbclanm(jz)
      eblam(iz) = eblam(jz)
      eeblam(iz) = eeblam(jz)                        
      enddo      
      return 
      end

      

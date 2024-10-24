c########################################################################
c
c       2D var velocity interpoaltion of HFR (Lyon configuration)
c                arbitrary data locations & weights
c                   X,Y-varying smoothness weights
c                         (u,v)=0 on LAND 
c  Observations:
c  AS: do not put iter along with gnorm in common
c
C########################################################################
       INCLUDE 'par.inc'
       PARAMETER (mdwork=LX*LY*15,lbuf=2*n,ldim=lbuf)
c     N - number of gridpoints in interpolation grid
c     ldim - dimension of the problem (number of unknowns)
c-----------------------------------------------------------------------
       Dimension izs(1),rzs(1)
       Integer Iwork(5)
       Real Rwork(mdwork),w2(lx,ly),gbuf0(lbuf),cbuf(lbuf),gbuf(lbuf)
       real hess(ldim,ldim),ara(nb),ada(nb)
       real*8 ssu,squ
c-----------------------------------------------------------------------
       real u(lx,ly),v(lx,ly),ut(lx,ly),vt(lx,ly),uu(nb)
       real rot0(lx,ly),div0(lx,ly),rot(lx,ly),div(lx,ly)
       integer msk1(lx,ly),msp(lx,ly),remPoints(nb)
       common /norm/gnorm
C=======================================================================
      external siml, EUCLID, CTONBE, CTCABE
      void=999.; voi=0.; pi=4.*atan(1.)
      open(1,file='ass100.go1')
c--------------------------------INPUT----------------------------------
      open(2,file='my.dat')
      read(2,13) wei_dat0   ! data weight
      read(2,13) wei_div0   ! divergence smoothing weight
      read(2,13) wei_rot0   ! curl smoothing weight
      read(2,13) WUV        ! velocity smoothing weight
      read(2,13) roh        ! search radius to calculate number of
!                             obs. points around the given gridpoint to
!                             compute spatial dependence of weights (dat,div,rot)
      read(2,14) jswd; read(2,14) Kiter
      read(2,14) Iprint; read(2,14) ihess; 
      read(2,'(10x,f5.3)') xnoise   ! noise v twin-data experiments
      close(2)
13    format(10x,e9.2)
14    format(10x,i5)
      write(6,'(a,4e9.2)') 'WEIGHTS-->',wei_div0,wei_rot0,WUV
c=================== read mask & set the domain =======================
      open(22,file='mask.dat'); open(23,file='mask11.dat')
      do j=ly,1,-1
      read(22,*) (msk(i,j),i=1,lx); read(23,*) (msk1(i,j),i=1,lx); 
      enddo
      msp=1
      do i=1,lx; do j=1,ly
       if(msk1(i,j).ne.1) msk1(i,j)=0
       if(msk (i,j).ne.1) msp(i,j)=0
      enddo; enddo
c--------- control & harmonic masks ------------------------
       msc=void; msh=0
       do i=1,lx; do j=1,ly
       ip=min(lx,i+1); im=max(1,i-1); jp=min(ly,j+1); jm=max(1,j-1)
        if(msk(i,j).eq.1) then
           msc(i,j)=1 
           if(msk(ip,j)*msk(im,j)*msk(i,jm)*msk(i,jp).eq.1) msh(i,j)=1       
        endif
        if(msk(i,j).eq.2) then
        if(msk(ip,j).eq.1.or.msk(im,j).eq.1.or.
     >     msk(i,jp).eq.1.or.msk(i,jm).eq.1) msc(i,j)=1
        endif
       enddo; enddo      
c--------------------- check the masks ------------------
       do i=1,lx; do j=1,ly
       if(msc(i,j).eq.void) msc(i,j)=0
       enddo; enddo
       nn=2*sum(msc)
       write(6,*) 'CONTROL  DIMENSION--->',nn
       write(6,*) 'div-curl DIMENSION--->',sum(msp)
       write(6,*) '   harm  DIMENSION--->',sum(msh)
       call map0(msk*1.,lx,ly,0.,'    mask')
       call map0(msc*1.,lx,ly,0.,'cntrmask')
c       call map0(msh*1.,lx,ly,0.,'harmmask')
c....... check control dimension ...............        
       if(nn.ne.lbuf) then
       write(6,*) 
     >'set N in PAR.INC equal to the number of control:  lbuf=',lbuf/2,
     >'  cntrl:',nn/2
       stop
      endif
c=============== read Tonkin data file ===============================
c   PROVERIT ZDESSSS koord
       x0=+0; y0=+0  ! Testing for COF site, origin [153.010561,-31.609578]
      !  dx=5.8    ! For NEWC
       dx=1.7468094825772     ! For COF
       dy=1.504578552427569   ! For COF
       open(22,file='out12.dat')
c       open(22,file='out540.dat')
       i=0; 
       do kk=1,25000
       read(22,*,end=700) sii,css,x7,y7,udata7
         if(udata7.ne.999.) then
         xt=(x7-x0)/dx; yt=(y7-y0)/dy
         jx=int(xt)+1; jxp=jx+1  
         jy=int(yt)+1; jyp=jy+1  
        IF(xt.gt.0..and.yt.gt.0..and.xt.lt.lx1*1..and.yt.lt.ly1*1.)THEN
        j1=msc(jx,jy)+msc(jxp,jy)+msc(jx,jyp)+msc(jxp,jyp)
        j2=msc(jx,jy)*msk(jx,jy) + msc(jxp,jy)*msk(jxp,jy) +
     >     msc(jx,jyp)*msk(jx,jyp) + msc(jxp,jyp)*msk(jxp,jyp)
             if(j1.eq.4.or.j2.le.4) then
          i=i+1
		  remPoints(i) = kk
          si(i)=sii; cs(i)=css; x(i)=xt; y(i)=yt; udata(i)=udata7
c.................. vortex field for check .......................
!         xc=18.; yc=18.  ! koordinaty centra vihrya
!         wid=6.          ! polushirina vihrya
!         xv=x(i)-xc; yv=y(i)-yc; ee=exp(-(xv*xv+yv*yv)/(2*wid*wid))*1.e3
!         uuk=yv/(2*wid*wid)*ee; vvk=-xv/(2*wid*wid)*ee
!          udata(i)=uuk*cs(i)+vvk*si(i)
!          ara(i)=uuk; ada(i)=vvk
!.................................................................
          endif; 
          endif; endif
       enddo
700    close(22)
       nb7=i
c====================================================================
       write(6,*) 'number of observation points -->',nb7
c===== compute data weights (inversely proportional to obs.point density) ====
       wei_dat=0.; 
       do k=1,nb7
        wei_dat(k)=wei_dat0*10; 
       !    kl=0
       !    do m=1,nb7
       !      r=sqrt( (x(k)-x(m))**2 + (y(k)-y(m))**2 )   !   count obs. points kl
       !      if(r.lt.1.) kl=kl+1                         !   within the serach radius "1"
       !    enddo
       !if(k.ne.0) wei_dat(k)=wei_dat(k)/kl               !   data weight inversely proportional to kl
       !write(78,*) wei_dat(k)          
       enddo
c===== compute smoothness weights (inversely proportional to obs.point density) ====
       wdiv=0.; wrot=0.
       do i=1,lx; do j=1,ly
       if(msh(i,j).eq.1) then
        wdiv(i,j)=wei_div0; wrot(i,j)=wei_rot0 
       !    k=0
       !    do m=1,nb7
       !      r=sqrt( (i-x(m))**2 + (j-y(m))**2 )   !   count obs. points
       !      if(r.lt.roh) k=k+1                    !   within the serach radius roh
       !    enddo
       !if(k.ne.0) then
       !   wrot(i,j)=wrot(i,j)/k; 
       !   wdiv(i,j)=wdiv(i,j)/k; endif
       endif
       enddo; enddo
c       call map0(wrot,lx,ly,0.,'wei_curl')
c       call map0(wdiv,lx,ly,0.,'wei_divv')
c=================== compute bilinear interpolation weights =====================
      do i=1,nb7
        ix(i)=int(x(i)); iy(i)=int(y(i)); 
        xx=x(i)-ix(i);   yy=y(i)-iy(i)
        wei(i,1)=(1-xx)*(1-yy)
        wei(i,2)=(1-xx)*yy
        wei(i,3)=xx*(1-yy)
        wei(i,4)=xx*yy
      enddo
c-------------------------- check coordinates -----------------------      
c      ixmin=10000; ixmax=-ixmin; iymin=ixmin; iymax=ixmax
c       do i=1,nb7
c       if(ix(i).lt.ixmin) ixmin=ix(i); if(ix(i).gt.ixmax) ixmax=ix(i)
c       if(iy(i).lt.iymin) iymin=iy(i); if(iy(i).gt.iymax) iymax=iy(i)
c       enddo
c      write(6,*)'xmin,xmax,ymin,ymax-->',ixmin,ixmax,iymin,iymax,lx1,ly1     
c======================= data statistics ============================
       um=sum(udata(1:nb7))/nb7
       ud=0.
       do i=1,nb7
       ud=ud+(udata(i)-um)**2
       enddo
       ud=sqrt(ud/nb7) !AS
       write(6,*) 'mean of the radial velocities--->',um
       write(6,*) 'RMS  of the radial velocities--->',ud
       write(6,*) 'anticipated relative error--->',.05/(ud*sqrt(2.))

c=========== output to draw radial velocities (for velbc MAtlab) =======================
c       open(22,file='v_obs.dat')
       open(23,file='u_obs.dat')
       do i=1,nb7
       write(23,'(4f10.4)') udata(i)*cs(i),udata(i)*si(i),x(i),y(i)
c       write(22,'(4f10.4)') ara(i),ada(i),x(i),y(i)      !   istinnye skorosti dlya "otladochnogo vihrya"
       enddo
       close(22); close(23)
c---------------- ~ ZERO FIRST GUESS FIELDS --------------------
      U=0.; V=0; a=.005
       do i=1,lx; do j=1,ly
        u(i,j)=sin(i*j*3.)*msk1(i,j)*a; !!!!! arbitrary
        v(i,j)=cos(i*j*2.)*msk1(i,j)*a  !!!!! fields for Lagr. Identity test
       enddo; enddo
c.................. Lagrangian identity test .................
!      udata=0.
c      do i=1,lx
c      do j=1,ly
c      if(msk(i,j).eq.1) then
c        ph(i,j)=1.
       call cbuf1(u,v,cbuf,lbuf,1)   ! u,v->cbuf         ! keep 'em
       call siml(indic,lbuf,cBuf,fun,gBuf,izs,rzs,dzs)   ! unchanged
c
c      CBUF1 - perevodit u,v --> dlinnyy odnomernyy vektor (1) i obratno (0)
c      SIML  - po dannomu dlinnomu vektoru CBUF rasschityvaet funkcional (FUN) i ego gradient (GBUF)
c
!       squ=0.
!       do k=1,lbuf
!       squ=squ+gbuf(k)*cbuf(k)*.5
!       enddo
!       write(6,*) 'first sum--->',squ
!       call VEL(u,v,uu)
!       ssu=0.d0
!       do k=1,nb7
!       ssu = ssu + udata(k)*wei_dat(k)*(uu(k)-udata(k))*.5
!       enddo
!       write(6,*) 'second sum--->',ssu
!       fun1=squ-ssu
!       err=2*abs(fun1-fun)/(abs(fun+fun1)+1.e-34)
c       if(abs(err).gt.1.e-9) then
!       write(6,*) 'they gonna be equal:',fun
!       write(6,*) '                    ',fun1
!       write(6,*) '                    ',err !,i,j; pause
c       endif
c         ph(i,j)=0.
c      endif
c      enddo
c      enddo
!       STOP
c-----------------------------------------------------------------------
      Iter=0
c---- optimizer call ---------------------------------------------------
      impres = 1;   io     = 1
      mode   = 0
      dxmin  = 1.0e-10
      epsg   = 5.0e-6
      df1    = 0.5*fun
      itera  = Kiter;    isimu  = Kiter*1.1
      ndwork = mdwork
c------------------------------------------------------------------
      call M1QN3(Siml, EUCLID, CTONBE, CTCABE,
     >           lbuf, cBuf,  Fun, gBuf, dxmin, df1, epsg, impres,
     >           io,  mode, itera, isimu, Iwork, Rwork, ndwork,
     >           izs, rzs, dzs )

      call cbuf1(u,v,cbuf,lbuf,0)   ! perevod rezul'tata v dvumernyy massiv U,V

c........... velocity in measurement points ................
	  open(138,file='err.dat')
      CALL VEL(u,v,uu)           ! raschet radialnyh skorostey v tochkah izmereniy
      u=u*msk1; v=v*msk1         !  zamaskirovali kraya          
        eru=0.; err=0.
        do k=1,nb7
        ii=ix(k); jj=iy(k); ip=ii+1; jp=jj+1
        i9=msk1(ii,jj)*msk1(ip,jj)*msk1(ii,jp)*msk1(ip,jp)
        if(i9.eq.1) then
        eru=eru+abs(uu(k)-udata(k)); err=err+abs(udata(k)); endif
	    write(138,'(i11,6f10.4)') remPoints(k),cs(k),si(k),x(k),y(k),udata(k),uu(k)
        enddo
        eru=eru/err
      write(6,*) 'error at obs. points--->',eru
	  close(138)
c............................
      call map0(u,lx,ly,void,'U_opt   ')
      call map0(v,lx,ly,void,'V_opt   ')
      write(6,*) 'total iterations-->',itera
C================ CHECK FINAL ERRORS ============================
c................ velocity errors .................
!      varu=(sum(abs(ut)*msk1)+sum(abs(vt)*msk1))*.5
!      eru=sum(abs(u-ut)*msk1)/varu;     erv=sum(abs(v-vt)*msk1)/varu
!       write(6,'(a,2f8.3,a,f5.2)')
!     >' VEL-errors:', eru, erv,' Noise level:',xnoise
c............ divergence curl ...................
      div0=0.; rot0=0.
      do i=2,lx1; do j=2,ly1
      if(msk(i,j).eq.1) then
      rot0(i,j)=(u(i,j-1)-u(i,j+1) + v(i+1,j)-v(i-1,j))*msk1(i,j)
      div0(i,j)=(u(i+1,j)-u(i-1,j) + v(i,j+1)-v(i,j-1))*msk1(i,j)
      endif
      enddo; enddo
c.........................................................
       call map0(rot0,lx,ly,void,'rotor_op')
       call map0(div0,lx,ly,void,'diver_op')
c------------------------ dump to draw in Matlab (velbc) ------------------------
       open(24,file='div_opt.dat'); open(25,file='rot_opt.dat')
       do i=1,lx; write(24,'(360f10.6)') (div0 (i,j),j=1,ly); enddo       
       do i=1,lx; write(25,'(360f10.6)') (rot0 (i,j),j=1,ly); enddo
       close(24); close(25)

       open(24,file='u_opt.dat'); open(25,file='v_opt.dat')
       do i=1,lx; write(24,'(360f10.6)') (u(i,j),j=1,ly); enddo       
       do i=1,lx; write(25,'(360f10.6)') (v(i,j),j=1,ly); enddo
       close(24); close(25)
c========================== HESSIAN ==========================
       IF(ihess.eq.1) THEN
       pert=.1; iprint=10000
       call cbuf1(u,v,cbuf,lbuf,1)        
       call siml(indic,lbuf,cBuf,fun,gBuf0,izs,rzs,dzs)
       gbuf0=gbuf
       do i=1,ldim
       if(i/iprint*iprint.eq.i) write(6,*) 'Hessian row:',i
         cbuf(i)=cbuf(i)+pert
         call siml(indic,lbuf,cBuf,fun,gBuf,izs,rzs,dzs)
         cbuf(i)=cbuf(i)-2*pert
         call siml(indic,lbuf,cBuf,fun,gBuf0,izs,rzs,dzs)
         cbuf(i)=cbuf(i)+pert
            do j=1,ldim
            hess(i,j)=(gbuf(j)-gbuf0(j))/(2*pert)
            enddo
       enddo
         a=0.
         do i=1,ldim
         do j=i+1,ldim
      o=(abs(hess(i,i))+abs(hess(j,j)))
          if(o.ne.0.) a=a + 2.*abs(hess(i,j)-hess(j,i))/o
         enddo
         enddo
        write(6,*) 
     >   'hessian assymmtry index-->',2*a/((ldim-1)*ldim)

       write(6,*) 'ldim-->',ldim
       open(22,file='hess.dat',form='unformatted')
       write(22) hess,msc,msk
c       write(22) msc,msk
       close(22)

c       open(22,file='msc.dat',form='unformatted')
c       write(22) msc,msk
c       close(22)

       ENDIF
       call map0(msk*1.,lx,ly,0.,'    mask')
       call map0(msc*1.,lx,ly,0.,'cntrmask')       
       STOP
       END
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
       Subroutine siml(indic,lbuf,cBuf,fun,gBuf,izs,rzs,dzs)
       INCLUDE 'par.inc'
       REAL cbuf(lbuf),gbuf(lbuf),uu(nb)
       REAL ps(lx,ly),ph(lx,ly),w1(lx,ly),w2(lx,ly)
       REAL smou(lx,ly),smov(lx,ly),wsu(lx,ly),wsv(lx,ly)
       REAL rot(lx,ly),div(lx,ly),rotc(lx,ly),divc(lx,ly),ws(lx,ly)
       data iter/0/
       real*8 fun_harm,fun_dat
       common /norm/gnorm
       wsu=0.; wsv=0.; w1=0.; w2=0.;
c------------------- func value ------------------------
       ps=void; ph=void;
       call cbuf1(ps,ph,cbuf,lbuf,0)
c------------------- data-func -------------------------
       CALL VEL(ps,ph,uu)      !  UU --> radial velocities v tochkah zmereniy
       fun_dat=0.d0
       do i=1,nb7
       fun_dat=fun_dat + wei_dat(i)*(uu(i)-udata(i))**2
       enddo
c------------------ div-rot & UV-smoothness func -----------------------
      div=0.; rot=0.; smou=0.; smov=0.; fdiv=0.; frot=0.
      do i=2,lx1; do j=2,ly1
      ip=i+1; im=i-1; jp=j+1; jm=j-1
      if(msk(i,j).eq.1) then
      rot(i,j)=ps(i,jm)-ps(i,jp) + ph(ip,j)-ph(im,j)
      div(i,j)=ps(ip,j)-ps(im,j) + ph(i,jp)-ph(i,jm)
      smou(i,j)=-4.*ps(i,j)+ps(ip,j)+ps(im,j)+ps(i,jp)+ps(i,jm)
      smov(i,j)=-4.*ph(i,j)+ph(ip,j)+ph(im,j)+ph(i,jp)+ph(i,jm)
      endif
      enddo; enddo
c------------------ div-rot & UV smoothness func ------------
      w1=0.; w2=0.; fdiv=0.; frot=0.; fsmo=0.
      do i=2,lx1; do j=2,ly1
      ip=i+1; im=i-1; jp=j+1; jm=j-1
      if(msh(i,j).eq.1) then
        w1(i,j)=-4.*rot(i,j)+rot(im,j)+rot(ip,j)+rot(i,jm)+rot(i,jp)
        frot=frot + wrot(i,j)*w1(i,j)**2
        w2(i,j)=-4.*div(i,j)+div(im,j)+div(ip,j)+div(i,jm)+div(i,jp)
        fdiv=fdiv + wdiv(i,j)*w2(i,j)**2
            
      endif
         if(msk(i,j).eq.1) Fsmo=Fsmo + WUV*(smou(i,j)**2 +smov(i,j)**2)
      enddo; enddo         
      FUN=(fun_dat+fdiv+frot+Fsmo)*0.5
c-----------------adjoint of smoothness --------------------------
      rot=0.; div=0.
      do i=2,lx1; do j=2,ly1
      ip=i+1; im=i-1; jp=j+1; jm=j-1
      if(msk(i,j).eq.1) then
        rot(i,j)=-4.*w1(i,j)*wrot(i,j) + w1(im,j)*wrot(im,j) +
     >    w1(ip,j)*wrot(ip,j)+w1(i,jm)*wrot(i,jm)+w1(i,jp)*wrot(i,jp)
        div(i,j)=-4.*w2(i,j)*wdiv(i,j) + w2(im,j)*wdiv(im,j) +
     >    w2(ip,j)*wdiv(ip,j)+w2(i,jm)*wdiv(i,jm)+w2(i,jp)*wdiv(i,jp)
      wsu(i,j)=WUV*
     >(-4.*smou(i,j)+smou(ip,j)+smou(im,j)+smou(i,jp)+smou(i,jm))
      wsv(i,j)=WUV*
     >(-4.*smov(i,j)+smov(ip,j)+smov(im,j)+smov(i,jp)+smov(i,jm))
      endif
      enddo; enddo  
c----------------- div-rot-gradient ------------------------------
      rotc=0.; divc=0.
      do i=1,lx; do j=1,ly
      ip=min(lx,i+1); im=max(1,i-1); jp=min(ly,j+1); jm=max(1,j-1)
      if(msc(i,j).eq.1) then
      rotc(i,j)=-rot(i,jm) + rot(i,jp) -div(ip,j) + div(im,j)
      divc(i,j)=-rot(ip,j) + rot(im,j) -div(i,jp) + div(i,jm)
c      rotc(i,j)=-rot(i,jm)*wrot(i,jm) + rot(i,jp)*wrot(i,jp)
c     >          -div(ip,j)*wdiv(ip,j) + div(im,j)*wdiv(im,j)
c      divc(i,j)=-rot(ip,j)*wrot(ip,j) + rot(im,j)*wrot(im,j) 
c     >          -div(i,jp)*wdiv(i,jp) + div(i,jm)*wdiv(i,jm)
      endif
      enddo; enddo      
c------------------ Observation gradient -----------------------  
       W1=0.; W2=0.     
       DO ii=1,nb7
             i0=ix(ii)+1; ip=i0+1  
             j0=iy(ii)+1; jp=j0+1  
c................ U - gradients .................
      du=(uu(ii)-udata(ii))*wei_dat(ii)*cs(ii)     
      w1(i0,j0)=w1(i0,j0) + du*wei(ii,1)
      w1(i0,jp)=w1(i0,jp) + du*wei(ii,2)
      w1(ip,j0)=w1(ip,j0) + du*wei(ii,3)
      w1(ip,jp)=w1(ip,jp) + du*wei(ii,4)
c................ V - gradients .................      
      dv=(uu(ii)-udata(ii))*wei_dat(ii)*si(ii)
      w2(i0,j0)=w2(i0,j0) + dv*wei(ii,1)
      w2(i0,jp)=w2(i0,jp) + dv*wei(ii,2)
      w2(ip,j0)=w2(ip,j0) + dv*wei(ii,3)
      w2(ip,jp)=w2(ip,jp) + dv*wei(ii,4)
      ENDDO
c---------------------------------------------------------------
100   w1=w1+rotc + wsu; w2=w2+divc + wsv        
      call cbuf1(w1,w2,gbuf,lbuf,1)
c-------------- diagnostics ------------------------------------
      if(iter.eq.0) gnorm=sqrt(sum(gbuf*gbuf))
      IF(iter/iprint*iprint.eq.iter) THEN
          gn=sqrt(sum(gbuf*gbuf))/gnorm
      write(6,12) 
     >' FUNC-->',fdiv,frot,Fsmo,fun_dat,'  GRAD',gn,' ITER-->',iter
      write(1,12) 
     >' FUNC-->',fdiv,frot,Fsmo,fun_dat,'  GRAD',gn,' ITER-->',iter
      ENDIF
12     format(a,4e11.3,a,f8.5,a,i5)        
c--------------------------------------------------------------
        iter=iter+1
        return
        end
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
       SUBROUTINE VEL(u7,v7,uu)
       include 'par.inc'
!       COMMON/RADAR/ bear(nb),radi(nb),u9(nb),v9(nb),nend(nradars)
       real u7(lx,ly),v7(lx,ly),uu(nb)
       do i=1,nb7
       ii=ix(i)+1; ip=ii+1  
       jj=iy(i)+1; jp=jj+1  
        u1= u7(ii,jj)*wei(i,1) + u7(ii,jp)*wei(i,2) +
     >      u7(ip,jj)*wei(i,3) + u7(ip,jp)*wei(i,4) 
        v1= v7(ii,jj)*wei(i,1) + v7(ii,jp)*wei(i,2) +
     >      v7(ip,jj)*wei(i,3) + v7(ip,jp)*wei(i,4)
       uu(i) = u1*cs(i) + v1*si(i)
       enddo
       return
       end
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
       subroutine cbuf1(a,b,cbuf,lbuf,index)
       INCLUDE 'par.inc'
       real a(lx,ly),b(lx,ly),cbuf(lbuf)
       data icall/0/
       ll=1
       IF(index.eq.1) THEN !=============== 2d-->1d ==================
          do j=1,ly; do i=1,lx
          if(msc(i,j).eq.1) then
          cbuf(ll)=a(i,j);  ll=ll+1
          cbuf(ll)=b(i,j);  ll=ll+1
          endif
          enddo; enddo
       ELSE                !=============== 1d-->2d ==================
          do j=1,ly; do i=1,lx
          if(msc(i,j).eq.1) then         
          a(i,j)=cbuf(ll);  ll=ll+1
          b(i,j)=cbuf(ll);  ll=ll+1
          else 
c          a(i,j)=void; b(i,j)=void
           a(i,j)=0.; b(i,j)=0.
          endif
          enddo; enddo
c.......... shoreline ........
       nn=0
       do i=1,lx; do j=1,ly
       if(msk(i,j).eq.0) then
          a(i,j)=0.; b(i,j)=0.; nn=nn+1
       endif
       enddo; enddo 
       ENDIF
       icall=1
       return
       end          

C******************************************************************
C               THIS ROUTINE PLOTS D1(LX,LY)
C         SIMB1 --> field name;  VOID --> data absence
C******************************************************************
      SUBROUTINE MAP0(D1,LX,LY,VOID,SIMB1)
      PARAMETER (LXX=160,LXP=2*LXX)
      DIMENSION D1(2),AA(LXP),SMB(12)
C
      CHARACTER*8 SIMB1
      CHARACTER*1 SMB,BLANC,AA,XMAT,COAST
      DATA SMB/'+','9','8','7','6','5','4','3','2','1','0','-'/
c********** to change size change figures only here ***************
      DATA BLANC/' '/,XMAT/'*'/,AA/LXP*' '/,VOI/0./,COAST/'M'/
      ISIZE=LXX
 20   FORMAT(1X, 160(A1))
c******************************************************************
      LXM=LX
      IF(LXM.GT.ISIZE) LXM=ISIZE
      LL=LX*LY
C_______ CALCULATION OF MIN & MAX ______________________________
      DO 1 L=1,LL
      IF(D1(L).NE.VOID) GO TO 2
 1    CONTINUE
 2    CONTINUE
C______________
      PMAX1=D1(L)
      PMIN1=PMAX1
C______
      DO 3 L=1,LL
      DL1=D1(L)
      IF(DL1.GT.PMAX1.AND.DL1.NE.VOID.AND.DL1.NE.VOI) PMAX1=DL1
      IF(DL1.LT.PMIN1.AND.DL1.NE.VOID.AND.DL1.NE.VOI) PMIN1=DL1
   3  CONTINUE
C________ PRINT MIN & MAX
      M=1
      write(1,10) SIMB1,PMAX1,PMIN1
   10 FORMAT(2X,A8,':',2G12.4)
      S1=(PMAX1-PMIN1)*0.1
C-----------------------------------------------------------
      DO 5 J=LY,1,-1
      DO 6 I1=1,LXM
      I2=ISIZE+I1
      L3=I1+(J-1)*LX
      AA(I1)=BLANC
      AA(I2)=BLANC
      IF(D1(L3).EQ.VOID) AA(I1)=COAST
      IF(D1(L3).EQ.VOI ) AA(I1)=XMAT
 6    CONTINUE
      DO 106 I1=1,LXM
      L3=I1+(J-1)*LX
      I2=I1+ISIZE
      DL1=D1(L3)
      IF(DL1.EQ.PMIN1.AND.DL1.NE.VOID.AND.DL1.NE.VOI) AA(I1)=SMB(12)
      P12=PMIN1+S1
      P11=PMIN1
      DO 7 K=1,10
      IF(DL1.GT.P11.AND.DL1.LE.P12.AND.DL1.NE.VOID.AND.DL1.NE.VOI) 
     >AA(I1)=SMB(12-K)
      P11=P12
      P12=P12+S1
    7 CONTINUE
      IF(DL1.GE.PMAX1.AND.DL1.LT.PMAX1+.01*S1
     >.AND.DL1.NE.VOID.AND.DL1.NE.VOI) AA(I1)=SMB(1)
 106  CONTINUE
      WRITE(1,20) (AA(I),I=1,ISIZE )
  5   CONTINUE
C---------------------------------------------------------------------
      DO 6666 I=1,ISIZE
 6666 AA(I)=BLANC
      RETURN
      END

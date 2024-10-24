      subroutine m1qn3 (simul,prosca,ctonb,ctcab,n,x,f,g,dxmin,df1, 
     /                  epsg,impres,io,mode,niter,nsim,iz,rz,nrz, 
     /                  izs,rzs,dzs) 
c---- 
c 
c     M1QN3, Version 2.0, March 1993 
c     Jean Charles Gilbert, Claude Lemarechal, INRIA. 
c 
c     M1qn3 has two running modes: the SID (Scalar Initial Scaling) mode 
c     and the DIS (Diagonal Initial Scaling) mode. Both do not require 
c     the same amount of storage, the same subroutines, ... 
c     In the description below, items that differ in the DIS mode with 
c     respect to the SIS mode are given in brakets. 
c 
c     Use the following subroutines: 
c         M1QN3A 
c         DD, DDS 
c         MLIS0 + ECUBE (Dec 88) 
c         MUPDTS, YSTBL. 
c 
c     The following routines are proposed to the user in case the 
c     Euclidean scalar product is used: 
c         EUCLID, CTONBE, CTCABE. 
c 
c     La sous-routine M1QN3 est une interface entre le programme 
c     appelant et la sous-routine M1QN3A, le minimiseur proprement dit. 
c 
c     Le module PROSCA est sense realiser le produit scalaire de deux 
c     vecteurs de Rn; le module CTONB est sense realiser le changement 
c     de coordonnees correspondant au changement de bases: base 
c     euclidienne -> base orthonormale (pour le produit scalaire 
c     PROSCA); le module CTBAB fait la transformation inverse: base 
c     orthonormale -> base euclidienne. 
c 
c     Iz is an integer working zone for M1QN3A, its dimension is 5. 
c     It is formed of 5 scalars that are set by the optimizer: 
c         - the dimension of the problem, 
c         - a identifier of the scaling mode, 
c         - the number of updates, 
c         - two pointers. 
c 
c     Rz est la zone de travail pour M1QN3A, de dimension nrz. 
c     Elle est subdivisee en 
c         3 [ou 4] vecteurs de dimension n: d,gg,[diag,]aux 
c         m scalaires: alpha 
c         m vecteurs de dimension n: ybar 
c         m vecteurs de dimension n: sbar 
c 
c     m est alors le plus grand entier tel que 
c         m*(2*n+1)+3*n .le. nrz [m*(2*n+1)+4*n .le. nrz)] 
c     soit m := (nrz-3*n) / (2*n+1) [m := (nrz-4*n) / (2*n+1)]. 
c     Il faut avoir m >= 1, donc nrz >= 5n+1 [nrz >= 6n+1]. 
c 
c     A chaque iteration la metrique est formee a partir d'un multiple 
c     de l'identite [d'une matrice diagonale] D qui est mise a jour m 
c     fois par la formule de BFGS en utilisant les m couples {y,s} les 
c     plus recents. 
c 
c---- 
c 
c         arguments 
c 
      integer n,impres,io,mode,niter,nsim,iz(5),nrz,izs(1) 
      real x(1),f,g(1),dxmin,df1,epsg,rz(1),rzs(1) 
      double precision dzs(1) 
      external simul,prosca,ctonb,ctcab 
c 
c         variables locales 
c 
      logical inmemo,sscale 
      integer ntravu,id,igg,idiag,iaux,ialpha,iybar,isbar,m,mmemo 
      real r1,r2 
      double precision ps 
c 
c---- impressions initiales et controle des arguments 
c 
         if (impres.ge.1) 
     /    write (io,900) n,dxmin,df1,epsg,niter,nsim,impres 
900   format (/' M1QN3 (Version 2.0, March 1993): entry point'/ 
     /    5x,'dimension of the problem (n):',i6/ 
     /    5x,'absolute precision on x (dxmin):',e9.2/ 
     /    5x,'expected decrease for f (df1):',e9.2/ 
     /    5x,'relative precision on g (epsg):',e9.2/ 
     /    5x,'maximal number of iterations (niter):',i6/ 
     /    5x,'maximal number of simulations (nsim):',i6/ 
     /    5x,'printing level (impres):',i4) 
      if (n.le.0.or.niter.le.0.or.nsim.le.0.or.dxmin.le.0..or.epsg.le.0. 
     /    .or.epsg.gt.1..or.mode.lt.0.or.mode.gt.3) then 
          mode=2 
          if (impres.ge.1) write (io,901) 
901       format (/' >>> m1qn3: inconsistent call') 
          return 
      endif 
c 
c---- what method 
c 
      if (mod(mode,2).eq.0) then 
          if (impres.ge.1) write (io,920) 
  920     format (/' m1qn3: Diagonal Initial Scaling mode') 
          sscale=.false. 
      else 
          if (impres.ge.1) write (io,921) 
  921     format (/' m1qn3: Scalar Initial Scaling mode') 
          sscale=.true. 
      endif 
c 
      if ((nrz.lt.5*n+1).or.((.not.sscale).and.(nrz.lt.6*n+1))) then 
          mode=2 
          if (impres.ge.1) write (io,902) 
902       format (/' >>> m1qn3: not enough memory allocated') 
          return 
      endif 
c 
c---- Compute m 
c 
      call mupdts (sscale,inmemo,n,m,nrz) 
c 
c     --- Check the value of m (if (y,s) pairs in core, m will be >= 1) 
c 
      if (m.lt.1) then 
          mode=2 
          if (impres.ge.1) write (io,9020) 
 9020     format (/' >>> m1qn3: m is set too small in mupdts') 
          return 
      endif 
c 
c     --- mmemo = number of (y,s) pairs in core memory 
c 
      mmemo=1 
      if (inmemo) mmemo=m 
c 
      ntravu=2*(2+mmemo)*n+m 
      if (sscale) ntravu=ntravu-n 
      if (impres.ge.1) write (io,903) nrz,ntravu,m 
903   format (/5x,'allocated memory (nrz) :',i8/ 
     /         5x,'used memory :           ',i8/ 
     /         5x,'number of updates :     ',i8) 
      if (nrz.lt.ntravu) then 
          mode=2 
          if (impres.ge.1) write (io,902) 
          return 
      endif 
c 
      if (impres.ge.1) then 
          if (inmemo) then 
              write (io,907) 
          else 
              write (io,908) 
          endif 
      endif 
907   format (5x,'(y,s) pairs are stored in core memory') 
908   format (5x,'(y,s) pairs are stored by the user') 
c 
c---- cold start or warm restart ? 
c     check iz: iz(1)=n, iz(2)=(0 if DIS, 1 if SIS), 
c               iz(3)=m, iz(4)=jmin, iz(5)=jmax 
c 
      if (mode/2.eq.0) then 
          if (impres.ge.1) write (io,922) 
      else 
          iaux=0 
          if (sscale) iaux=1 
          if (iz(1).ne.n.or.iz(2).ne.iaux.or.iz(3).ne.m.or.iz(4).lt.1 
     &        .or.iz(5).lt.1.or.iz(4).gt.iz(3).or.iz(5).gt.iz(3)) then 
              mode=2 
              if (impres.ge.1) write (io,923) 
              return 
          endif 
          if (impres.ge.1) write (io,924) 
      endif 
  922 format (/' m1qn3: cold start') 
  923 format (/' >>> m1qn3: inconsistent iz for a warm restart') 
  924 format (/' m1qn3: warm restart') 
      iz(1)=n 
      iz(2)=0 
      if (sscale) iz(2)=1 
      iz(3)=m 
c 
c---- split the working zone rz 
c 
      idiag=1 
      iybar=idiag+n 
      if (sscale) iybar=1 
      isbar=iybar+n*mmemo 
      id=isbar+n*mmemo 
      igg=id+n 
      iaux=igg+n 
      ialpha=iaux+n 
c 
c---- call the optimization code 
c 
      call m1qn3a (simul,prosca,ctonb,ctcab,n,x,f,g,dxmin,df1,epsg, 
     /             impres,io,mode,niter,nsim,inmemo,iz(3),iz(4),iz(5), 
     /             rz(id),rz(igg),rz(idiag),rz(iaux),rz(ialpha), 
     /             rz(iybar),rz(isbar),izs,rzs,dzs) 
c 
c---- impressions finales 
c 
      if (impres.ge.1) write (io,905) mode,niter,nsim,epsg 
905   format (/,79('-')/ 
     /        /' m1qn3: output mode is ',i2 
     /        /5x,'number of iterations: ',i4 
     /        /5x,'number of simulations: ',i6 
     /        /5x,'realized relative precision on g: ',e9.2) 
      call prosca (n,x,x,ps,izs,rzs,dzs) 
      r1=sngl(ps) 
      r1=sqrt(r1) 
      call prosca (n,g,g,ps,izs,rzs,dzs) 
      r2=sngl(ps) 
      r2=sqrt(r2) 
      if (impres.ge.1) write (io,906) r1,f,r2 
906   format (5x,'norm of x = ',e15.8 
     /       /5x,'f         = ',e15.8 
     /       /5x,'norm of g = ',e15.8) 
      return 
      end 
c 
c----------------------------------------------------------------------- 
c 
      subroutine m1qn3a (simul,prosca,ctonb,ctcab,n,x,f,g,dxmin,df1, 
     /                   epsg,impres,io,mode,niter,nsim,inmemo,m,jmin, 
     /                   jmax,d,gg,diag,aux,alpha,ybar,sbar,izs,rzs,dzs) 
c---- 
c 
c     Code d'optimisation proprement dit. 
c 
c---- 
c 
c         arguments 
c 
      logical inmemo 
      integer n,impres,io,mode,niter,nsim,m,jmin,jmax,izs(1) 
      real x(n),f,g(n),dxmin,df1,epsg,d(n),gg(n),diag(n),aux(n), 
     /     alpha(m),ybar(n,1),sbar(n,1),rzs(1) 
      double precision dzs(1) 
      external simul,prosca,ctonb,ctcab 
c 
c         variables locales 
c 
      logical sscale,cold,warm 
      integer i,iter,moderl,isim,jcour,indic 
      real r1,t,tmin,tmax,gnorm,eps1,ff,preco,precos,ys,den,dk,dk1 
      double precision ps,ps2,hp0 
c 
c         parametres 
c 
      real rm1,rm2 
      parameter (rm1=0.0001,rm2=0.9) 
      real pi 
      parameter (pi=3.1415927) 
      real rmin 
c 
c---- initialisation 
c 
      rmin=1.e-20 
c 
      sscale=.true. 
      if (mod(mode,2).eq.0) sscale=.false. 
c 
      warm=.false. 
      if (mode/2.eq.1) warm=.true. 
      cold=.not.warm 
c 
      iter=0 
      isim=1 
      eps1=1. 
c 
      call prosca (n,g,g,ps,izs,rzs,dzs) 
      gnorm=sngl(ps) 
      gnorm=sqrt(gnorm) 
      if (impres.ge.1) write (io,900) f,gnorm 
  900 format (5x,'f         = ',e15.8 
     /       /5x,'norm of g = ',e15.8) 
      if (gnorm.lt.rmin) then 
          mode=2 
          if (impres.ge.1) write (io,901) 
          goto 1000 
      endif 
  901 format (/' >>> m1qn3a: initial gradient is too small') 
c 
c     --- initialisation pour dd 
c 
      if (cold) then 
          jmin=1 
          jmax=0 
      endif 
      jcour=1 
      if (inmemo) jcour=jmax 
c 
c     --- mise a l'echelle de la premiere direction de descente 
c 
      if (cold) then 
c 
c         --- use Fletcher's scaling and initialize diag to 1. 
c 
          precos=2.*df1/gnorm**2 
          do 10 i=1,n 
              d(i)=-g(i)*precos 
              diag(i)=1. 
   10     continue 
          if (impres.ge.5) write(io,902) precos 
  902     format (/' m1qn3a: descent direction -g: precon = ',e10.3) 
      else 
c 
c         --- use the matrix stored in [diag and] the (y,s) pairs 
c 
          if (sscale) then 
              call prosca (n,ybar(1,jcour),ybar(1,jcour),ps,izs,rzs,dzs) 
              precos=1./sngl(ps) 
          endif 
          do 11 i=1,n 
              d(i)=-g(i) 
  11      continue 
          if (inmemo) then 
              call dd (prosca,ctonb,ctcab,n,sscale,m,d,aux,jmin,jmax, 
     /                 precos,diag,alpha,ybar,sbar,izs,rzs,dzs) 
          else 
              call dds (prosca,ctonb,ctcab,n,sscale,m,d,aux,jmin,jmax, 
     /                  precos,diag,alpha,ybar,sbar,izs,rzs,dzs) 
          endif 
      endif 
c 
      if (impres.eq.3) then 
          write(io,903) 
          write(io,904) 
      endif 
      if (impres.eq.4) write(io,903) 
  903 format (/,79('-')) 
  904 format (1x) 
c 
c     --- initialisation pour mlis0 
c 
      tmax=1.e+20 
      call prosca (n,d,g,hp0,izs,rzs,dzs) 
      if (hp0.ge.0.d+0) then 
          mode=7 
          if (impres.ge.1) write (io,905) iter,hp0 
          goto 1000 
      endif 
  905 format (/' >>> m1qn3 (iteration ',i2,'): ' 
     /        /5x,' the search direction d is not a ', 
     /         'descent direction: (g,d) = ',d12.5) 
c 
c     --- compute the angle (-g,d) 
c 
      if (warm.and.impres.ge.5) then 
          call prosca (n,g,g,ps,izs,rzs,dzs) 
          ps=dsqrt(ps) 
          call prosca (n,d,d,ps2,izs,rzs,dzs) 
          ps2=dsqrt(ps2) 
          ps=hp0/ps/ps2 
          ps=dmin1(-ps,1.d+0) 
          ps=dacos(ps) 
          r1=sngl(ps) 
          r1=r1*180./pi 
          write (io,906) r1 
      endif 
  906 format (/' m1qn3: descent direction d: ', 
     /        'angle(-g,d) = ',f5.1,' degrees') 
c 
c---- Debut de l'iteration. on cherche x(k+1) de la forme x(k) + t*d, 
c     avec t > 0. On connait d. 
c 
c         Debut de la boucle: etiquette 100, 
c         Sortie de la boucle: goto 1000. 
c 
100   iter=iter+1 
      if (impres.lt.0) then 
          if (mod(iter,-impres).eq.0) then 
              indic=1 
              call simul (indic,n,x,f,g,izs,rzs,dzs) 
              goto 100 
          endif 
      endif 
      if (impres.ge.5) write(io,903) 
      if (impres.ge.4) write(io,904) 
      if (impres.ge.3) write (io,910) iter,isim,f,hp0 
  910 format (' m1qn3: iter ',i3,', simul ',i3, 
     /        ', f=',e15.8,', h(0)=',d12.5) 
      do 101 i=1,n 
          gg(i)=g(i) 
101   continue 
      ff=f 
c 
c     --- recherche lineaire et nouveau point x(k+1) 
c 
      if (impres.ge.5) write (io,911) 
  911 format (/' m1qn3: line search') 
c 
c         --- calcul de tmin 
c 
      tmin=0. 
      do 200 i=1,n 
          tmin=amax1(tmin,abs(d(i))) 
200   continue 
      tmin=dxmin/tmin 
      t=1. 
      r1=sngl(hp0) 
c 
      call mlis0 (n,simul,prosca,x,f,r1,t,tmin,tmax,d,g,rm2,rm1, 
     /           impres,io,moderl,isim,nsim,aux,izs,rzs,dzs) 
c 
c         --- mlis0 renvoie les nouvelles valeurs de x, f et g 
c 
      if (moderl.ne.0) then 
          if (moderl.lt.0) then 
c 
c             --- calcul impossible 
c                 t, g: ou les calculs sont impossibles 
c                 x, f: ceux du t_gauche (donc f <= ff) 
c 
              mode=moderl 
          elseif (moderl.eq.1) then 
c 
c             --- descente bloquee sur tmax 
c                 [sortie rare (!!) d'apres le code de mlis0] 
c 
              mode=3 
              if (impres.ge.1) write(io,912) iter 
  912         format (/' >>> m1qn3 (iteration ',i3, 
     /                '): line search blocked on tmax: ' 
     /                'decrease the scaling') 
          elseif (moderl.eq.4) then 
c 
c             --- nsim atteint 
c                 x, f: ceux du t_gauche (donc f <= ff) 
c 
              mode=5 
          elseif (moderl.eq.5) then 
c 
c             --- arret demande par l'utilisateur (indic = 0) 
c                 x, f: ceux en sortie du simulateur 
c 
              mode=0 
          elseif (moderl.eq.6) then 
c 
c             --- arret sur dxmin ou appel incoherent 
c                 x, f: ceux du t_gauche (donc f <= ff) 
c 
              mode=6 
          endif 
          goto 1000 
      endif 
c 
c NOTE: stopping tests are now done after having updated the matrix, so 
c that update information can be stored in case of a later warm restart 
c 
c     --- mise a jour de la matrice 
c 
      if (m.gt.0) then 
c 
c         --- mise a jour des pointeurs 
c 
          jmax=jmax+1 
          if (jmax.gt.m) jmax=jmax-m 
          if ((cold.and.iter.gt.m).or.(warm.and.jmin.eq.jmax)) then 
              jmin=jmin+1 
              if (jmin.gt.m) jmin=jmin-m 
          endif 
          if (inmemo) jcour=jmax 
c 
c         --- y, s et (y,s) 
c 
          do 400 i=1,n 
              sbar(i,jcour)=t*d(i) 
              ybar(i,jcour)=g(i)-gg(i) 
400       continue 
          if (impres.ge.5) then 
              call prosca (n,sbar(1,jcour),sbar(1,jcour),ps,izs,rzs,dzs) 
              dk1=sqrt(sngl(ps)) 
              if (iter.gt.1) write (io,930) dk1/dk 
  930         format (/' m1qn3: convergence rate, s(k)/s(k-1) = ', 
     /                e12.5) 
              dk=dk1 
          endif 
          call prosca (n,ybar(1,jcour),sbar(1,jcour),ps,izs,rzs,dzs) 
          ys=sngl(ps) 
          if (ys.le.0.) then 
              mode=7 
              if (impres.ge.1) write (io,931) iter,ys 
  931         format (/' >>> m1qn3 (iteration ',i2, 
     /                '): the scalar product (y,s) = ',e12.5 
     /                /27x,'is not positive') 
              goto 1000 
          endif 
c 
c         --- ybar et sbar 
c 
          r1=sqrt(1./ys) 
          do 410 i=1,n 
              sbar(i,jcour)=r1*sbar(i,jcour) 
              ybar(i,jcour)=r1*ybar(i,jcour) 
  410     continue 
          if (.not.inmemo) call ystbl (.true.,ybar,sbar,n,jmax) 
c 
c         --- compute the scalar or diagonal preconditioner 
c 
          if (impres.ge.5) write(io,932) 
  932     format (/' m1qn3: matrix update:') 
c 
c             --- Here is the Oren-Spedicato factor, for scalar scaling 
c 
          if (sscale) then 
              call prosca (n,ybar(1,jcour),ybar(1,jcour),ps,izs,rzs,dzs) 
              precos=1./sngl(ps) 
c 
              if (impres.ge.5) write (io,933) precos 
  933         format (5x,'Oren-Spedicato factor = ',e10.3) 
c 
c             --- Scale the diagonal to Rayleigh's ellipsoid. 
c                 Initially (iter.eq.1) and for a cold start, this is 
c                 equivalent to an Oren-Spedicato scaling of the 
c                 identity matrix. 
c 
          else 
              call ctonb (n,ybar(1,jcour),aux,izs,rzs,dzs) 
              ps=0.d0 
              do 420 i=1,n 
                  ps=ps+dble(diag(i))*dble(aux(i))*dble(aux(i)) 
  420         continue 
              r1=sngl(1.d0/ps) 
              if (impres.ge.5) then 
                  write (io,934) r1 
  934             format(5x,'fitting the ellipsoid: factor = ',e10.3) 
              endif 
              do 421 i=1,n 
                  diag(i)=diag(i)*r1 
  421         continue 
c 
c             --- update the diagonal 
c                 (gg is used as an auxiliary vector) 
c 
              call ctonb (n,sbar(1,jcour),gg,izs,rzs,dzs) 
              ps=0.d0 
              do 430 i=1,n 
                  ps=ps+dble(gg(i))*dble(gg(i))/dble(diag(i)) 
  430         continue 
              den=sngl(ps) 
              do 431 i=1,n 
                  diag(i)=1./ 
     &                    (1./diag(i)+aux(i)**2-(gg(i)/diag(i))**2/den) 
                  if (diag(i).le.0.) then 
                      if (impres.ge.5) write (io,935) i,diag(i),rmin 
                      diag(i)=rmin 
                  endif 
  431         continue 
  935         format (/' >>> m1qn3-WARNING: diagonal element ',i8, 
     &                 ' is negative (',e10.3,'), reset to ',e10.3) 
c 
              if (impres.ge.5) then 
                  ps=0. 
                  do 440 i=1,n 
                      ps=ps+dble(diag(i)) 
  440             continue 
                  ps=ps/n 
                  preco=sngl(ps) 
c 
                  ps2=0. 
                  do 441 i=1,n 
                      ps2=ps2+(dble(diag(i))-ps)**2 
  441             continue 
                  ps2=dsqrt(ps2/n) 
                  write (io,936) preco,sngl(ps2) 
  936             format (5x,'updated diagonal: average value = ',e10.3, 
     &                   ', sqrt(variance) = ',e10.3) 
              endif 
          endif 
      endif 
c 
c     --- tests d'arret 
c 
      call prosca(n,g,g,ps,izs,rzs,dzs) 
      eps1=sngl(ps) 
      eps1=sqrt(eps1)/gnorm 
c 
      if (impres.ge.5) write (io,940) eps1 
  940 format (/' m1qn3: stopping criterion on g: ',e12.5) 
      if (eps1.lt.epsg) then 
          mode=1 
          goto 1000 
      endif 
      if (iter.eq.niter) then 
          mode=4 
          if (impres.ge.1) write (io,941) iter 
  941     format (/' >>> m1qn3 (iteration ',i3, 
     /            '): maximal number of iterations') 
          goto 1000 
      endif 
      if (isim.gt.nsim) then 
          mode=5 
          if (impres.ge.1) write (io,942) iter,isim 
  942     format (/' >>> m1qn3 (iteration ',i3,'): ',i6, 
     /            ' simulations (maximal number reached)') 
          goto 1000 
      endif 
c 
c     --- calcul de la nouvelle direction de descente d = - H.g 
c 
      if (m.eq.0) then 
          preco=2.*(ff-f)/(eps1*gnorm)**2 
          do 500 i=1,n 
              d(i)=-g(i)*preco 
  500     continue 
      else 
          do 510 i=1,n 
              d(i)=-g(i) 
  510     continue 
          if (inmemo) then 
              call dd (prosca,ctonb,ctcab,n,sscale,m,d,aux,jmin,jmax, 
     /                 precos,diag,alpha,ybar,sbar,izs,rzs,dzs) 
          else 
              call dds (prosca,ctonb,ctcab,n,sscale,m,d,aux,jmin,jmax, 
     /                  precos,diag,alpha,ybar,sbar,izs,rzs,dzs) 
          endif 
      endif 
c 
c         --- test: la direction d est-elle de descente ? 
c             hp0 sera utilise par mlis0 
c 
      call prosca (n,d,g,hp0,izs,rzs,dzs) 
      if (hp0.ge.0.d+0) then 
          mode=7 
          if (impres.ge.1) write (io,905) iter,hp0 
          goto 1000 
      endif 
      if (impres.ge.5) then 
          call prosca (n,g,g,ps,izs,rzs,dzs) 
          ps=dsqrt(ps) 
          call prosca (n,d,d,ps2,izs,rzs,dzs) 
          ps2=dsqrt(ps2) 
          ps=hp0/ps/ps2 
          ps=dmin1(-ps,1.d+0) 
          ps=dacos(ps) 
          r1=sngl(ps) 
          r1=r1*180./pi 
          write (io,906) r1 
      endif 
c 
c---- on poursuit les iterations 
c 
      goto 100 
c 
c---- retour 
c 
 1000 niter=iter 
      nsim=isim 
      epsg=eps1 
      return 
      end 
c 
c----------------------------------------------------------------------- 
c 
      subroutine dd (prosca,ctonb,ctcab,n,sscale,nm,depl,aux,jmin,jmax, 
     &               precos,diag,alpha,ybar,sbar,izs,rzs,dzs) 
c---- 
c 
c     calcule le produit H.g ou 
c         . H est une matrice construite par la formule de bfgs inverse 
c           a nm memoires a partir de la matrice diagonale diag 
c           dans un espace hilbertien dont le produit scalaire 
c           est donne par prosca 
c           (cf. J. Nocedal, Math. of Comp. 35/151 (1980) 773-782) 
c         . g est un vecteur de dimension n (en general le gradient) 
c 
c     la matrice diag apparait donc comme un preconditionneur diagonal 
c 
c     depl = g (en entree), = H g (en sortie) 
c 
c     la matrice H est memorisee par les vecteurs des tableaux 
c     ybar, sbar et les pointeurs jmin, jmax 
c 
c     alpha(nm) est une zone de travail 
c 
c     izs(1),rzs(1),dzs(1) sont des zones de travail pour prosca 
c 
c---- 
c 
c         arguments 
c 
      logical sscale 
      integer n,nm,jmin,jmax,izs(1) 
      real depl(n),precos,diag(n),alpha(nm),ybar(n,1),sbar(n,1),rzs(1), 
     &     aux(n) 
      double precision dzs(1) 
      external prosca,ctonb,ctcab 
c 
c         variables locales 
c 
      integer jfin,i,j,jp 
      real r 
      double precision ps 
c 
      jfin=jmax 
      if (jfin.lt.jmin) jfin=jmax+nm 
c 
c         phase de descente 
c 
      do 100 j=jfin,jmin,-1 
          jp=j 
          if (jp.gt.nm) jp=jp-nm 
          call prosca (n,depl,sbar(1,jp),ps,izs,rzs,dzs) 
          r=sngl(ps) 
          alpha(jp)=r 
          do 20 i=1,n 
              depl(i)=depl(i)-r*ybar(i,jp) 
20        continue 
100   continue 
c 
c         preconditionnement 
c 
      if (sscale) then 
          do 150 i=1,n 
              depl(i)=depl(i)*precos 
  150     continue 
      else 
          call ctonb (n,depl,aux,izs,rzs,dzs) 
          do 151 i=1,n 
              aux(i)=aux(i)*diag(i) 
  151     continue 
          call ctcab (n,aux,depl,izs,rzs,dzs) 
      endif 
c 
c         remontee 
c 
      do 200 j=jmin,jfin 
          jp=j 
          if (jp.gt.nm) jp=jp-nm 
          call prosca (n,depl,ybar(1,jp),ps,izs,rzs,dzs) 
          r=alpha(jp)-sngl(ps) 
          do 120 i=1,n 
              depl(i)=depl(i)+r*sbar(i,jp) 
120       continue 
200   continue 
      return 
      end 
c 
c----------------------------------------------------------------------- 
c 
      subroutine dds (prosca,ctonb,ctcab,n,sscale,nm,depl,aux,jmin,jmax, 
     &                precos,diag,alpha,ybar,sbar,izs,rzs,dzs) 
c---- 
c 
c     This subroutine has the same role as dd (computation of the 
c     product H.g). It supposes however that the (y,s) pairs are not 
c     stored in core memory, but on a devise chosen by the user. 
c     The access to this devise is performed via the subroutine ystbl. 
c 
c---- 
c 
c         arguments 
c 
      logical sscale 
      integer n,nm,jmin,jmax,izs(1) 
      real depl(n),precos,diag(n),alpha(nm),ybar(n),sbar(n),rzs(1), 
     &     aux(n) 
      double precision dzs(1) 
      external prosca,ctonb,ctcab 
c 
c         variables locales 
c 
      integer jfin,i,j,jp 
      real r 
      double precision ps 
c 
      jfin=jmax 
      if (jfin.lt.jmin) jfin=jmax+nm 
c 
c         phase de descente 
c 
      do 100 j=jfin,jmin,-1 
          jp=j 
          if (jp.gt.nm) jp=jp-nm 
          call ystbl (.false.,ybar,sbar,n,jp) 
          call prosca (n,depl,sbar,ps,izs,rzs,dzs) 
          r=sngl(ps) 
          alpha(jp)=r 
          do 20 i=1,n 
              depl(i)=depl(i)-r*ybar(i) 
20        continue 
100   continue 
c 
c         preconditionnement 
c 
      if (sscale) then 
          do 150 i=1,n 
              depl(i)=depl(i)*precos 
  150     continue 
      else 
          call ctonb (n,depl,aux,izs,rzs,dzs) 
          do 151 i=1,n 
              aux(i)=aux(i)*diag(i) 
  151     continue 
          call ctcab (n,aux,depl,izs,rzs,dzs) 
      endif 
c 
c         remontee 
c 
      do 200 j=jmin,jfin 
          jp=j 
          if (jp.gt.nm) jp=jp-nm 
          call ystbl (.false.,ybar,sbar,n,jp) 
          call prosca (n,depl,ybar(1),ps,izs,rzs,dzs) 
          r=alpha(jp)-sngl(ps) 
          do 120 i=1,n 
              depl(i)=depl(i)+r*sbar(i) 
120       continue 
200   continue 
      return 
      end 
c 
c----------------------------------------------------------------------- 
c 
      subroutine mlis0 (n,simul,prosca,xn,fn,fpn,t,tmin,tmax,d,g, 
     1                  amd,amf,imp,io,logic,nap,napmax,x,izs,rzs,dzs) 
c ---- 
c 
c     mlis0 + minuscules + commentaires 
c     + version amelioree (XII 88): interpolation cubique systematique 
c       et anti-overflows 
c     + declaration variables (II/89, JCG). 
c 
c     ---------------------------------------------------------------- 
c 
c        en sortie logic = 
c 
c        0          descente serieuse 
c        1          descente bloquee 
c        4          nap > napmax 
c        5          retour a l'utilisateur 
c        6          fonction et gradient pas d'accord 
c        < 0        contrainte implicite active 
c 
c ---- 
c 
c --- arguments 
c 
      external simul,prosca 
      integer n,imp,io,logic,nap,napmax,izs(*) 
      real xn(n),fn,fpn,t,tmin,tmax,d(n),g(n),amd,amf,x(n),rzs(*) 
      double precision dzs(*) 
c 
c --- variables locales 
c 
      integer i,indic,indica,indicd 
      real tesf,tesd,tg,fg,fpg,td,ta,fa,fpa,d2,f,fp,ffn,fd,fpd, 
     1 z,test,barmin,barmul,barmax,barr,gauche,droite, 
     2 taa 
      double precision ps 
c 
 1000 format (/4x,9h mlis0   ,4x,4hfpn=,e10.3,4h d2=,e9.2, 
     1 7h  tmin=,e9.2,6h tmax=,e9.2) 
 1001 format (/4x,6h mlis0,3x,'stop on tmin',8x, 
     1   'step',11x,'functions',5x,'derivatives') 
 1002 format (4x,6h mlis0,37x,e10.3,2e11.3) 
 1003 format (4x,6h mlis0,e14.3,2e11.3) 
 1004 format (4x,6h mlis0,37x,e10.3,7h indic=,i3) 
 1005 format (4x,6h mlis0,14x,2e18.8,e11.3) 
 1006 format (4x,6h mlis0,14x,e18.8,12h      indic=,i3) 
 1007 format (/4x,6h mlis0,10x,'tmin forced to tmax') 
 1008 format (/4x,6h mlis0,10x,'inconsistent call') 
      if (n.gt.0 .and. fpn.lt.0. .and. t.gt.0. 
     1 .and. tmax.gt.0. .and. amf.gt.0. 
     1 .and. amd.gt.amf .and. amd.lt.1.) go to 5 
      logic=6 
      go to 999 
    5 tesf=amf*fpn 
      tesd=amd*fpn 
      barmin=0.01 
      barmul=3. 
      barmax=0.3 
      barr=barmin 
      td=0. 
      tg=0. 
      fg=fn 
      fpg=fpn 
      ta=0. 
      fa=fn 
      fpa=fpn 
      call prosca (n,d,d,ps,izs,rzs,dzs) 
      d2=sngl(ps) 
c 
c               elimination d'un t initial ridiculement petit 
c 
      if (t.gt.tmin) go to 20 
      t=tmin 
      if (t.le.tmax) go to 20 
      if (imp.gt.0) write (io,1007) 
      tmin=tmax 
   20 if (fn+t*fpn.lt.fn+0.9*t*fpn) go to 30 
      t=2.*t 
      go to 20 
   30 indica=1 
      logic=0 
      if (t.gt.tmax) then 
          t=tmax 
          logic=1 
      endif 
      if (imp.ge.4) write (io,1000) fpn,d2,tmin,tmax 
c 
c     --- nouveau x 
c 
      do 50 i=1,n 
          x(i)=xn(i)+t*d(i) 
   50 continue 
c 
c --- boucle 
c 
  100 nap=nap+1 
      if(nap.gt.napmax) then 
          logic=4 
          fn=fg 
          do 120 i=1,n 
              xn(i)=xn(i)+tg*d(i) 
  120     continue 
          go to 999 
      endif 
      indic=4 
c 
c     --- appel simulateur 
c 
      call simul(indic,n,x,f,g,izs,rzs,dzs) 
      if(indic.eq.0) then 
c 
c         --- arret demande par l'utilisateur 
c 
          logic=5 
          fn=f 
          do 170 i=1,n 
              xn(i)=x(i) 
  170     continue 
          go to 999 
      endif 
      if(indic.lt.0) then 
c 
c         --- les calculs n'ont pas pu etre effectues par le simulateur 
c 
          td=t 
          indicd=indic 
          logic=0 
          if (imp.ge.4) write (io,1004) t,indic 
          t=tg+0.1*(td-tg) 
          go to 905 
      endif 
c 
c     --- les tests elementaires sont faits, on y va 
c 
      call prosca (n,d,g,ps,izs,rzs,dzs) 
      fp=sngl(ps) 
c 
c     --- premier test de Wolfe 
c 
      ffn=f-fn 
      if(ffn.gt.t*tesf) then 
          td=t 
          fd=f 
          fpd=fp 
          indicd=indic 
          logic=0 
          if(imp.ge.4) write (io,1002) t,ffn,fp 
          go to 500 
      endif 
c 
c     --- test 1 ok, donc deuxieme test de Wolfe 
c 
      if(imp.ge.4) write (io,1003) t,ffn,fp 
      if(fp.gt.tesd) then 
          logic=0 
          go to 320 
      endif 
      if (logic.eq.0) go to 350 
c 
c     --- test 2 ok, donc pas serieux, on sort 
c 
  320 fn=f 
      do 330 i=1,n 
          xn(i)=x(i) 
  330 continue 
      go to 999 
c 
c 
c 
  350 tg=t 
      fg=f 
      fpg=fp 
      if(td.ne.0.) go to 500 
c 
c              extrapolation 
c 
      taa=t 
      gauche=(1.+barmin)*t 
      droite=10.*t 
      call ecube (t,f,fp,ta,fa,fpa,gauche,droite) 
      ta=taa 
      if(t.lt.tmax) go to 900 
      logic=1 
      t=tmax 
      go to 900 
c 
c              interpolation 
c 
  500 if(indica.le.0) then 
          ta=t 
          t=0.9*tg+0.1*td 
          go to 900 
      endif 
      test=barr*(td-tg) 
      gauche=tg+test 
      droite=td-test 
      taa=t 
      call ecube (t,f,fp,ta,fa,fpa,gauche,droite) 
      ta=taa 
      if (t.gt.gauche .and. t.lt.droite) then 
          barr=barmin 
        else 
          barr=amin1(barmul*barr,barmax) 
      endif 
c 
c --- fin de boucle 
c     - t peut etre bloque sur tmax 
c       (venant de l'extrapolation avec logic=1) 
c 
  900 fa=f 
      fpa=fp 
  905 indica=indic 
c 
c --- faut-il continuer ? 
c 
      if (td.eq.0.) go to 950 
      if (td-tg.lt.tmin) go to 920 
c 
c     --- limite de precision machine (arret de secours) ? 
c 
      do 910 i=1,n 
          z=xn(i)+t*d(i) 
          if (z.ne.xn(i).and.z.ne.x(i)) go to 950 
  910 continue 
c 
c --- arret sur dxmin ou de secours 
c 
  920 logic=6 
c 
c     si indicd<0, derniers calculs non faits par simul 
c 
      if (indicd.lt.0) logic=indicd 
c 
c     si tg=0, xn = xn_depart, 
c     sinon on prend xn=x_gauche qui fait decroitre f 
c 
      if (tg.eq.0.) go to 940 
      fn=fg 
      do 930 i=1,n 
  930 xn(i)=xn(i)+tg*d(i) 
  940 if (imp.le.0) go to 999 
      write (io,1001) 
      write (io,1005) tg,fg,fpg 
      if (logic.eq.6) write (io,1005) td,fd,fpd 
      if (logic.eq.7) write (io,1006) td,indicd 
      go to 999 
c 
c               recopiage de x et boucle 
c 
  950 do 960 i=1,n 
  960 x(i)=xn(i)+t*d(i) 
      go to 100 
  999 return 
      end 
c 
c----------------------------------------------------------------------- 
c 
      subroutine ecube(t,f,fp,ta,fa,fpa,tlower,tupper) 
c 
c --- arguments 
c 
      real t,f,fp,ta,fa,fpa,tlower,tupper 
c 
c --- variables locales 
c 
      real sign,den,anum 
      double precision z1,b,discri 
c 
c           Using f and fp at t and ta, computes new t by cubic formula 
c           safeguarded inside [tlower,tupper]. 
c 
      z1=dble(fp)+dble(fpa)-3.d0*dble(fa-f)/dble(ta-t) 
      b=z1+dble(fp) 
c 
c              first compute the discriminant (without overflow) 
c 
      if (abs(z1).le.1.) then 
          discri=z1*z1-dble(fp)*dble(fpa) 
        else 
          discri=dble(fp)/z1 
          discri=discri*dble(fpa) 
          discri=z1-discri 
          if (z1.ge.0.d0 .and. discri.ge.0.d0) then 
              discri=dsqrt(z1)*dsqrt(discri) 
              go to 120 
          endif 
          if (z1.le.0.d0 .and. discri.le.0.d0) then 
              discri=dsqrt(-z1)*dsqrt(-discri) 
              go to 120 
          endif 
          discri=-1.d0 
      endif 
      if (discri.lt.0.d0) then 
          if (fp.lt.0.) t=tupper 
          if (fp.ge.0.) t=tlower 
          go to 900 
      endif 
c 
c  discriminant nonnegative, compute solution (without overflow) 
c 
      discri=dsqrt(discri) 
  120 if (t-ta.lt.0.) discri=-discri 
      sign=(t-ta)/abs(t-ta) 
      if (sngl(b)*sign.gt.0.) then 
          t=t+fp*(ta-t)/sngl(b+discri) 
        else 
          den=z1+b+dble(fpa) 
          anum=b-discri 
          if (abs((t-ta)*anum).lt.(tupper-tlower)*abs(den)) then 
              t=t+anum*(ta-t)/den 
            else 
              t=tupper 
          endif 
      endif 
  900 t=amax1(t,tlower) 
      t=amin1(t,tupper) 
      return 
      end 
c 
c----------------------------------------------------------------------- 
c 
      subroutine mupdts (sscale,inmemo,n,m,nrz) 
c 
c         arguments 
c 
      logical sscale,inmemo 
      integer n,m,nrz 
c---- 
c 
c     On entry: 
c       sscale: .true. if scalar initial scaling, 
c               .false. if diagonal initial scaling 
c       n:      number of variables 
c 
c     This routine has to return: 
c       m:      the number of updates to form the approximate Hessien H, 
c       inmemo: .true., if the vectors y and s used to form H are stored 
c                  in core memory, 
c               .false. otherwise (storage of y and s on disk, for 
c                  instance). 
c     When inmemo=.false., the routine `ystbl', which stores and 
c     restores (y,s) pairs, has to be rewritten. 
c 
c---- 
c 
      if (sscale) then 
          m=(nrz-3*n)/(2*n+1) 
      else 
          m=(nrz-4*n)/(2*n+1) 
      endif 
      inmemo=.true. 
      return 
      end 
c 
c----------------------------------------------------------------------- 
c 
      subroutine ystbl (store,ybar,sbar,n,j) 
c---- 
c 
c     This subroutine should store (if store = .true.) or restore 
c     (if store = .false.) a pair (ybar,sbar) at or from position 
c     j in memory. Be sure to have 1 <= j <= m, where m in the number 
c     of updates specified by subroutine mupdts. 
c 
c     The subroutine is used only when the (y,s) pairs are not 
c     stored in core memory in the arrays ybar(.,.) and sbar(.,.). 
c     In this case, the subroutine has to be written by the user. 
c 
c---- 
c 
c         arguments 
c 
      logical store 
      integer n,j 
      real ybar(n),sbar(n) 
c 
      return 
      end 
c 
c----------------------------------------------------------------------- 
c 
      subroutine ctonbe (n,u,v,izs,rzs,dzs) 
      integer n,izs(1) 
      real u(1),v(1),rzs(1) 
      double precision dzs(1) 
c 
      integer i 
c 
      do 1 i=1,n 
          v(i)=u(i) 
 1    continue 
      return 
      end 
c 
c----------------------------------------------------------------------- 
c 
      subroutine ctcabe (n,u,v,izs,rzs,dzs) 
      integer n,izs(1) 
      real u(1),v(1),rzs(1) 
      double precision dzs(1) 
c 
      integer i 
c 
      do 1 i=1,n 
          v(i)=u(i) 
 1    continue 
      return 
      end 
c 
c----------------------------------------------------------------------- 
c 
      subroutine euclid (n,x,y,ps,izs,rzs,dzs) 
      integer n,izs(1) 
      real x(1),y(1),rzs(1) 
      double precision ps,dzs(1) 
c 
      integer i 
c 
      ps=0.d0 
      do 10 i=1,n 
   10 ps=ps+dble(x(i))*dble(y(i)) 
      return 
      end 

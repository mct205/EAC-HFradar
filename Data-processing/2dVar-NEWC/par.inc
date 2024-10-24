       PARAMETER (LX=71,LY=71,LX1=LX-1,LY1=LY-1,N=3502,NB=3502) ! 
       PARAMETER (nradars=2,nvec=2*n)
       common/b/x(nb),y(nb),wei(nb,4),cs(nb),si(nb),ix(nb),iy(nb)
       common /wei/ wei_harm(lx,ly),wei_dat(nb),void,iharm,iprint
       common/wei_di/ wdiv(lx,ly),wrot(lx,ly),WUV
       common/mask/ msc(lx,ly),msk(lx,ly),msh(lx,ly),udata(nb),nb7
       common/rigid/ in(lx),iyb(lx),bb(lx),psic,psig
       
 
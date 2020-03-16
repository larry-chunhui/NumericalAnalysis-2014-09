!!!!!!!!!!!!!!!!!!!!!!!!!
!hw 2.12.9
!��xyΪ����ֵ���������Է�����v,w,t,u
! ţ�ٷ��������Է�����
!�������շ�Ƭ���β�ֵ
!��С���˷��������
!!!!!!!!!!!!!!!!!!!!!!!!!

      program step
	implicit none
!      integer,parameter::k=5
      integer,parameter::n=6
      integer,parameter::pe=11
      integer,parameter::qe=21
	integer,parameter::it=1000  !��������������
	real,parameter::er=1.e-12  !��������������
	real x(pe),y(qe),u(pe,qe),t(pe,qe),z(pe,qe)  !������Է������õ��ı���
	real ub(n),tb(n),zb(n,n)
	real c(20,20)  !���ϵ��,��϶���ʽ��ֵ
	real v,w    
	real error                  !�в�
	integer i,j,m,l,k   !���ɱ���
	do i=1,pe
	x(i)=0.08*(i-1)
	end do

	do j=1,qe
	y(j)=0.5+0.05*(j-1)
	end do

	call solvaleq(x,y,u,t,pe,qe)    !��11*21�������Է��̲��ش���֤
	call output(u,t,z,pe,qe)        !���utz��ֵ

	call fenpianchazhi(u,t,z,pe,qe,n)  !��Ƭ���β�ֵ���ut��Ӧ��z
	call outputxyz(x,y,z,pe,qe)

	call solvek(c,x,y,z,pe,qe,k)
!	call lastquestion(c,k)
!	call nihe(c,x,y,z,pe,qe,k)
!	call totalerror(error,c,x,y,z,k,pe,qe)
!	write(*,*)error
	stop 
	end 

	subroutine lastquestion (c,k)
	integer,parameter::n=6
      integer,parameter::pe=8 !�涨x geshu
      integer,parameter::qe=5 !�涨y geshu
	integer i,j
	real x(pe),y(qe),z(pe,qe),u(pe,qe),t(pe,qe),p(pe,qe)
	real ub(n),tb(n),zb(n,n)
	real c(20,20)
	
	open(100,FILE='last xyzp.txt',STATUS='UNKNOWN')
	do i=1,pe
	x(i)=0.1*i
	end do
	do j=1,qe
	y(j)=0.5+0.2*j
	end do
!	write(*,*)(x(i),i=1,pe)
!������Է�����
	call solvaleq(x,y,u,t,pe,qe)    !��8*5�������Է��̲��ش���֤
	call output(u,t,z,pe,qe)        !���utz��ֵ

	call fenpianchazhi(u,t,z,pe,qe,n)  !��Ƭ���β�ֵ���ut��Ӧ��z
	call outputxyz(x,y,z,pe,qe)        
	do i=1,pe                            !����֮ǰ�������Cij,��p��ֵͬz�Ա�
	do j=1,qe
	call fpvalue(c,x(i),y(j),k,fp)
	p(i,j)=fp
	write(100,10)x(i),y(j),z(i,j),p(i,j)
	end do
	end do
10    FORMAT (F4.1,F4.1,E20.12,E20.12)
	close(100)




	end 

	subroutine solvek(c,x,y,z,p,q,k)
	integer i,j
	integer p,q,k
	real x(p),y(q),z(p,q)
      integer,parameter::m=20 !�涨����������
	real c(m,m)
	do k=1,19
	call nihe(c,x,y,z,p,q,k)
	call totalerror(error,c,x,y,z,k,p,q)
!���Բ�ͬ��k����Ͼ��ȵ�Ӱ��
!	if (error.lt.(1.e-7)) then
!	write(*,*)k,error
!	exit
!	else 
	write(*,*) k,error
!	end if

	end do
	open(6,FILE='k error c.txt',STATUS='UNKNOWN')
	write(6,*)'k=',k
	write(6,*)'error=',error
	write(6,*)'r=      s=       Crs='
	do i=1,k+1
	do j=1,k+1
	write(6,*) i,j,c(i,j)
	end do
	end do

	end 

	subroutine totalerror(err,c,x,y,z,k,p,q)
	integer p,q
	real c(20,20),x(p),y(q),z(p,q),fp
	integer i,j
	real err
	err=0
	do i=1,p
	do j=1,q
	call fpvalue(c,x(i),y(j),k,fp)
	err=err+(z(i,j)-fp)**2
	end do
	end do

	end 


	subroutine fpvalue(c,x,y,k,fp)
	real c(20,20),x,y
	integer i,j
	real s
	s=0
	do i=1,k+1
	do j=1,k+1
	s=s+c(i,j)*x**(i-1)*y**(j-1)
	end do
	end do
	fp=s
	end

	subroutine nihe(c,x,y,z,p,q,k)
	integer p,q,k
	integer i,j,l
	real c(20,20),x(p),y(q),z(p,q)
	real b(p,k+1),b1(k+1,k+1),b2(k+1,q),a(k+1,q)
	real g(q,k+1),g1(k+1,k+1),d(q,k+1),e(k+1)
	real s

	do i=1,p
	do j=1,k+1
	b(i,j)=x(i)**(j-1)
	end do
	end do
!��b1
      do i=1,k+1
	do j=1,k+1
	s=0
	do l=1,p
	s=s+b(l,i)*b(l,j)
	end do
	b1(i,j)=s
	end do
	end do
!��b1���棬���ظ�b1
      do i=1,k+1
	e(i)=1
	end do
	call gaussj(b1,k+1,e)
!��b2  B2=Bt*z
      do i=1,k+1
	do j=1,q
	s=0
	do l=1,p
	s=s+b(l,i)*z(l,j)
	end do
	b2(i,j)=s
	end do
	end do
!��A   A=B1*B2
      do i=1,k+1
	do j=1,q
	s=0
	do l=1,k+1
	s=s+b1(i,l)*b2(l,j)
	end do
	a(i,j)=s
	end do
	end do
!X������Ͻ�������ʼY�������
!��g
      do i=1,q
	do j=1,k+1
	g(i,j)=y(i)**(j-1)
	end do
	end do
!��g1 G1=Gt*G
      do i=1,k+1
	do j=1,k+1
	s=0
	do l=1,q
	s=s+g(l,i)*g(l,j)
	end do
	g1(i,j)=s
	end do
	end do
!��G1���棬���ظ�G1
      CALL gaussj(g1,k+1,e)
!��D  D=G*G1
      do i=1,q
	do j=1,k+1
	s=0
	do l=1,k+1
	s=s+g(i,l)*g1(l,j)
	end do
	d(i,j)=s
	end do
	end do
!��c   C=A*D
      DO i=1,k+1
	do j=1,k+1
	s=0
	do l=1,q
	s=s+a(i,l)*d(l,j)
	end do
	c(i,j)=s
	end do
	end do

!check
      do i=1,k+1
	do j=1,k+1
!	write(*,*)i,j,c(i,j)
	end do
	end do
      
	end 


	subroutine outputxyz(x,y,z,p,q)
	integer i,j,p,q
	real x(p),y(q),z(p,q)
	open(4,FILE='xyz.txt',STATUS='UNKNOWN')
	write(4,*) 'x           y             z'
	do i=1,p
	do j=1,q
	write(4,*) x(i),y(j),z(i,j)
	end do
	end do

	end 
	
	subroutine fenpianchazhi(u,t,z,p,q,n)  !��Ƭ���β�ֵ���ut��Ӧ��z
	integer p,q
	real r(4)
	real u(p,q),t(p,q),z(p,q)
	real ub(n),tb(n),zb(n,n)
	integer i,j,k,l,m,s
!	write(*,*) p,q,n

	call readdata(ub,tb,zb,n)      !�������е�����

	do i=1,p             !p
	do j=1,q           !q

	do k=1,n-1                     !�ҵ�u�������䣬��¼k
	if (u(i,j).GE.ub(k).AND.u(i,j).LT.ub(k+1)) then
	exit
	end if
	end do

	do l=1,n-1                      !�ҵ�t�������䣬��¼l
	if (t(i,j).GE.tb(l).AND.t(i,j).LT.tb(l+1)) then
	exit
	end if
	end do
      
!�����������Ľڵ���Ϊ��ֵ���Ľڵ�
	call distance(u(i,j),t(i,j),ub(k),tb(l),r(1))
	call distance(u(i,j),t(i,j),ub(k),tb(l+1),r(2))
	call distance(u(i,j),t(i,j),ub(k+1),tb(l),r(3))
	call distance(u(i,j),t(i,j),ub(k+1),tb(l+1),r(4))
	call nearestp(r,s)
	if (s.eq.1) then
	k=k
	l=l
	end if
	if(s.eq.2) then
	k=k
	l=l+1
	end if
	if(s.eq.3) then
	k=k+1
	l=l
	end if
	if(s.eq.4) then
	k=k+1
	l=l+1
	end if
!�������Ľڵ�λ�ڱ߽����������ڲ�ƽ��һ��
      if (k.eq.1) then
	    k=k+1
	end if
      if(k.eq.n) then
	    k=k-1
	end if

	if (l.eq.1) then
	    l=l+1
	end if
      if(l.eq.n) then
	    l=l-1
	end if
!������Ľڵ����½ǽڵ���k��l
	k=k-1
	l=l-1
!	write(*,*),k,l
	call interpolate(ub,tb,zb,u(i,j),t(i,j),z(i,j),k,l,n)
!	write(*,*)'z(i,j)',i,j,z(i,j)        !�����������z����0

	end do
	end do
	end
!��Ƭ���β�ֵ
	subroutine interpolate(ub,tb,zb,u,t,z,k,l,n)
	real ub(n),tb(n),zb(n,n)
	real u,t,z
	integer k,l
	dimension x1a(3),x2a(3),ya(3,3)
	real x1,x2,y,dy
	integer i,j
	integer m
	m=3
	x1=u
	x2=t

      do i=1,3
	  x1a(i)=ub(k+i-1)
	do j=1,3
	  x2a(j)=tb(l+j-1)
	  ya(i,j)=zb(k+i-1,l+j-1)
	end do
	end do
      
!	write(*,*) x1,x2
!	write(*,*) (x1a(i),i=1,3)
!	write(*,*) (x2a(i),i=1,3)
	do i=1,3
	write(*,*)(ya(i,j),j=1,3)
	end do
	call POLIN2(X1A,X2A,YA,m,m,X1,X2,Y,DY)
	z=y
!	write(*,*)'DY=',dy
!	write(*,*)'Y=',y,'z=',z
	end 
!��Ƭ���β�ֵ�����㷨 �������ղ�ֵ

       SUBROUTINE polin2(x1a,x2a,ya,m,n,x1,x2,y,dy)
       INTEGER m,n,NMAX,MMAX
       REAL dy,x1,x2,y,x1a(m),x2a(n),ya(m,n)
        PARAMETER (NMAX=20,MMAX=20)
        !USES polint
       INTEGER j,k
       REAL ymtmp(MMAX),yntmp(NMAX)
       do j=1,m
       do k=1,n
       yntmp(k)=ya(j,k)
        end do
        call polint(x2a,yntmp,n,x2,ymtmp(j),dy)
        end do
        call polint(x1a,ymtmp,m,x1,y,dy)
       END SUBROUTINE polin2

       SUBROUTINE polint(xa,ya,n,x,y,dy)
          INTEGER n,NMAX
        REAL dy,x,y,xa(n),ya(n)
       PARAMETER (NMAX=10)
       INTEGER i,m,ns
      REAL den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
        ns=1
      dif=abs(x-xa(1))
      do i=1,n
      dift=abs(x-xa(i))
       if (dift<dif) then
       ns=i
       dif=dift
      endif
      c(i)=ya(i)
      d(i)=ya(i)
       end do
      y=ya(ns)
       ns=ns-1
       do m=1,n-1
      do i=1,n-m
      ho=xa(i)-x
      hp=xa(i+m)-x
      w=c(i+1)-d(i)
      den=ho-hp
      if(den==0.) pause 'failure in polint'
      den=w/den
      d(i)=hp*den
      c(i)=ho*den
      end do
      if (2*ns<n-m) then
      dy=c(ns+1)
      else
      dy=d(ns)
      ns=ns-1
      endif
      y=y+dy
      end do
      END SUBROUTINE polint


!�ҵ������ֵ�ڵ�����Ľڵ㣬����p
      subroutine nearestp(r,s)
	integer i,s
	real r(4),t
	s=1
      t=r(1)

	do i=2,4
	if(t.GT.r(i)) then
	t=r(i)
	s=i
	end if
	end do

	
	end


!������֮��ľ���,�����
	subroutine distance(x1,y1,x2,y2,d)
	real x1,y1,x2,y2,d
	d=sqrt((x1-x2)**2+(y1-y2)**2)
	end

	subroutine checkdistance()
	real x1,y1,x2,y2,r
	x1=1
	x2=2
	y1=1
	y2=2
	call distance(x1,y1,x2,y2,r)
!	write(*,*)'-------',r
	end

!�������ݿα��еı��A1
	subroutine readdata(u,t,z,n)
	real u(n),t(n),z(n,n)
	integer i,j,k
	character*65 DUMMY
      open(1,FILE='shuju.txt',STATUS='UNKNOWN')

	read(1,300) DUMMY
	write(*,300)DUMMY

	read(1,*) (u(i),i=1,6)
	WRITE(*,*) (u(i),i=1,6)

	read(1,*) DUMMY
	write(*,*)DUMMY

	read(1,*)  (t(i),i=1,6)
	write(*,*) (t(i),i=1,6)
	read(1,300) DUMMY
	write(*,300)DUMMY

	do j=1,n

	read(1,*) (z(i,j),i=1,6)

	end do
	close(1)


!      read(1,*) z(1,:),z(2,:),z(3,:),z(4,:),z(5,:),z(6,:)

      write(*,*) z(1,:),z(2,:),z(3,:),z(4,:),z(5,:),z(6,:)


300   FORMAT(A65)

	end 
	subroutine output(u,t,z,p,q)
	integer p,q
	integer i,j
	real u(p,q),t(p,q),z(p,q)
	open(2,FILE='utz.txt',STATUS='UNKNOWN')
	write(2,*) 'u     t      z'
	do i=1,p
	do j=1,q
	write(2,*) u(i,j),t(i,j),z(i,j)
	end do
	end do

	end 

	subroutine 	solvaleq(x,y,u,t,p,q)
	integer p,q
	real x(p),y(q),u(p,q),t(p,q)
	real v,w
	integer i,j
	
	open(3,FILE='solerror.txt',STATUS='UNKNOWN')


	do i=1,p
	do j=1,q
	call solnleq(v,w,x(i),y(j),u(i,j),t(i,j))

	write(3,*) u(i,j),t(i,j)
	write(3,*) 0.5*cos(t(i,j))+u(i,j)+v+w-x(i)-2.67
	write(3,*) t(i,j)+0.5*sin(u(i,j))+v+w-y(j)-1.07
	write(3,*) 0.5*t(i,j)+u(i,j)+cos(v)+w-x(i)-3.74
	write(3,*) t(i,j)+0.5*u(i,j)+v+sin(w)-y(j)-0.79
	end do
	end do


	end



      subroutine solnleq(v,w,x,y,u,t)
	real v,w,x,y,u,t
	real fanshu(2)
	real d(4,100000),deltd(4)
	integer i,k,j
	real er
	real error
!12.9���
	integer,parameter::m=4
	real A(m,m),s(m),b(m)  !�ⷽ��As=b
	er=1.E-12
	v=1
	w=1
	u=1
	t=0


	do i=2,100000
	
	do k=1,m
	do j=1,m
	A(k,j)=1
	end do 
	end do
	A(3,4)=0.5
	A(4,3)=0.5
	A(1,4)=-0.5*sin(t)
	A(2,3)=0.5*cos(u)
	A(3,1)=-sin(v)
	A(4,2)=cos(w)

	b(1)=-(v+w+u+0.5*cos(t)-(2.67+x))
	b(2)=-(v+w+0.5*sin(u)+t-(y+1.07))
	b(3)=-(cos(v)+w+u+0.5*t-(x+3.74))
	b(4)=-(v+sin(w)+0.5*u+t-(y+0.79))
	
	call soleq(m,A,b,s)
!	write(*,*)('--------')
!	write(*,*)(s(j),j=1,m)
!	write(*,*)('--------')


	v=v+s(1)
	w=w+s(2)
	u=u+s(3)
	t=t+s(4)
	
	error=max(s(1),s(2),s(3),s(4),fanshu(1))/
     >            max(v,w,u,t,fanshu(2))

	 
	if(error.LT.er) then
	exit
	end if
      
	if (i.eq.100000.and.error.GT.er) then
	write(*,*)'netown failed i=',i,error
	end if

	end do
	END

	subroutine max(a,b,c,d,e)
	real a,b,c,d,e
	real t,u(4)
	integer i

	u(1)=abs(a)
	u(2)=abs(b)
	u(3)=abs(c)
	u(4)=abs(d)
	t=u(1)

	do i=2,4
	if(t.LT.u(i)) then
	t=u(i)
	end if
	end do
	e=t

	end 

	subroutine soleq(m,c,b,x)
	real c(m,m),b(m),x(m)
	integer i
	CALL gaussj(c,m,b)
	do i=1,m
	x(i)=b(i)
	end do
	end 
		

      
      SUBROUTINE gaussj(a,n,b)
       INTEGER n,NMAX
      REAL a(n,n),b(n)
       PARAMETER (NMAX=50)
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX)
      INTEGER  indxr(NMAX),ipiv(NMAX)
      REAL big,dum,pivinv
      do j=1,n
         ipiv(j)=0
      end do
       do i=1,n
         big=0.
      do j=1,n
       if(ipiv(j)/=1) then
      do k=1,n
        if (ipiv(k)==0) then
          if (abs(a(j,k))>=big) then
            big=abs(a(j,k))
            irow=j
            icol=k
          endif
        else if (ipiv(k)>1) then
          pause 'singular matrix in gaussj'
        endif
      end do
      endif
      end do
      ipiv(icol)=ipiv(icol)+1
      if (irow/=icol) then
      do l=1,n
      dum=a(irow,l)
      a(irow,l)=a(icol,l)
      a(icol,l)=dum
      end do
      dum=b(irow)
	b(irow)=b(icol)
	b(icol)=dum
      endif
      indxr(i)=irow
      indxc(i)=icol
      if (a(icol,icol)==0.) pause 'singular matrix in gaussj'
      pivinv=1./a(icol,icol)
      a(icol,icol)=1.
      do l=1,n
      a(icol,l)=a(icol,l)*pivinv
      end do
      b(icol)=b(icol)*pivinv
      do ll=1,n
      if(ll/=icol) then
      dum=a(ll,icol)
      a(ll,icol)=0.
      do l=1,n
        a(ll,l)=a(ll,l)-a(icol,l)*dum
      end do
      b(ll)=b(ll)-b(icol)*dum
       endif
      end do
      end do
      do l=n,1,-1
      if(indxr(l)/=indxc(l)) then
      do k=1,n
      dum=a(k,indxr(l))
      a(k,indxr(l))=a(k,indxc(l))
      a(k,indxc(l))=dum
      end do
      endif
      end do
       END SUBROUTINE gaussj





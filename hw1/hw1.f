!!!!!!!!!!!!!!!!!!!!!!!!!
! 1.1幂法求按模最大特征值c(m,m)
! 1.2LU Doolittle solove eq c(m,m)*x(m)=b(m)
! E(M,M)为单位矩阵
!矩阵0元素不存储第一问求lamb(1)，lamb(501)正确，LUDOOLITTLE(m,c,lu)修改正确，反幂法中逆矩阵不能表示为A的形式，求不下去了11.15
!反幂法有问题，不用去求A 的逆矩阵，参考课本对反幂法算法说明11.16
!经测试修改后的LU分解正确
!经测试反幂法正确，1.1求解完毕
!全部搞定 detAproblem
!!!!!!!!!!!!!!!!!!!!!!!!!
      program step
	integer,parameter::m=501
	real c(m+2),ff(m+2),cre(m,m),d(m+2),lu(5,m),b(m),x(m)
	real u(39),condA2  !question 1.2
	real ev,minabsev,maxabsev
	real lamb(m)       !question 1.1
	integer i,j,iter
	OPEN(1,FILE='shufen hw.1.txt',STATUS='UNKNOWN')

	do i=1,m
	c(i)=(1.64-0.024*i)*sin(0.2*i)-0.64*exp(0.1/i)
!      c(i)=i
	end do
c	
	c(m+1)=0.16
	c(m+2)=-0.064

! question 1.1	
	CALL eighenvalue(m,c,ev)
	lamb(1)=ev
	maxabsev=ev           !A的按模最大特征值
c	write(*,*)ev

      do i=1,m
	d(i)=c(i)-lamb(1)
	end do 
      d(m+1)=c(m+1)
	d(m+2)=c(m+2)
	CALL eighenvalue(m,d,ev)
c	write(*,*)ev
	lamb(m)=ev+lamb(1)
c	write(*,*)lamb(m)

	if(lamb(m).LT.lamb(1)) then 
	ev=lamb(m)
	lamb(m)=lamb(1)
	lamb(1)=ev
	end if 
	write(*,*)'lamb(1)=,lamb(m)=',lamb(1),lamb(m)
	write(1,*)'lamb(1)=',lamb(1),  'lamb(m)=',lamb(m)
	CALL  reeighenvalue(m,c,ev,iter)
	minabsev=ev          !A的按模最小特征值
	write(*,*)'minabsev=',1/ev
	write(1,*)'minabsev=',1/ev
!question 1.2 
      write(*,*)'question1.2-----------------'
	write(1,*)'question1.2-----------------'
	write(1,*)'iter times=           ', 'eigenvalue='
	ff(m+1)=c(m+1)
	ff(m+2)=c(m+2)
	do k=1,39
	u(k)=lamb(1)+k*(lamb(m)-lamb(1))/40

	do i=1,m
	ff(i)=c(i)-u(k)
	end do
	CALL  reeighenvalue(m,ff,ev,iter)
	write(*,*)1/ev+u(k)
	write(1,*)iter,1/ev+u(k)
	end do
!question 1.3
      write(1,*)'question1.3-----------------'
	CALL eighenvalue(m,c,ev)
	maxabsev=ev 
	WRITE(*,*)EV
	CALL  reeighenvalue(m,c,ev,iter)
	minabsev=1/ev          !A的按模最小特征值
	condA2=(maxabsev/minabsev)
	write(*,*)  'cond(A)2=' ,condA2
	write(1,*)  'cond(A)2=' ,condA2

	call LUDOOLITTLE(m,c,lu)
	s=1
	do i=1,m
	s=s*lu(3,i)
	end do
	write(*,*) 'detA=',s
	write(1,*) 'detA=',s

	stop
      end
!幂法
	SUBROUTINE eighenvalue(m,c,ev)
	real,parameter::er=1.e-12
      real c(m+2),u(m),y(m),b(1000000)
	real s,a
      integer i,j,k
	real error
	error=1


	do i=1,m          !initial vector
	u(i)=1
	end do
	u(m)=9

	k=2
	b(1)=10000

	do while (error.GT.er)


	s=0
	do i=1,m
	s=s+u(i)**2
	end do
	a=sqrt(s)

	do i=1,m
	y(i)=u(i)/a
	end do

	do i=1,m
	s=0.
	if(i.EQ.1) then
	s=c(i)*y(i)+c(m+1)*y(i+1)+c(m+2)*y(i+2)
	else if(i.EQ.2) then
	s=c(i)*y(i)+c(m+1)*(y(i+1)+y(i-1))+c(m+2)*y(i+2)
	else if(i.GT.2.AND.i.LT.(m-1)) then
	s=c(i)*y(i)+c(m+1)*(y(i+1)+y(i-1))+c(m+2)*(y(i+2)+y(i-2))
	else if(i.EQ.(m-1)) then
	s=c(i)*y(i)+c(m+1)*(y(i+1)+y(i-1))+c(m+2)*(0+y(i-2))
	else if(i.EQ.(m)) then
	s=c(i)*y(i)+c(m+1)*(0+y(i-1))+c(m+2)*(0+y(i-2))
	end if
	u(i)=s
	end do

	s=0.
	do i=1,m
	s=s+y(i)*u(i)
	end do 
	b(k)=s

	error=abs(b(k)-b(k-1))
	k=k+1
	end do
		
c	write(*,*)k
	k=k-1
	ev=b(k)
c	write(*,*) ev
	

	END


      
	SUBROUTINE LUDOOLITTLE(m,c,lu)
	real c(m+2),lu(5,m)
      integer i,j
	do i=1,5
	do j=1,m
	lu(i,j)=0
	end do
	end do

!lu decomposition
      i=1
c    compute u(i,j)
	lu(3,i)=c(i)
	lu(4,i+1)=c(m+1)
	lu(5,i+2)=c(m+2)
c    compute l(i,j)
      lu(2,i+1)=(c(m+1)-0)/lu(3,i)
	lu(1,i+2)=c(m+2)/lu(3,i)

	i=2
c    compute u(i,j)
	lu(3,i)=c(i)-lu(2,i)*lu(4,i)  
	lu(4,i+1)=c(m+1)-lu(2,i)*lu(5,i+1)
      lu(5,i+2)=c(m+2)
c    compute l(i,j)
	 lu(2,i+1)=(c(m+1)-lu(1,i-1)*lu(4,i))/lu(3,i)
	 lu(1,i+2)=c(m+2)/lu(3,i)

	do i=3,m-2
c    compute u(i,j)
	lu(3,i)=c(i)-lu(1,i)*lu(5,i)-lu(2,i)*lu(4,i)  
	lu(4,i+1)=c(m+1)-lu(2,i)*lu(5,i+1)
	lu(5,i+2)=c(m+2)
c    compute l(i,j)
      lu(2,i+1)=(c(m+1)-lu(1,i+1)*lu(4,i))/lu(3,i)
	lu(1,i+2)=c(m+2)/lu(3,i)
	end do

	i=m-1
c    compute u(i,j)
c	lu(i,i)=c(i)-lu(i,i-2)*lu(i-2,i)-lu(i,i-1)*lu(i-1,i)
	lu(3,i)=c(i)-lu(1,i)*lu(5,i)-lu(2,i)*lu(4,i)
	lu(4,i+1)=c(m+1)-lu(2,i)*lu(5,i+1)
c    compute l(i,j)
      lu(2,i+1)=(c(m+1)-lu(1,i+1)*lu(4,i))/lu(3,i)

      i=m
c    compute u(i,j)
	lu(3,i)=c(i)-lu(1,i)*lu(5,i)-lu(2,i)*lu(4,i)


c	do i=1,5
c	do j=1,m
c	write(*,*)i,j,lu(i,j)
c	end do
c	end do

	END

      subroutine INVERSELU(m,lu,b,x)
	real lu(5,m),x(m),y(m),b(m)
	real s,a
      integer i,j,k,t

	y(1)=b(1)
	i=2
	y(i)=b(i)-y(i-1)*lu(2,i)
	do i=3,m
	y(i)=b(i)-y(i-2)*lu(1,i)-y(i-1)*lu(2,i)
	end do

c	do i=1,m
c	write(*,*) 'i=',i,y(i)
c	end do
c right
	x(m)=y(m)/lu(3,m)
	i=m-1
      x(i)=(y(i)-x(i+1)*lu(4,i+1))/lu(3,i)
	do i=m-2,1,-1
 	x(i)=(y(i)-x(i+2)*lu(5,i+2)-x(i+1)*lu(4,i+1))/lu(3,i)
	end do


	END

	subroutine reeighenvalue(m,c,ev,iter)
	real,parameter::er=1.e-12
      real c(m+2),lu(5,m),u(m),y(m),b(1000000)
	real s,a
      integer i,j,k,iter
	real error
	error=1


	do i=1,m          !initial vector
	u(i)=1
	end do
	u(m)=1

	k=2
	b(1)=10000
		
	do while (error.GT.er)
      call LUDOOLITTLE(m,c,lu)  !将矩阵c分解为LU(M,M)	
	s=0
	do i=1,m
	s=s+u(i)**2
	end do
	a=sqrt(s)

	do i=1,m
	y(i)=u(i)/a
	end do
	
	call INVERSELU(m,lu,y,u)

	s=0.
	do i=1,m
	s=s+y(i)*u(i)
	end do 
	b(k)=s

	error=abs(b(k)-b(k-1))
	k=k+1
	end do
		
c	write(*,*)k
	k=k-1
	iter=k
	ev=b(k)
c	write(*,*) ev
	
	END


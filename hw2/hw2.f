 !!!!!!!!!!!!!!!!!!!!!!!!!
!hw 2双位移法求特征值
!!!!!!!!!!!!!!!!!!!!!!!!!
      program step
	implicit none
      integer,parameter::n=4
	integer,parameter::it=10000  !控制最大迭代步数
	real,parameter::er=1.E-12  !控制最大迭代步数
	real B(n,n),A(n,n),Q(n,n)
	REAL EE(n,n)                !用于检查 
!A迭代过程中主对角线元素,lamb(i,0)表示复数，lamb(i,1)表示实数
	real lamb(n,2)             
	integer i,j,k,m,iter,t
	real s(3),pa,pb,pc          !一元二次方程三个系数及两个根
	WRITE(*,*)'ER',er
	
	OPEN(1,FILE='CHECKhw.2.txt',STATUS='UNKNOWN')

	do i=1,n
	do j=1,n

	B(i,j)=0
	if(i.eq.j)  then
	B(i,j)=1
	end if
	if (i.gt.(j)) then 
	B(i,j)=1000
	end if

	end do 
	end do





	do i=1,n
	do j=1,n
	if (i.NE.j) then 
!	B(i,j)=sin(0.5*i+0.2*j)
	else 
!	B(i,j)=1.52*cos(i+1.2*j)
	end if
	end do
	end do



	do i=1,n
	do j=1,n
	A(i,j)=B(i,j)
	end do
	end do
!	CALL nishangsanjiaohua(A,n)

      write(1,*)'initial A'
	DO i=1,n
	WRITE(1,*)(A(i,j),j=1,n)
	end do

	k=1
	m=n
	t=3

      do iter=1,it
	write(*,*) iter

	if (m.eq.3) then
	write(*,*)'kankan'
	end if
!第三步
	if(t.EQ.3) THEN    
	   if (abs(A(m,m-1)).LT.er) then 
          	lamb(m,1)=A(m,m)
	        lamb(m,2)=0
		   write(*,*)'three '
	       write(*,*)'m=',m
	       write(*,*)lamb(m,1),lamb(m,2)

	        m=m-1
	        t=4
	   else 
	        t=5
	   end if
	end if
!第四步
	if(t.EQ.4) THEN
	   if (m.LT.1) then 
               exit
	   end if

	   if (m.EQ.1) then 
	         lamb(m,1)=A(m,m)
	         lamb(m,2)=0
		     exit
         end if

	   if(m.GT.1) then 
		     t=3
         end if
      end if
!第五步
      if(t.EQ.5) THEN
	   pa=1
	   pb=-(A(m-1,m-1)+A(m,m)) 
	   pc=A(m-1,m-1)*A(m,m)-A(m,m-1)*A(m-1,m)
	   CALL soleq(pa,pb,pc,s)
	   t=6
	END IF
!第六步
	if(t.EQ.6) THEN
	   if (m.EQ.2) then
	      if(s(3).EQ.1) then
	        lamb(m,1)=s(1)
	        lamb(m,2)=0
	        lamb(m-1,1)=s(2)
              lamb(m-1,2)=0
	      else 
              lamb(m,1)=s(1)
              lamb(m,2)=s(2)
              lamb(m-1,1)=s(1)
              lamb(m-1,2)=-s(2)
	      end if

		   write(*,*)'six '
	       write(*,*)'m=',m
	       write(*,*)lamb(m,1),lamb(m,2)
	       write(*,*)lamb(m-1,1),lamb(m-1,2)

	      exit
	    else 
		  t=7
	    end if
      end if
!第七步
      if(t.EQ.7) THEN
	    if(abs(A(m-1,m-2)).LT.er) then
       	   if(s(3).EQ.1) then
	         lamb(m,1)=s(1)
	         lamb(m,2)=0
	         lamb(m-1,1)=s(2)
               lamb(m-1,2)=0
	       else 
               lamb(m,1)=s(1)
               lamb(m,2)=s(2)
               lamb(m-1,1)=s(1)
               lamb(m-1,2)=-s(2)
	       end if

	       write(*,*)'seven '
	       write(*,*)'m=',m
	       write(*,*)lamb(m,1),lamb(m,2)
	       write(*,*)lamb(m-1,1),lamb(m-1,2)

	       m=m-2
	       t=4
	    else 
             t=9
	    end if
	end if
!第九步
      if(t.EQ.9) THEN
     	    CALL DQR(A,m,n)
       	k=k+1
          t=3
	end if

	end do

5     write(*,*)'end '
      
	write(*,*)'end ----------'
	do i=1,n
	write(*,*)(lamb(i,j),j=1,2)
	end do

	if (k.eq.200) then
!	write(1,*)'k=',k
	DO i=1,n
!	WRITE(1,*)(A(i,j),j=1,n)
	end do
	end if






 
!	do k=1,1
!	call DQR(A,m,n)
!	end do
	
!	do i=1,n
!	write(*,*) (A(i,j),j=1,n)
!	end do






!	WRITE(*,*) S(1),S(2),S(3)


	stop 
	end 

      subroutine MULT(A,Q,n)
	real A(n,n),Q(n,n),R(n,n),EE(n,n)
	integer i,j,k
	real s
	
	do i=1,n
	do j=1,n
	s=0
	do k=1,n
	s=s+A(i,k)*Q(k,j)         
	end do 
	R(i,j)=s
	end do
	end do

	do i=1,n
	do j=1,n
	A(i,j)=R(i,j)
	end do 
	end do

c	write(1,*) 'neibukankankan'
c	do j=1,n
c	WRITE(1,*) (A(j,k),k=1,n)
c	end do 

	end 


	subroutine QR(A,Q,n)

	real A(n,n),Q(n,n),B(n,n)
      real EE(n,n)   !EE=Q*R to check if QR DECOMPOSITION IF RIght
	real u(n),p(n),w(n),h,c       
	real s                 !for plus and 
	integer i,j,r,k

	
	
	do i=1,n
	do j=1,n
	IF(i.NE.j) then
	               Q(i,j)=0
	           else 
	               Q(i,j)=1
	END IF
	end do 
	end do



	do r=1,n-1      !n-1
	s=0
	do i=r+1,n
	s=s+A(i,r)**2
	end do
!	WRITE(*,*)'S=',S

	if(s.NE.0) THEN 
	s=0
	do i=r,n             !9:42
	s=s+A(i,r)**2
	end do
	d=sqrt(s)
c	WRITE(*,*)'D=',d
	if(A(r,r).EQ.0) THEN 
	c=d
	else if(A(r,r).GT.0) THEN 
	c=-d
	else if(A(r,r).LT.0) THEN 
	c=d
	end if

	h=c**2-c*A(r,r)
!对u(r)赋值
      do i=1,r-1
	u(i)=0
	end do 
	u(r)=A(r,r)-c
	do i=r+1,n
	u(i)=A(i,r)
	end do

!w=Q*u
      do i=1,n
	s=0
	do j=1,n
	s=s+Q(i,j)*u(j)
	end do
	w(i)=s
	end do

!Q=Q-W*u/h
      do i=1,n
	do j=1,n
	Q(i,j)=Q(i,j)-w(i)*u(j)/h
	end do
	end do

! p=At*u/h
      do i=1,n
	s=0
	do j=1,n
	s=s+A(j,i)*u(j)
	end do
	p(i)=s/h
	end do

!A=A-u*p
      do i=1,n
	do j=1,n
	A(i,j)=A(i,j)-u(i)*p(j)
	end do
	end do

	end if
	end do

! CHECK AND OUTPUT
!      WRITE(*,*)'A(I,J)'
!	do i=1,n
!	write(*,*)(B(i,j),j=1,n)
!	end do

c      WRITE(*,*)'Q(I,J)'
c	do i=1,n
c	write(*,*)(Q(i,j),j=1,n)
c	end do

c     WRITE(*,*)'R(I,J)'
c	do i=1,n
c	write(*,*)(A(i,j),j=1,n)
c	end do
	
c     WRITE(*,*)'A=Q*R'          !NOT RIGHT
c	do i=1,n
c	do j=1,n
c	s=0
	
c	do k=1,n
c	s=s+Q(i,k)*A(k,j)
c	end do 

c	EE(i,j)=s
c	end do
c	end do

c	do i=1,n
c	write(*,*)(EE(i,j),j=1,n)
c	end do

	end

	subroutine soleq(a,b,c,s)
	real a,b,c,s(3),d
	d=b**2-4*a*c
!	WRITE(*,*) a,b,c,D

	if(d.GE.0) then
	s(1)=(-b+sqrt(d))/(2*a)
	s(2)=(-b-sqrt(d))/(2*a)
	s(3)=1
	else 
	s(1)=-b/(2*a)
	s(2)=sqrt(-d)/(2*a)
	s(3)=0
	end if
	END 

	subroutine DQR(A,m,n)
	real s,t,A(n,n),MM(n,n),Q(n,n),R(n,n)
	real B(n,n),E(n,n),C(n,n)
	real p       !用于累加
	integer i,j,k,m,n

	s=A(m-1,m-1)+A(m,m)
	t=A(m-1,m-1)*A(m,m)-A(m-1,m)*A(m,m-1)
!	write(*,*)'s=',s,'t=',t


	do i=1,n
      do j=1,n
	if(i.NE.j) then
	E(i,j)=0
	else 
	E(i,j)=1
	end if 
	end do
	end do

	do i=1,n
	do j=1,n
	p=0
	do k=1,n
	p=p+A(i,k)*A(k,j)
	end do
	MM(i,j)=p-s*A(i,j)+t*E(i,j)
	end do
	end do

!	do i=1,n
!	write(*,*) (MM(i,j),j=1,n)
!	end do
!22:23全部正确
	call QR(MM,Q,n)

!	do i=1,n
!	write(*,*) (MM(i,j),j=1,n)
!	end do

!B=Qt*A
      do i=1,n
	do j=1,n
	p=0
	do k=1,n
	p=p+Q(k,i)*A(k,j)
	end do
	B(i,j)=p
	end do 
	end do
!C=B*Q
      do i=1,n
	do j=1,n
	p=0
	do k=1,n
	p=p+B(i,k)*Q(k,j)
	end do
	C(i,j)=p
	end do
	end do

	do i=1,n
	do j=1,n
	A(i,j)=C(i,j)
	end do
	end do

	END 

	
	subroutine nishangsanjiaohua(A,n)
	real A(n,n),p(n),q(n),w(n),u(n)
	real d,s,c,h,t
	integer i,j,r

	do r=1,n-2

	s=0
	do i=r+2,n
	s=s+A(i,r)**2
	end do

	if(s.NE.0) THEN 
	d=0
	do i=r+1,n
	d=d+A(i,r)**2
	end do
	d=sqrt(d)

	if(A(r+1,r).EQ.0) THEN 
	c=d
	else if(A(r+1,r).GT.0) THEN 
	c=-d
	else if(A(r+1,r).LT.0) THEN 
	c=d
	end if

	h=c**2-c*A(r+1,r)
!对u(r)赋值
      do i=1,r
	u(i)=0
	end do

	do i=r+2,n
	u(i)=A(i,r)
	end do

	u(r+1)=A(r+1,r)-c
!p=At*u/h
      do i=1,n
	s=0
	do j=1,n
	s=s+A(j,i)*u(j)
	end do

	p(i)=s/h
	end do
!q=A*u/h
      do i=1,n
	s=0
	do j=1,n
	s=s+A(i,j)*u(j)
	end do
	q(i)=s/h
	end do
!t=Pt*u/h
      s=0
	do i=1,n
	s=s+p(i)*u(i)
	end do
	t=s/h
!w=q-t*u
      do i=1,n
	w(i)=q(i)-t*u(i)
	end do
!A=A-w*Ut-u*Pt
      do i=1,n
	do j=1,n
	A(i,j)=A(i,j)-w(i)*u(j)-u(i)*p(j)
	end do
	end do

!	write(1,*)'r=',r,'A='
!	DO i=1,n
!	write(1,*) (A(i,j),j=1,n)
!	end do

	END IF
	END DO

	END






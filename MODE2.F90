PROGRAM MAIN
!Matheus Araujo Souza ra:184145 Modelagem Matemática: ms480
implicit none
INTEGER , PARAMETER :: n =4
DOUBLE PRECISION, PARAMETER :: MG = 0.45
DOUBLE PRECISION, PARAMETER :: t1=0.75
DOUBLE PRECISION, PARAMETER :: t2=1.00
DOUBLE PRECISION, PARAMETER :: t3=1.25
DOUBLE PRECISION, PARAMETER :: t4=1.50
INTEGER, PARAMETER :: ti1=15
INTEGER,PARAMETER :: ti2=20
INTEGER, PARAMETER :: ti3=25
INTEGER,PARAMETER :: ti4=30
DOUBLE PRECISION :: x(n)
DOUBLE PRECISION, DIMENSION(35,n) :: guard !Matriz de armazenamento de dados         
DOUBLE PRECISION, DIMENSION(35,n) :: ERRO
DOUBLE PRECISION, DIMENSION(35,n+1) :: Resp2
DOUBLE PRECISION , PARAMETER :: dx = 0.1
INTEGER :: icont,g,i,j,F,maxiter  
DOUBLE PRECISION , PARAMETER :: dt = 0.05 
DOUBLE PRECISION, DIMENSION(n,n) :: yacobian       
DOUBLE PRECISION :: fminus(n), fplus(n), xnew(n),h(n),xtol




guard(:,:)=0.0

!x(1)=0;x(2)=0;x(3)=0;x(4)=0;x(5)=0;x(6)=0;x(7)=0;x(8)=0;x(9)=0;x(10)=0
x(:)=0.0             !atribuindo valores para todo vetor como zero 
!ainda preciso mandar a Matriz principal para dentro do codigo e depois rodar o loop indo ate 2 para resolver na segunda linha do tempo
xtol=0.1E-15; maxiter=100
DO i=1,35
icont=i
CALL DAMPED_NEWTON(x,n,xtol,maxiter,icont,dx,dt,guard,MG,t1,t2,t3,t4,ti1,ti2,ti3,ti4,yacobian,fminus,fplus,xnew,h)
!armazenando dados da saida na Matriz resultados finais

	DO j=1,n
		guard(i,j)=x(j)
	END DO




	if (i>=0 .AND. i< ti1) THEN
		Resp2(i,1)=0
	END IF
	if (i>=ti1 .AND. i< ti2) THEN
		Resp2(i,1)=(MG*((i*dt-t1)/(t2-t1)))
		
	END IF
	if (i>=ti2 .AND. i<= ti3) THEN
		Resp2(i,1)=MG
	END IF
	if (i>ti3 .AND. i<= ti4) THEN
		Resp2(i,1)=(MG*((t4-(i*dt))/(t4-t3)))
	END IF
	if (i>ti4) THEN
		Resp2(i,1)=0
	END IF

END DO

DO F=1,35
	DO J=1,n
		Resp2(F,J+1)=guard(F,J)
	END DO
END DO



!PRINT*, (x(i), i=1,n) 
!PRINT*, (x(i), i=1,n) 

DO g=1,35
PRINT*, Resp2(g,:)
END DO

!DO i=1,100
!	DO j=1,n
!		ERRO(i,j)=ABS(sin((j*vdx)*(i*vdt))  - Resp(i,j))
!	END DO
!END DO

!DO g=1,100
!PRINT*, ERRO(g,:)
!END DO

!PRINT *,'Maior erro encontrado', MAXVAL(ERRO)
PRINT *,'dx=',dx,'dt=',dt,'tmax=1 e dmax=1'
PRINT*,'valor de M',MG,'intervalos da rampa T1, T2, T3, T4',t1,t2,t3,t4
END
!****************************************
SUBROUTINE fn(x,f,n,i,dx,dt,guard,MG,gt1,gt2,gt3,gt4,ti1,ti2,ti3,ti4)
implicit none
DOUBLE PRECISION :: x(n), f(n), guard(35,n)
DOUBLE PRECISION :: dx, dt
DOUBLE PRECISION :: MG, gt1, gt2, gt3, gt4
INTEGER :: ti1, ti2, ti3, ti4, i, j, n
!PRINT *, MG, ti1,ti2,ti3,ti4
!PRINT *, MG, gt1, gt2, gt3, gt4
!PAUSE
!gerando resultados para o primeiro nivel do tempo
!f(i) = (u(i,t+1) - u(i,t))/dt + ((u(i+1,t+1) - u(i+1,t-1) )/(2.d0*dx))**2 
!- x(i)*dcos(x(i)*t)  - t**2*dcos(x(i)*t)**2
if(i==1) THEN
	 F(1)= (x(1))*4*((dx)**2)   +  dt*(x(2))**2   - dt*4*((dx)**2) * (1*dx)*cos(dx*dt) &
	 		 - dt*4*((dx)**2)*((dt**2)*(cos(dx*dt))**2)
        DO j=2,n - 1
        F(j)= (x(j))*4*((dx)**2)   +  dt*(x(j+1))**2 - dt*2*(x(j+1)*(x(j-1))) &
        + dt*((x(j-1))**2) - dt*4*((dx)**2) *(j*dx)*cos((j*dx)*dt)- dt*4*((dx)**2)*((dt**2)*(cos((j*dx)*dt))**2) 
        END DO
        F(n)= (x(n))*((dx)**2)   +  dt*(x(n))**2 - dt*2*(x(n)*(x(n-1))) &
        + dt*(x(n-1))**2   - dt*(dx)**2 * (n*dx)*cos((n*dx)*dt) - dt*((dx)**2)*((dt**2)*(cos((n*dx)*dt))**2) 
END IF



if(i > 1 .AND. i < (ti1)) THEN
	

	F(1)= (x(1))*4*((dx)**2)   +  dt*((x(2))**2)   - dt*4*((dx)**2) * (((1*dx)*cos((1*dx)*((i)*dt))) &
		+ ( (((i)*dt)**2)*((cos((1*dx)*((i)*dt)))**2)) ) - (4*(dx**2) * guard(i-1,1))
     
    DO j=2,n - 1 

     F(j)= x(j)*4*(dx**2)   +  dt*((x(j+1))**2) - dt*2*(x(j+1)*(x(j-1))) &
     + dt*((x(j-1))**2)  - dt*4*((dx)**2) * (((j*dx)*cos((j*dx)*((i)*dt))) &
     	+ ( (((i)*dt)**2)*((cos((j*dx)*((i)*dt)))**2)) )  -4*(dx**2)*(guard(i-1,j))

               
    END DO 
    F(n)=(x(n))*((dx)**2)   +  dt*((x(n))**2) - dt*2*(x(n)*x(n-1)) + dt*((x(n-1))**2)   &
    - dt*(dx)**2 * ((((n*dx)*cos((n*dx)*((i)*dt))) +  (((i)*dt)**2)*(cos(n*dx*(i)*dt))**2 )) - (dx**2)*(guard(i-1,n))
    
END IF
!como tenho que meu i é um contador de tempo vou criar varias restrições e colocar minhas novas distretizações            

if(i >= (ti1) .AND. i < (ti2)) THEN


	F(1)= (x(1))*4*(dx**2) + dt*((x(2))**2) - dt*2*(x(2))*(MG*((i*dt-gt1)/(gt2-gt1)))  &
	    +dt*((MG*((i*dt-gt1)/(gt2-gt1)))**2) &
		- dt*4*(dx**2) * (((1*dx)*cos((1*dx)*((i)*dt))) &
		+ ( (((i)*dt)**2)*((cos((1*dx)*((i)*dt)))**2)) ) - (4*(dx**2) * guard(i-1,1))
     
    DO j=2,n - 1 

     F(j)= x(j)*4*(dx**2) + dt*((x(j+1))**2) - dt*2*(x(j+1)*(x(j-1))) &
     + dt*((x(j-1))**2)  - dt*4*((dx)**2) * (((j*dx)*cos((j*dx)*((i)*dt))) &
     	+ ( (((i)*dt)**2)*((cos((j*dx)*((i)*dt)))**2)) )  -4*(dx**2)*(guard(i-1,j))

               
    END DO 
    F(n)=(x(n))*((dx)**2)   +  dt*((x(n))**2) - dt*2*(x(n)*x(n-1)) + dt*((x(n-1))**2)   &
    - dt*(dx)**2 * ((((n*dx)*cos((n*dx)*((i)*dt))) +  (((i)*dt)**2)*(cos(n*dx*(i)*dt))**2 )) - (dx**2)*(guard(i-1,n))
    
END IF




if(i >= (ti2) .AND. i <=(ti3)) THEN

	F(1)= (x(1))*4*(dx**2)   +  dt*(x(2)**2) - dt*2*x(2)*(MG) + dt*(MG**2)  &
		- dt*4*((dx)**2) * (((1*dx)*cos((1*dx)*((i)*dt))) &
		+ ( (((i)*dt)**2)*((cos((1*dx)*((i)*dt)))**2)) ) - (4*(dx**2) * guard(i-1,1))
     
    DO j=2,n - 1 

     F(j)= x(j)*4*(dx**2)   +  dt*((x(j+1))**2) - dt*2*(x(j+1)*(x(j-1))) &
     + dt*((x(j-1))**2)  - dt*4*((dx)**2) * (((j*dx)*cos((j*dx)*((i)*dt))) &
     	+ ( (((i)*dt)**2)*((cos((j*dx)*((i)*dt)))**2)) )  -4*(dx**2)*(guard(i-1,j))

               
    END DO 
    F(n)=(x(n))*((dx)**2)   +  dt*((x(n))**2) - dt*2*(x(n)*x(n-1)) + dt*((x(n-1))**2)   &
    - dt*(dx)**2 * ((((n*dx)*cos((n*dx)*((i)*dt))) +  (((i)*dt)**2)*(cos(n*dx*(i)*dt))**2 )) - (dx**2)*(guard(i-1,n))
    
END IF





if(i > (ti3) .AND. i <= (ti4)) THEN

	F(1)= (x(1))*4*((dx)**2) + dt*((x(2))**2) - dt*2*x(2)*(MG*((gt4 - (i*dt))/(gt4-gt3)))  &
	      +dt*((MG*((gt4 - (i*dt))/(gt4-gt3)))**2) &
		- dt*4*((dx)**2) * (((1*dx)*cos((1*dx)*((i)*dt))) &
		+ ( (((i)*dt)**2)*((cos((1*dx)*((i)*dt)))**2)) ) - (4*(dx**2) * guard(i-1,1))
     
    DO j=2,n - 1 

     F(j)= x(j)*4*(dx**2)   +  dt*((x(j+1))**2) - dt*2*(x(j+1)*(x(j-1))) &
     + dt*((x(j-1))**2)  - dt*4*((dx)**2) * (((j*dx)*cos((j*dx)*((i)*dt))) &
     	+ ( (((i)*dt)**2)*((cos((j*dx)*((i)*dt)))**2)) )  -4*(dx**2)*(guard(i-1,j))

               
    END DO 
    F(n)=(x(n))*((dx)**2)   +  dt*((x(n))**2) - dt*2*(x(n)*x(n-1)) + dt*((x(n-1))**2)   &
    - dt*(dx)**2 * ((((n*dx)*cos((n*dx)*((i)*dt))) +  (((i)*dt)**2)*(cos(n*dx*(i)*dt))**2 )) - (dx**2)*(guard(i-1,n))
    
END IF

if(i > (ti4) ) THEN

	F(1)= (x(1))*4*((dx)**2)   +  dt*((x(2))**2)   - dt*4*((dx)**2) * (((1*dx)*cos((1*dx)*((i)*dt))) &
		+ ( (((i)*dt)**2)*((cos((1*dx)*((i)*dt)))**2)) ) - (4*(dx**2) * guard(i-1,1))
     
    DO j=2,n - 1 

     F(j)= x(j)*4*(dx**2)   +  dt*((x(j+1))**2) - dt*2*(x(j+1)*(x(j-1))) &
     + dt*((x(j-1))**2)  - dt*4*((dx)**2) * (((j*dx)*cos((j*dx)*((i)*dt))) &
     	+ ( (((i)*dt)**2)*((cos((j*dx)*((i)*dt)))**2)) )  -4*(dx**2)*(guard(i-1,j))

               
    END DO 
    F(n)=(x(n))*((dx)**2)   +  dt*((x(n))**2) - dt*2*(x(n)*x(n-1)) + dt*((x(n-1))**2)   &
    - dt*(dx)**2 * ((((n*dx)*cos((n*dx)*((i)*dt))) +  (((i)*dt)**2)*(cos(n*dx*(i)*dt))**2 )) - (dx**2)*(guard(i-1,n))
    
END IF




RETURN
END SUBROUTINE fn
!***************************************

SUBROUTINE JACOBIAN(x,yacobian,n,ic,dx,dt,guard,MG,gt1,gt2,gt3,gt4,ti1,ti2,ti3,ti4,fminus,fplus)
implicit none
DOUBLE PRECISION :: x(n), yacobian(n,n), fplus(n), fminus(n) ,guard(35,n)
DOUBLE PRECISION :: dx, dt 
DOUBLE PRECISION :: MG, gt1, gt2, gt3, gt4,xj,hj  
INTEGER :: ti1, ti2, ti3, ti4, ic, i, j, n
!PRINT*,x,fplus,fminus                          
!PRINT *, dx, dt, ic
fplus=0
fminus=0
hj=0.0001
DO i=1,n
DO j=1,n
xj=x(j); x(j)=xj+hj; CALL fn(x,fplus,n,ic,dx,dt,guard,MG, gt1, gt2, gt3, gt4,ti1,ti2,ti3,ti4)
x(j)=xj-hj; CALL fn(x,fminus,n,ic,dx,dt,guard,MG, gt1, gt2, gt3, gt4,ti1,ti2,ti3,ti4); x(j)=xj
yacobian(i,j)=(fplus(i)-fminus(i))/(2*hj)
END DO
END DO
END SUBROUTINE JACOBIAN
!***************************************


SUBROUTINE DAMPED_NEWTON(x,n,xtol,maxiter,ic,dx,dt,guard,MG,gt1,gt2,gt3,gt4,ti1,ti2,ti3,ti4,yacobian,fminus,fplus,xnew,h)
implicit none
DOUBLE PRECISION :: x(n), f(n), h(n), yacobian(n,n), xnew(n),guard(35,n),fminus(n),fplus(n)
INTEGER :: ic, i, j, k, n, ti1, ti2, ti3, ti4, maxiter
DOUBLE PRECISION :: dx, dt
DOUBLE PRECISION :: MG, gt1, gt2, gt3, gt4,absdiff,xtol

DO k=1,maxiter
	

	CALL fn(x,f,n,ic,dx,dt,guard,MG, gt1, gt2, gt3, gt4,ti1,ti2,ti3,ti4)
	
	CALL JACOBIAN(x,yacobian,n,ic,dx,dt,guard,MG,gt1,gt2,gt3,gt4,ti1,ti2,ti3,ti4,fminus,fplus)

	CALL GAUSS(n,yacobian,f)
	
        
	DO i=1,n
	h(i)= - f(i)
	END DO

	xnew= x + h
    
	DO i=1,n
		absdiff=ABS(xnew(i)-x(i))                               
	END DO    
	IF(absdiff<xtol) EXIT

	
	DO i=1,n
		x(i)=xnew(i)
	END DO


END DO
 
END SUBROUTINE DAMPED_NEWTON

!************************************
SUBROUTINE GAUSS(n,a,b)
! n=number of unknowns and equations. (Input)
! a=n×n Matrix of the coefficients. (Input)
! b=right hand side vector. (Input)
! The solution is returned in vector b. (Output)
! s(i)=size of the row i of Matrix a.
!=======================================================================================
! Solucao de um sistema linear A*xsol = b usando
! eliminacao gaussiana com pivotamento parcial
!---------------------------------------------------------------------------------------
! Entrada <----
! a(n,n)   :: matriz A
! b(n)     :: vetor b do sistema A*xsol = b
! n        :: numero de indeterminadas (tamanho da matriz A)
! saida   ---->
! xsol(n)  :: vetor de solucoes 
!=======================================================================================
implicit none 
INTEGER         ::  n, i, j, k, l
DOUBLE PRECISION::  a(n,n), b(n), s(n), c, pivot, m


!do i=1,4
!	PRINT*,a(i,:)
!end do

! paso 1: eliminacao
do k=1, n-1

! paso 2: escalonamento

  do i=k,n                       
    s(i) = 0.0
    do j=k,n                    
      s(i) = max(s(i),abs(a(i,j)))
    end do
  end do

! paso 3: pivotando 

  pivot = abs(a(k,k)/s(k))
  l = k
  do j=k+1,n
    if(abs(a(j,k)/s(j)) > pivot) then
      pivot = abs(a(j,k)/s(j))
      l = j
    end if
  end do

! verificando se a matriz do sistema é singular
  if(pivot == 0.0) then
    PRINT*, ' A matriz e singular '
    return
  end if

! step 4: pivotando e permutando linhas (se necessario )
if (l /= k) then
  do j=k,n
          m = a(k,j)
     a(k,j) = a(l,j)
     a(l,j) = m
  end do
     m = b(k)
  b(k) = b(l)
  b(l) = m
end if

! paso 5: eliminacao (depois do escalonamento e pivotamento)
   do i=k+1,n
      c=a(i,k)/a(k,k)
      a(i,k) = 0.0
      b(i)=b(i)- c*b(k)
      do j=k+1,n
         a(i,j) = a(i,j)-c*a(k,j)
      end do
   end do
end do

! paso 6: back substitution 
b(n) = b(n)/a(n,n)
do i=n-1,1,-1
   c=0.0
   do j=i+1,n
     c= c + a(i,j)*b(j)
   end do 
   b(i) = (b(i)- c)/a(i,i)
end do



return

END SUBROUTINE GAUSS
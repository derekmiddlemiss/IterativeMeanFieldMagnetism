PROGRAM iterative_mean_field

IMPLICIT none

INTEGER :: i,j,k,l,z,num_sites,max_shell,max_neigh,count,Tcount,Tstop
LOGICAL :: debug,converged,low_x_brill,momkeep
DOUBLE PRECISION :: bril,bril_low,B,T,converge,damping,max_diff,max_diff_last,g
DOUBLE PRECISION :: new_mom,diff,arg,Tstart,Tstep
DOUBLE PRECISION, ALLOCATABLE :: mom(:),mom_last(:),exch_field(:)
DOUBLE PRECISION, ALLOCATABLE :: spin(:),start_mom(:),J_list(:,:) 
INTEGER, ALLOCATABLE :: neigh_list(:,:,:),num_shells(:),num_neigh(:,:)

! Read initial data from INPUT file
OPEN (unit=10,file='control.in',status='old')
READ(10,*) debug
READ(10,*) low_x_brill
READ(10,*) momkeep
READ(10,*) converge,damping
READ(10,*) B,Tstart,Tstop,Tstep,g
READ(10,*) num_sites

ALLOCATE(mom(num_sites),mom_last(num_sites))
ALLOCATE(start_mom(num_sites),spin(num_sites))

READ(10,*) max_shell,max_neigh

ALLOCATE(neigh_list(num_sites,max_shell,max_neigh),num_shells(num_sites),num_neigh(num_sites,max_shell))
ALLOCATE(J_list(num_sites,max_shell),exch_field(num_sites))

neigh_list = 0
num_shells = 0
num_neigh = 0
J_list = 0d0

DO i=1,num_sites,1
  READ(10,*) spin(i), num_shells(i)

  IF (num_shells(i) .gt. max_shell) THEN
    WRITE(*,*) "Problem. Site",i,"num_shells exceeds max_shell. Program stopping..."
    STOP
  END IF

  IF (debug .eqv. .true.) THEN
    WRITE(*,*) "--------------------------------------------------------------------------------------------"
    WRITE(*,*) "****Site",i,"****"
    WRITE(*,*) "Spin,num_shells--",spin(i),num_shells(i)
  END IF

  DO j=1,num_shells(i),1
    READ(10,*) num_neigh(i,j), (neigh_list(i,j,k),k=1,num_neigh(i,j),1), J_list(i,j)

    IF (num_neigh(i,j) .gt. max_neigh) THEN
      WRITE(*,*) "Problem. Site",i,"shell",j,"num_neigh exceeds max_neigh. Program stopping..."
      STOP
    END IF

    IF (debug .eqv. .true.) THEN
      WRITE(*,*)
      WRITE(*,*) "Shell",j,"num_neigh--",num_neigh(i,j)
      WRITE(*,*) "neigh_list--",(neigh_list(i,j,k),k=1,num_neigh(i,j),1)
      WRITE(*,*) "J_list--",J_list(i,j)
    END IF

  END DO
END DO

IF(debug .eqv. .true.) THEN
  WRITE(*,*) "--------------------------------------------------------------------------------------------"
END IF

DO i=1,num_sites,1
  READ(10,*) start_mom(i)
END DO

OPEN(unit=11,file="output",status='replace')

! Begin temperature loop
DO Tcount=1,Tstop,1
  
  T=Tstart+(DBLE(Tcount-1)*Tstep)
  
  converged= .false.

  IF (Tcount .gt. 1 .AND. momkeep .eqv. .true.) THEN
  ELSE
    mom=start_mom
  END IF

  WRITE(*,*)
  WRITE(*,*) "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
  WRITE(*,*) "TTTT Temperature",T,"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
  WRITE(*,*) "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
  WRITE(*,*) "Starting moments"
  DO i=1,num_sites,1
    WRITE(*,*) "Site",i,"--",mom(i)
  END DO
  
  count=0

! Begin main loop
  DO WHILE (converged .eqv. .false.)

    count=count+1
    max_diff=0d0

    WRITE(*,*)
    WRITE(*,*) "********************************************************************************************"
    WRITE(*,*) "****Iteration number", count,"***********************************************************"
    WRITE(*,*) "********************************************************************************************"

    DO i=1,num_sites,1
      exch_field(i)=0d0
      DO j=1,num_shells(i),1
        DO k=1,num_neigh(i,j),1
          exch_field(i)=exch_field(i)+(2d0*J_list(i,j)*mom(neigh_list(i,j,k))*1.3806488d-23/(g*g*9.27400968d-24))
        END DO
      END DO

      arg=g*9.27400968d-24*spin(i)*(exch_field(i)+B)/(1.3806488d-23*T)

      IF (debug .eqv. .true.) THEN
        WRITE(*,*) "Full arg site",i,"--",arg
      END IF

      IF (low_x_brill .eqv. .true.) THEN
        new_mom=g*spin(i)*bril_low(arg,spin(i))
      ELSE
        new_mom=g*spin(i)*bril(arg,spin(i))
      END IF

      diff=new_mom-mom(i)    

      IF (DABS(diff) .gt. max_diff) THEN
        max_diff=DABS(diff)
      END IF

      IF (debug .eqv. .true.) THEN
        WRITE(*,*) "New_mom,mom,diff,max_diff",new_mom,mom(i),diff,max_diff
      END IF

      mom(i)=mom(i)+damping*diff

      WRITE(*,*) "Site",i,"moment--",mom(i),", exchange field--",exch_field(i),"Tesla"
     
    END DO
    WRITE(*,*) "Max_diff--",max_diff
    IF (max_diff .lt. converge) THEN
      converged=.true.
      WRITE(*,*) "Max_diff converged to within tolerance",converge
    END IF 
  END DO

! Write to output file
  WRITE(11,*) "Temperature=",T
  DO z=1,num_sites,1
    WRITE(11,*) "Moment at site",z,"=",mom(z)
  END DO

END DO


! Calculate Helmholtz energy

END PROGRAM


!--------------------------------------------
DOUBLE PRECISION FUNCTION bril_low(x,S)
IMPLICIT none

DOUBLE PRECISION :: x,S

bril_low=(S+1d0)*x/(3d0*S)

END FUNCTION

!--------------------------------------------
DOUBLE PRECISION FUNCTION bril(x,S)
! Codes the Brillouin function for spin S
IMPLICIT none

DOUBLE PRECISION :: x,S,dcoth,bril1

bril1=(2d0*S+1d0)/(2d0*S)
bril1=bril1*dcoth(x*(2d0*S+1d0)/(2d0*S)) 
bril=bril1-(1d0/(2d0*S))*dcoth(x/(2d0*S))

END FUNCTION

!--------------------------------------------
DOUBLE PRECISION FUNCTION dcoth(x)
! Codes the double precision hyperbolic cotangent
IMPLICIT none

DOUBLE PRECISION :: x

dcoth=dcosh(x)/dsinh(x)

END FUNCTION

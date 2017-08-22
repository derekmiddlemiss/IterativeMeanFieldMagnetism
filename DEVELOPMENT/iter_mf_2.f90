PROGRAM iterative_mean_field

! A program to solve the MF equations iteratively for arbitrary spin lattices
! Can use the linearised or fully nonlinear Brillouin functions
! Calculates thermodynamic potentials
! Derek S. Middlemiss, Cambridge, 2013

IMPLICIT none

INTEGER :: i,j,k,l,z,num_sites,max_shell,max_neigh,count,Tcount,Tstop
LOGICAL :: debug,converged,low_x_brill,momkeep
DOUBLE PRECISION :: bril,bril_low,B,T,converge,damping,max_diff,max_diff_last,g
DOUBLE PRECISION :: new_mom,diff,arg,Tstart,Tstep,x,part,A,U
DOUBLE PRECISION, ALLOCATABLE :: mom(:),exch_field(:)
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

! Allocate arrays with dimension num_sites
ALLOCATE(mom(num_sites))
ALLOCATE(start_mom(num_sites))
ALLOCATE(spin(num_sites))
ALLOCATE(num_shells(num_sites))
ALLOCATE(exch_field(num_sites))

READ(10,*) max_shell,max_neigh

! Allocate higher dimensional arrays
ALLOCATE(neigh_list(num_sites,max_shell,max_neigh))
ALLOCATE(num_neigh(num_sites,max_shell))
ALLOCATE(J_list(num_sites,max_shell))

! Initialize arrays
mom = 0d0
start_mom = 0d0
spin = 0d0
num_shells = 0
exch_field = 0d0
neigh_list = 0
num_neigh = 0
J_list = 0d0

DO i=1,num_sites,1
  READ(10,*) spin(i), num_shells(i)

  IF (num_shells(i) .gt. max_shell) THEN
    WRITE(*,*) "Problem! Site",i,"num_shells exceeds max_shell. Program stopping..."
    STOP
  END IF

  IF (debug .eqv. .true.) THEN
    WRITE(*,*) "--------------------------------------------------------------------------------------------"
    WRITE(*,*) "Site--",i,"spin,num_shells--",spin(i),num_shells(i)
  END IF

  DO j=1,num_shells(i),1
    READ(10,*) num_neigh(i,j), (neigh_list(i,j,k),k=1,num_neigh(i,j),1), J_list(i,j)

    IF (num_neigh(i,j) .gt. max_neigh) THEN
      WRITE(*,*) "Problem! Site",i,"shell",j,"num_neigh exceeds max_neigh. Program stopping..."
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

! Open output file
OPEN(unit=11,file="output",status='replace')

! Begin temperature loop
DO Tcount=1,Tstop,1
  
  T=Tstart+(DBLE(Tcount-1)*Tstep)
  
  converged= .false.

  WRITE(*,*)
  WRITE(*,*) "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
  WRITE(*,*) "TTTT Temperature",T,"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
  WRITE(*,*) "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"

  IF (Tcount .gt. 1 .AND. momkeep .eqv. .true.) THEN
    WRITE(*,*) "Keepmom active, initialising from converged moments at last temperature step."
  ELSE
    mom=start_mom
  END IF

  WRITE(*,*) "Starting moments"
  DO i=1,num_sites,1
    WRITE(*,*) "Site",i,"=",mom(i)
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
         exch_field(i)=exch_field(i)+(2.97745818182066d0*J_list(i,j)*mom(neigh_list(i,j,k))/(g*g))
        END DO
      END DO

      arg=g*6.71713884081165d-1*spin(i)*(exch_field(i)+B)/T

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
        WRITE(*,*) "New_mom,mom,diff,max_diff --",new_mom,mom(i),diff,max_diff
      END IF

      mom(i)=mom(i)+damping*diff

      WRITE(*,*) "Site",i,"moment=",mom(i),", exchange field=",exch_field(i),"Tesla"
     
    END DO
    WRITE(*,*) "Max_diff=",max_diff
    IF (max_diff .lt. converge) THEN
      converged=.true.
      WRITE(*,*) "Max_diff converged to within tolerance",converge
    END IF 
  END DO

! Write to output file
  WRITE(11,*) "Temperature=",T,"External field=",B,"Iterations=",count
  DO z=1,num_sites,1
    WRITE(11,*) "Moment at site",z,"=",mom(z)
  END DO

! Calculate Helmholtz and internal energies
  A=0d0
  U=0d0
  DO z=1,num_sites,1
    x=g*9.27400968d-24*(exch_field(z)+B)/(1.3806488d-23*T)
    part=dsinh((2d0*spin(z)+1d0)*(x/2d0))/dsinh(x/2d0)
    IF (debug .eqv. .true.) THEN
      WRITE(*,*) "Site",z,"x=",x,"Partition fn=",part,"A contribution=",-1d0*T*1.3806488d-23*dlog(part)/1.602176565d-19,"eV"
      WRITE(*,*) "Site",z,"U contribution=",-1d0*mom(z)*9.27400968d-24*(exch_field(z)+B)/1.602176565d-19,"eV"
      WRITE(*,*) "Site",z,"S contribution=",((-1d0*mom(z)*9.27400968d-24*(exch_field(z)+B)/1.602176565d-19)-&
 &    (-1d0*T*1.3806488d-23*dlog(part)/1.602176565d-19))/T
    END IF
    A=A-(T*1.3806488d-23*dlog(part))
    U=U-mom(z)*9.27400968d-24*(exch_field(z)+B)
    
  END DO
  WRITE(11,*) "Helmholtz energy A=",A/1.602176565d-19,"eV"
  WRITE(11,*) "Internal energy  U=",U/1.602176565d-19,"eV"
  WRITE(11,*) "Entropy          S=",(U-A)/(T*1.602176565d-19),"eV/K"

END DO

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

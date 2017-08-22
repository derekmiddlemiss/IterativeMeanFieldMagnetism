PROGRAM iterative_mean_field

! A program to solve the MF equations iteratively for arbitrary spin lattices
! Can use the linearised or fully nonlinear Brillouin functions
! Calculates thermodynamic potentials
! Altered so that it applies updates to moments simultaneously at end of each SCF, if required
! Added stoichiometric weights
! All exchange interactions defined within a -2J SUM_i>j S_i . S_j spin Hamiltonian

! Derek S. Middlemiss, Cambridge, 2013

IMPLICIT none

INTEGER :: i,j,k,l,z,num_sites,max_shell,max_neigh,count,Tcount,Tstop
LOGICAL :: debug,detail,converged,low_x_brill,momkeep,simupdate
DOUBLE PRECISION :: bril,bril_low,B,T,converge,damping,max_diff,g
DOUBLE PRECISION :: new_mom,arg,Tstart,Tstep,x,part,A,U,mom_sum
DOUBLE PRECISION, ALLOCATABLE :: mom(:),exch_field(:),diff(:),weights(:),exch_weights(:)
DOUBLE PRECISION, ALLOCATABLE :: spin(:),start_mom(:),J_list(:,:) 
INTEGER, ALLOCATABLE :: neigh_list(:,:,:),num_shells(:),num_neigh(:,:)

! Read initial data from INPUT file
OPEN (unit=10,file='control.in',status='old')
READ(10,*) debug
READ(10,*) detail
READ(10,*) low_x_brill
READ(10,*) momkeep
READ(10,*) simupdate
READ(10,*) converge,damping
READ(10,*) B,Tstart,Tstop,Tstep,g
READ(10,*) num_sites

! Allocate arrays with dimension num_sites
ALLOCATE(mom(num_sites))
ALLOCATE(start_mom(num_sites))
ALLOCATE(spin(num_sites))
ALLOCATE(num_shells(num_sites))
ALLOCATE(exch_field(num_sites))
ALLOCATE(diff(num_sites))
ALLOCATE(weights(num_sites))
ALLOCATE(exch_weights(num_sites))

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
weights = 0d0
diff = 0d0

READ(10,*) (weights(i),i=1,num_sites,1)
READ(10,*) (exch_weights(i),i=1,num_sites,1)

DO i=1,num_sites,1
  READ(10,*) spin(i), num_shells(i)

  IF (num_shells(i) .gt. max_shell) THEN
    WRITE(*,*) "Problem! Site",i,"num_shells exceeds max_shell. Program stopping..."
    STOP
  END IF

  IF (debug .eqv. .true.) THEN
    WRITE(*,*) "----------------------------------------------------------------------------------------------------"
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

IF (MAXVAL(neigh_list) .gt. num_sites) THEN
  WRITE(*,*) "Problem! One of the neighbour numbers exceeds num_sites. Program stopping..."
  STOP
END IF

IF(debug .eqv. .true.) THEN
  WRITE(*,*) "----------------------------------------------------------------------------------------------------"
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

  IF (detail .eqv. .true.) THEN
    WRITE(*,*)
    WRITE(*,*) "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
    WRITE(*,*) "TTTT Temperature",T,"Kelvin TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
    WRITE(*,*) "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
  END IF

  IF (Tcount .gt. 1 .AND. momkeep .eqv. .true.) THEN
    IF (detail .eqv. .true.) THEN
      WRITE(*,*) "Keepmom active, initialising from converged moments at last temperature step."
    END IF
  ELSE
    mom=start_mom
  END IF

  IF (detail .eqv. .true.) THEN
    WRITE(*,*) "Starting moments"
    DO i=1,num_sites,1
      WRITE(*,*) "Site",i,"=",mom(i)
    END DO
  END IF
  
  count=0

! Begin main loop
  DO WHILE (converged .eqv. .false.)

    count=count+1
    max_diff=0d0

    IF (detail .eqv. .true.) THEN
      WRITE(*,*)
      WRITE(*,*) "****************************************************************************************************"
      WRITE(*,*) "****Iteration number", count,"*******************************************************************"
      WRITE(*,*) "****************************************************************************************************"
    END IF

    DO i=1,num_sites,1
      exch_field(i)=0d0
      DO j=1,num_shells(i),1
        DO k=1,num_neigh(i,j),1
         exch_field(i)=exch_field(i)+(2.97745818182066d0*J_list(i,j)*mom(neigh_list(i,j,k))*&
 &       exch_weights(neigh_list(i,j,k))/(g*g))
         exch_field(i)=exch_field(i)*exch_weights(i)
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

      diff(i)=(new_mom-mom(i))*damping    

      IF (simupdate .eqv. .false.) THEN
        mom(i)=mom(i)+diff(i)
        IF (detail .eqv. .true.) THEN
          WRITE(*,*) "Non-simult. update. Site",i,"moment=",mom(i),", exchange field=",exch_field(i),"Tesla"
        END IF
      END IF

      IF (DABS(diff(i)/damping) .gt. max_diff) THEN
        max_diff=DABS(diff(i))/damping
      END IF

      IF (debug .eqv. .true.) THEN
          WRITE(*,*) "New_mom,mom,(diff*damping),(max_diff/damping) --",new_mom,mom(i),diff(i),max_diff
      END IF
     
    END DO

    IF (simupdate .eqv. .true.) THEN
    DO i=1,num_sites,1
      mom(i)=mom(i)+diff(i)
      IF (detail .eqv. .true.) THEN
        WRITE(*,*) "Simult. update. Site",i,"moment=",mom(i),", exchange field=",exch_field(i),"Tesla"
      END IF
    END DO
    END IF

    IF (detail .eqv. .true.) THEN
      WRITE(*,*) "(Max_diff/damping)=",max_diff
    END IF
    IF (max_diff .lt. converge) THEN
      converged=.true.
      IF (detail .eqv. .true.) THEN
        WRITE(*,*) "(Max_diff/damping) converged to within tolerance",converge
      END IF
    END IF 
  END DO

! Write to output file
  WRITE(11,*) "Temperature=",T,"External field=",B,"SCF iterations=",count
  DO z=1,num_sites,1
    WRITE(11,*) "T=",T,"Moment at site",z,"=",mom(z),"Weight=",weights(z),"Mom/sat_mom=",mom(z)/(g*spin(z))
  END DO

! Calculate total moment, uncorrelated spin partition function, Helmholtz and internal energies, and entropy.
  A=0d0
  U=0d0
  mom_sum=0d0
  DO z=1,num_sites,1
    mom_sum=mom_sum+(mom(z)*weights(z))
    x=g*9.27400968d-24*(exch_field(z)+B)/(1.3806488d-23*T)
    part=dsinh((2d0*spin(z)+1d0)*(x/2d0))/dsinh(x/2d0)
    IF (debug .eqv. .true.) THEN
      WRITE(*,*) "Site",z,"x=",x,"Partition fn=",part,"A contribution=",(-1d0*T*1.3806488d-23*dlog(part)/1.602176565d-19)&
 &    *weights(z),"eV"
      WRITE(*,*) "Site",z,"U contribution=",(-1d0*mom(z)*9.27400968d-24*(exch_field(z)+B)/1.602176565d-19)*weights(z),"eV"
      WRITE(*,*) "Site",z,"S contribution=",((-1d0*mom(z)*9.27400968d-24*(exch_field(z)+B)/1.602176565d-19)-&
 &    (-1d0*T*1.3806488d-23*dlog(part)/1.602176565d-19))*weights(z)/T
    END IF
    A=A-((T*1.3806488d-23*dlog(part))*weights(z))
    U=U-(mom(z)*9.27400968d-24*(exch_field(z)+B)*weights(z))
    
  END DO
  WRITE(11,*) "T=",T,"Total moment     M=",mom_sum,"mu_B"
  WRITE(11,*) "T=",T,"Helmholtz energy A=",A/1.602176565d-19,"eV"
  WRITE(11,*) "T=",T,"Internal energy  U=",U/1.602176565d-19,"eV"
  WRITE(11,*) "T=",T,"Entropy          S=",(U-A)/(T*1.602176565d-19),"eV/K"

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

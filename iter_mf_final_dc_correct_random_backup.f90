PROGRAM iterative_mean_field

! A program to solve the MF equations iteratively for arbitrary spin lattices
! Can use the linearised or fully nonlinear Brillouin functions
! Calculates thermodynamic potentials
! Altered so that it applies updates to moments simultaneously at end of each SCF, if required
! Added stoichiometric weights
! All exchange interactions and exchange fields defined within a -2J SUM_i>j S_i . S_j spin Hamiltonian

! Derek S. Middlemiss, Cambridge, 2013

IMPLICIT none

INTEGER :: i,j,k,l,z,num_sites,max_shell,max_neigh,count,Tcount,Tstop,min_step,rand_iters,rand_count
INTEGER :: max_count,rand_count_found
INTEGER*4 :: time_array(3)
LOGICAL :: debug,detail,converged,low_x_brill,momkeep,simupdate,bril_warned,rand_init
DOUBLE PRECISION :: bril,bril_low,B,T,converge,damping,max_diff,g,A_min,A_diff,A_threshold
DOUBLE PRECISION :: new_mom,arg,Tstart,Tstep,x,part,A,U_exch,U_field,mom_sum
DOUBLE PRECISION, ALLOCATABLE :: mom(:),exch_field(:),diff(:),weights(:),exch_weights(:),mom_min(:)
DOUBLE PRECISION, ALLOCATABLE :: spin(:),start_mom(:),J_list(:,:),exch_field_min(:) 
INTEGER, ALLOCATABLE :: neigh_list(:,:,:),num_shells(:),num_neigh(:,:)

! Read initial data from INPUT file
OPEN (unit=10,file='control.in',status='old')
READ(10,*) debug
READ(10,*) detail
READ(10,*) low_x_brill
READ(10,*) momkeep
READ(10,*) simupdate
READ(10,*) rand_init

READ(10,*) converge,damping,min_step
READ(10,*) rand_iters,A_threshold
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
ALLOCATE(mom_min(num_sites))
ALLOCATE(exch_field_min(num_sites))

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
mom_min = 0d0
exch_field_min = 0d0

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

IF (rand_init .eqv. .false.) THEN
  IF (detail .eqv. .true.) THEN
    WRITE(*,*) "Reading starting moments."
  END IF
  DO i=1,num_sites,1
    READ(10,*) start_mom(i)
  END DO
END IF

! Open output file
OPEN(unit=11,file="output",status='replace')

DO z=1,num_sites,1
  WRITE(11,*) "Site weight at site",z,"=",weights(z),"Exch. weight=",exch_weights(z)
END DO

! Seed random number generator from system clock (different initialisation each run)
call itime(time_array) 
call srand(time_array(1)+time_array(2)+time_array(3))
call itime(time_array)
i=rand(time_array(1)+time_array(2)+time_array(3))

bril_warned = .false.

! Begin temperature loop
DO Tcount=1,Tstop,1
  
  T=Tstart+(DBLE(Tcount-1)*Tstep)
  
  IF (detail .eqv. .true.) THEN
    WRITE(*,*)
    WRITE(*,*) "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
    WRITE(*,*) "TTTT Temperature",T,"Kelvin TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
    WRITE(*,*) "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
  END IF

  A_min=0d0
  max_count=0

! Begin random initial moments loop
  DO rand_count=1,rand_iters,1

  IF ((detail .eqv. .true.) .AND. (rand_init .eqv. .true.)) THEN
    WRITE(*,*)
    WRITE(*,*) "-RI-RI-RI-RI-RI-RI-RI-RI-RI-RI-RI-RI-RI-RI-RI-RI-RI-RI-RI-RI-RI-RI-RI-RI-RI-RI-RI-RI-RI-RI-RI-RI-RI-"
  END IF

  IF (rand_init .eqv. .false.) THEN
    IF ((Tcount .gt. 1) .AND. (momkeep .eqv. .true.)) THEN
      IF (detail .eqv. .true.) THEN
        WRITE(*,*) "Keepmom active, initialising from converged moments at last temperature step."
      END IF
    ELSE
      mom=start_mom
    END IF
  ELSE
    IF ((Tcount .gt. 1) .AND. (rand_count .eq. 1)) THEN
      IF (detail .eqv. .true.) THEN
        WRITE(*,*) "Trying the lowest free energy configuration from the previous temperature."
      END IF
      mom=mom_min
    ELSE
      IF (detail .eqv. .true.) THEN
        WRITE(*,*) "Random starting moments."
      END IF
      DO i=1,num_sites,1
        mom(i)=(4d0*spin(i)*DBLE(rand(0))-2d0*spin(i))*exch_weights(i)
      END DO
    END IF
  END IF

  IF (detail .eqv. .true.) THEN
    WRITE(*,*) "Starting moments"
    DO i=1,num_sites,1
      WRITE(*,*) "Site",i,"=",mom(i)
    END DO
  END IF

  converged = .false.
  count=0

! Begin SCF loop
  DO WHILE (converged .eqv. .false.)

    count=count+1
    max_diff=0d0

    IF (detail .eqv. .true.) THEN
      WRITE(*,*)
      WRITE(*,*) "****************************************************************************************************"
      WRITE(*,*) "****SCF iteration number", count,"***************************************************************"
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

      IF ((arg .gt. 1d0) .AND. (low_x_brill .eqv. .true.) .AND. (bril_warned .eqv. .false.)) THEN
        WRITE(*,*) "Caution! At least somewhere in your run, the argument of the Brillouin function exceeded 1."
        WRITE(*,*) "Argument =",arg
        WRITE(*,*) "You are using the low x approximation to the Brillouin function."
        WRITE(*,*) "Consider re-running using the fully non-linear Brillouin function."
        bril_warned = .true.
      END IF

      IF (debug .eqv. .true.) THEN
        WRITE(*,*) "Full arg site",i,"--",arg
      END IF

      IF (low_x_brill .eqv. .true.) THEN
        new_mom=g*spin(i)*bril_low(arg,spin(i))
      ELSE
        new_mom=g*spin(i)*bril(arg,spin(i))
      END IF

      diff(i)=(new_mom-mom(i))*damping*exch_weights(i)    

      IF (debug .eqv. .true.) THEN
        WRITE(*,*) "New_mom,mom,(diff*damping*exch_weight) --",new_mom,mom(i),diff(i)
      END IF

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
        WRITE(*,*) "max_diff/damping --",max_diff
      END IF
     
    END DO

    IF (simupdate .eqv. .true.) THEN
      DO i=1,num_sites,1
        IF (detail .eqv. .true.) THEN
          WRITE(*,*) "Simult. update. Site",i,"moment=",mom(i),", exchange field=",exch_field(i),"Tesla"
        END IF
      END DO
    END IF

    IF (detail .eqv. .true.) THEN
      WRITE(*,*) "(Max_diff/damping)=",max_diff
    END IF
    IF ((max_diff .lt. converge) .AND. (count .ge. min_step)) THEN
      converged=.true.
      IF (detail .eqv. .true.) THEN
        WRITE(*,*) "(Max_diff/damping) converged to within tolerance",converge,"and n(iterations) >=",min_step
      END IF
    END IF 

! End SCF loop
  END DO

! Final update of exchange fields to be consistent with final moments
  DO z=i,num_sites,1
    exch_field(i)=0d0
    DO j=1,num_shells(i),1
      DO k=1,num_neigh(i,j),1
       exch_field(i)=exch_field(i)+(2.97745818182066d0*J_list(i,j)*mom(neigh_list(i,j,k))*&
 &     exch_weights(neigh_list(i,j,k))/(g*g))
       exch_field(i)=exch_field(i)*exch_weights(i)
      END DO
    END DO
  END DO

! Calculate total moment, uncorrelated spin partition function, Helmholtz and internal energies from this random iteration
  A=0d0
  U_exch=0d0
  U_field=0d0
  DO z=1,num_sites,1
    x=g*9.27400968d-24*(exch_field(z)+B)/(1.3806488d-23*T)
    part=dsinh((2d0*spin(z)+1d0)*(x/2d0))/dsinh(x/2d0)
    A=A-((T*1.3806488d-23*dlog(part))*weights(z))
    U_exch=U_exch-(mom(z)*9.27400968d-24*exch_field(z)*weights(z))
    U_field=U_field-(mom(z)*9.27400968d-24*B*weights(z))
  END DO

! Correct U_exch for double counting
  U_exch=U_exch/2d0
! Correct A for double counting
  A=A-U_exch

  IF (detail .eqv. .true.) THEN
    WRITE(*,*)
    IF (rand_init .eqv. .true.) THEN
      WRITE(*,*) "******** Results of random iteration",rand_count,"********"
    END IF
    WRITE(*,*) "Temperature           T=",T
    WRITE(*,*) "Internal energy  U_exch=",U_exch/1.602176565d-19
    WRITE(*,*) "Internal energy U_field=",U_field/1.602176565d-19
    WRITE(*,*) "Total internal energy U=",(U_exch+U_field)/1.602176565d-19
    WRITE(*,*) "Helmholtz energy      A=",A/1.602176565d-19
  END IF

! Check to see if we have a new low energy structure. If so, save free energy and moments
  IF ((A .lt. A_min) .AND. ((DABS(A-A_min)/1.602176565d-19) .gt. A_threshold)) THEN
    IF (detail .eqv. .true.) THEN
      WRITE(*,*) "New minimum free energy configuration found and saved."
    END IF
    A_diff=A-A_min
    A_min=A
    rand_count_found=rand_count
    DO z=1,num_sites,1
      mom_min(z)=mom(z)
      exch_field_min(z)=exch_field(z)
    END DO
  END IF

! Track maximum number of SCF iterations
  IF (count .gt. max_count) THEN
    max_count=count
  END IF

! End random initialisation loop
END DO

! Write to output file
  WRITE(11,*) 
  WRITE(11,*) "Temperature=",T,"External field=",B,"Max SCF iterations=",max_count
  IF (rand_init .eqv. .true.) THEN
    WRITE(11,*) "Lowest energy configuration found at random iteration",rand_count_found,"of",rand_iters,&
 & "with A diff (eV)",A_diff/1.602176565d-19
  END IF
  DO z=1,num_sites,1
    WRITE(11,*) "T=",T,"Moment at site",z,"=",mom_min(z),"Mom/sat_mom=",mom_min(z)/(g*spin(z)),"Exch. field=",exch_field_min(z)
  END DO
  
  A=0d0
  U_exch=0d0
  U_field=0d0
  mom_sum=0d0
  DO z=1,num_sites,1
    mom_sum=mom_sum+(mom_min(z)*weights(z))
    x=g*9.27400968d-24*(exch_field_min(z)+B)/(1.3806488d-23*T)
    part=dsinh((2d0*spin(z)+1d0)*(x/2d0))/dsinh(x/2d0)
    A=A-((T*1.3806488d-23*dlog(part))*weights(z))
    U_exch=U_exch-(mom_min(z)*9.27400968d-24*exch_field_min(z)*weights(z))
    U_field=U_field-(mom_min(z)*9.27400968d-24*B*weights(z))
  END DO

! Correct U_exch for double counting
  U_exch=U_exch/2d0
! Correct A for double counting
  A=A-U_exch

  WRITE(11,*) "T=",T,"Total moment          M=",mom_sum,"mu_B"
  WRITE(11,*) "T=",T,"Helmholtz energy      A=",A/1.602176565d-19,"eV"
  WRITE(11,*) "T=",T,"Internal energy  U_exch=",U_exch/1.602176565d-19,"eV"
  WRITE(11,*) "T=",T,"Internal energy U_field=",U_field/1.602176565d-19,"eV"
  WRITE(11,*) "T=",T,"Total internal energy U=",(U_exch+U_field)/1.602176565d-19,"eV"
  WRITE(11,*) "T=",T,"Entropy               S=",(U_exch+U_field-A)/(T*1.602176565d-19),"eV/K"

!End temperature loop
END DO

CLOSE(unit=10)
CLOSE(unit=11)

END PROGRAM

!--------------------------------------------
DOUBLE PRECISION FUNCTION bril_low(x,S)
! Codes the low x approximation to the Brillouin function
! Caution, use only at high temperature and low field
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

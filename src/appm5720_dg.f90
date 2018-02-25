program dg
  use type_defs
  use quad_element
  use problemsetup
  use legendre_module
  implicit none
  ! list of all elements
  type(quad), dimension(:), allocatable :: qds
  integer :: i,j,iy,k,it,ix,l,ivar
  integer :: num_quads
  ! Quadrature
  real(kind=dp) :: weights(0:nint), xnodes(0:nint),diffmat(0:nint,0:nint)
  real(kind=dp) :: BFWeights(0:nint,2)
  ! Basis at quadrature points
  real(kind=dp), dimension(:,:),   allocatable :: P,DERP

  ! Weights for quadrature and differentiation on the elements.
  call lglnodes(xnodes,weights,nint)
  ! Differentiation matrix for the metric.
  do i = 0,nint
   call weights1(xnodes(i),xnodes,nint,nint,1,BFWEIGHTS)
   DiffMat(i,:) = BFWEIGHTS(:,2)
  end do
  ! Basis and derivative of basis at the quadrature nodes.
  allocate(P(0:nint,0:q),DERP(0:nint,0:q))
  do i = 0,nint
   do k = 0,q
    P(i,k)   = legendre(xnodes(i),k)
    DERP(i,k)  = derlegendre(xnodes(i),k)
   end do
  end do

  ! Set up a grid here
  ! Try not to assume anything about the connectivity of the
  ! elements
  ! Do allocation of all elements and fill in various information
  ! that is based on the grid.
  ! Vertices, neighbours, type of BC / internal face etc.

  ! Compute and store the metric on each quad.
  ! do i = 1,num_quads
  !   call set_metric(qds(i),xnodes,diffmat,nint)
  ! end do

  ! ! Rought estimate of the smallest h from the area of each cell
  ! h_eff = 1.d16
  ! do i = 1,num_quads
  !    do iy = 0,nint
  !       fint(iy) = sum(weights*qds(i)%jac(1:nint,iy))
  !    end do
  !    area = sum(weights*fint)
  !    if (sqrt(area) .lt. h_eff) h_eff = sqrt(area)
  ! end do

  ! Compute timestep, cfl may depend on q
  ! dt = cfl*h_eff
  ! nsteps = ceiling(tend/dt)
  ! dt = tend/dble(nsteps)
  ! write(*,*) 'The timestep and h_eff are: ',dt,h_eff

  ! Assemble mass and stiffness matrices
  ! then factor the mass matrices
  ! Also, set the initial data.
  !  do i = 1,num_quads
  !     call assemble(qds(i),nint,P,DERP,weights)
  !     m_size_u = (qds(i)%q+1)**2*qds(i)%nvar
  !     CALL DGETRF(m_size_u,m_size_u,qds(i)%M,&
  !          m_size_u,qds(i)%IPIV,INFO)
  !     if (info .ne. 0) write(*,*) "Warning LU factorization of M in element "&
  !          , i, " did not work. INFO = ",INFO
  !       stop 123
  !     end if
  !     call set_initial_data(qds(i),nint)
  !  end do

  !  it = 0
  !  call print_solution... we will talk more about this.
  !

  ! !
  ! ! Timestep Stuff
  ! !
  !  do it = 1,nsteps
  !     t = real(it-1,dp)*dt
  !     ! Now evolve the solution using RK4.
  !     do i = 1,num_quads
  !       do ivar = 1,nvar
  !           do l = 0,q
  !              do k = 0,q
  !                 um(k,l,ivar,i) = qds(i)%u(k,l,ivar)
  !              end do
  !           end do
  !        end do
  !     end do
  !    ! Set the interior states.
  !    do i = 1,num_quads
  !        call compute_my_face_states(qds(i),nint,P,DERP)
  !     end do
  !    ! This sets the exterior states based on adjacent elements
  !    ! or physical boundary conditions.
  !    call compute_outside_face_states(qds,num_quads,t)

  !     ! On exit the right hand sides will be computed.
  !     do i = 1,num_quads
  !        call compute_surface_integrals(qds(i),nint,chebmat,dchebmat,weights)
  !     end do
  !     ! Now we are ready to compute time-derivatives.
  !     do i = 1,num_quads
  !        call DGEMV ('N',m_size_u,m_size_u,-1.0_dp,qds(i)%S,&
  !             m_size_u,qds(i)%u,1,1.0_dp,qds(i)%fu,1)
  !        CALL DGETRS('N',m_size_u,1,qds(i)%M,m_size_u,&
  !             qds(i)%IPIV,qds(i)%fu,m_size_u,INFO)
  !        ! The time derivative is now in the right hand side array.
  !    end do
  !     ! Store stage 1 and prepare for stage 2.
  !     do i = 1,num_quads
  !        do ivar = 1,nvar
  !           do l = 0,q_u
  !              do k = 0,q_u
  !                 ku(k,l,ivar,1,i) = dt*qds(i)%fu(k,l,ivar)
  !                 qds(i)%u(k,l,ivar) = um(k,l,ivar,i) + 0.5_dp*ku(k,l,ivar,1,i)
  !              end do
  !           end do
  !        end do
  !     end do
  ! DEEA DO MORE STAGES IN RK4

  ! DEAA Sum up the stages.

  ! DEAA Print and possibly compute errors at various times

  ! end do
contains

!!$ subroutine compute_my_face_states(qd,nint,chebmat,dchebmat)
!!$  ! This routine computes and assigns the interior face states
!!$  ! for u
!!$  !
!!$  ! We assume the modes are ordered in column major order:
!!$  ! u_00, u_10, u_20,..., u_01,..., u_(q)(q).
!!$  ! Number of modes for u.
!!$  m_u  = qd%q
!!$  nvar = qd%nvar
!!$
!!$  qd%u_in = 0.d0
!!$
!!$  do ivar = 1,nvar
!!$     ! Assign the u states
!!$     do l = 0,m_u
!!$        do k = 0,m_u
!!$           mode_c = qd%u(k,l,ivar)
!!$           ! Face 1. s = 1, r \in [-1,1].
!!$           qd%u_in(:,1,ivar) = qd%u_in(:,1,ivar) &
!!$                + mode_c*P(0:nint,k)*P(nint,l)
!!$           ! Face 2. r = -1, s \in [-1,1].
!!$           qd%u_in(:,2,ivar) = qd%u_in(:,2,ivar) &
!!$                + mode_c*P(1,k)*P(0:nint,l)
!!$           ! Face 3. s = -1, r \in [-1,1].
!!$           qd%u_in(:,3,ivar) = qd%u_in(:,3,ivar) &
!!$                + mode_c*P(0:nint,k)*P(1,l)
!!$           ! Face 4. r = 1, s \in [-1,1].
!!$           qd%u_in(:,4,ivar) = qd%u_in(:,4,ivar) &
!!$                + mode_c*P(nint,k)*P(0:nint,l)
!!$        end do
!!$     end do
!!$  end do
!!$end subroutine compute_my_face_states

!!$subroutine compute_outside_face_states(qds,num_quads,t)

  ! Now set the outside states.

  ! Interior boundary, supply outside states from
  ! the adjacent element. KEEP TRACK OF ORIENTATION!


  !
  ! Physical boundary conditions based on the
  ! types set in the problemsetup module.
  ! We need to compute u_out
  ! call compute_bc_forcing ... it is suitable to put this in the problemsetup module

!!$end subroutine compute_outside_face_states

!!$subroutine compute_surface_integrals(qd,nint,chebmat,dchebmat,weights)

  ! ! Compute numerical fluxes!
  ! ! Central part
  !    ustar = 0.5d0*(qd%u_in(:,j,ivar)+qd%u_out(:,j,ivar))
  ! ADD ON UPWINDING

  ! Compute the contribution from the integrals and add to qd%fu(k,l,ivar)

!!$end subroutine compute_surface_integrals


!!$subroutine set_initial_data(qd,nint)

!!! DEAA, from hwk

!!$end subroutine set_initial_data


!!! DEAA, this is left as is! Rewrite to fit with quadrature points 0:nint....

!!$subroutine set_metric(qd,xnodes,diffmat,nint)
!!$  use doublePrecision
!!$  use quad_element
!!$  use problemsetup, only: pis
!!$  implicit none
!!$  type(quad) :: qd
!!$  integer :: nint
!!$  integer :: ix,iy
!!$  real(kind=dp) :: xnodes(nint)
!!$  real(kind=dp) :: x_coord_elem(nint,nint) ,y_coord_elem(nint,nint),diffmat(nint,nint)
!!$  real(kind=dp) :: pi1(2),pi2(2),pi3(2),pi4(2),pi2_m(2),pi2_p(2),pi4_m(2),pi4_p(2)
!!$  real(kind=dp) :: xy_loc(2),xy_s(2),xy_e(2),eta,xi
!!$
!!$  ! Compute metric
!!$  ! We use a Gordon-Hall mapping
!!$  ! The soubroutine pis must contain the approproate information
!!$  ! for the parametrization of the curvilinear elements
!!$  !
!!$  xy_s = qd%xy(3,1:2)
!!$  xy_e = qd%xy(2,1:2)
!!$  eta = 1.d0
!!$  call pis(pi2_p,eta,xy_s,xy_e,qd%bc_type(2))
!!$  eta = -1.d0
!!$  call pis(pi2_m,eta,xy_s,xy_e,qd%bc_type(2))
!!$  !
!!$  xy_s = qd%xy(4,1:2)
!!$  xy_e = qd%xy(1,1:2)
!!$  eta = 1.d0
!!$  call pis(pi4_p,eta,xy_s,xy_e,qd%bc_type(4))
!!$  eta = -1.d0
!!$  call pis(pi4_m,eta,xy_s,xy_e,qd%bc_type(4))
!!$  !
!!$  do iy = 1,nint
!!$     eta = xnodes(iy)
!!$     !
!!$     xy_s = qd%xy(3,1:2)
!!$     xy_e = qd%xy(2,1:2)
!!$     call pis(pi2,eta,xy_s,xy_e,qd%bc_type(2))
!!$     !
!!$     xy_s = qd%xy(4,1:2)
!!$     xy_e = qd%xy(1,1:2)
!!$     call pis(pi4,eta,xy_s,xy_e,qd%bc_type(4))
!!$     do ix = 1,nint
!!$        xi  = xnodes(ix)
!!$        !
!!$        xy_s = qd%xy(2,1:2)
!!$        xy_e = qd%xy(1,1:2)
!!$        call pis(pi1,xi,xy_s,xy_e,qd%bc_type(1))
!!$        !
!!$        xy_s = qd%xy(3,1:2)
!!$        xy_e = qd%xy(4,1:2)
!!$        call pis(pi3,xi,xy_s,xy_e,qd%bc_type(3))
!!$        xy_loc = (1.d0-eta)/2.d0*pi3+(1.d0+eta)/2.d0*pi1&
!!$             +(1.d0-xi)/2.d0*(pi2-(1.d0+eta)/2.d0*pi2_p-(1.d0-eta)/2.d0*pi2_m)&
!!$             +(1.d0+xi)/2.d0*(pi4-(1.d0+eta)/2.d0*pi4_p-(1.d0-eta)/2.d0*pi4_m)
!!$        x_coord_elem(ix,iy) = xy_loc(1)
!!$        y_coord_elem(ix,iy) = xy_loc(2)
!!$     end do
!!$  end do
!!$
!!$  qd%x = x_coord_elem
!!$  qd%y = y_coord_elem
!!$
!!$  call compute_curve_metric(qd%rx,qd%sx,qd%ry,qd%sy,qd%jac,&
!!$       x_coord_elem,y_coord_elem,Diffmat,nint)
!!$  ! Compute normals and line elements on all sides
  ! DEAA LEFT FOR YOU TO DO!!!
!!$
!!$  ! Face 1. corresponds to s = 1 and r \in [-1,1].
!!$  ! Thus the normal is (s_x,s_y) / \sqrt(s_x^2+s_y^2).
!!$  ! The line integral element is dl = \sqrt(x_r^2+y_r)^2 = J * \sqrt(s_x^2+s_y^2).
!!$  ! Compute the norm of the metric.
!!$  qd%dl_face(1:nint,1) = sqrt(qd%sx(1:nint,nint)**2+qd%sy(1:nint,nint)**2)
!!$  ! Compute outward pointing unit normal.
!!$  qd%nx_in(1:nint,1)   = qd%sx(1:nint,nint)/qd%dl_face(1:nint,1)
!!$  qd%ny_in(1:nint,1)   = qd%sy(1:nint,nint)/qd%dl_face(1:nint,1)
!!$  ! Scale by Jacobian to get the metric.
!!$  qd%dl_face(1:nint,1) = qd%dl_face(1:nint,1)*qd%jac(1:nint,nint)
!!$

  ! DEAA fill in the rest

!!$end subroutine set_metric

  subroutine compute_curve_metric(r_x,s_x,r_y,s_y,jac,X,Y,D,n)
    use type_defs
    implicit none
    integer :: n
    real(kind=dp), dimension(0:n,0:n) :: r_x,s_x,r_y,s_y,jac,X,Y,D
    real(kind=dp), dimension(0:n,0:n) :: x_r, x_s, y_r, y_s

    integer :: i
    !% Compute the derivatives w.r.t r & s
    do i = 1,n
     x_r(:,i) = matmul(D,X(:,i))
     y_r(:,i) = matmul(D,Y(:,i))
     x_s(i,:) = matmul(D,X(i,:))
     y_s(i,:) = matmul(D,Y(i,:))
    end do
    jac = x_r*y_s-y_r*x_s
    r_x =  y_s/jac
    r_y = -x_s/jac
    s_x = -y_r/jac
    s_y =  x_r/jac

  end subroutine compute_curve_metric

end program dg

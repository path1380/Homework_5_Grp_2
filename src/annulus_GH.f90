program annulus_GH
!This function attempts to numerically compute the area of the
!half annulus via a Gordon-Hall mapping from the reference 
!square [-1,1]^2 to a curvilinear quadrilateral.

  use type_defs
  use quad_element
  use problemsetup, only : q, nint, pis, nt, nr
  use legendre_module
  implicit none

  ! list of all elements
  type(quad), dimension(:), allocatable :: qds
  integer :: i,j,k,ind
  integer :: num_quads
  ! Quadrature
  real(kind=dp) :: weights(0:nint), xnodes(0:nint),diffmat(0:nint,0:nint)
  real(kind=dp) :: BFWeights(0:nint,2), pi, stepsize, r_start, r_end, theta_start, theta_end
  real(kind=dp), dimension(:), allocatable :: rvec, tvec
  real(kind=dp) :: area_approx, true_val

  ! Weights for quadrature and differentiation on the elements.
  call lglnodes(xnodes,weights,nint)

  ! Differentiation matrix for the metric.
  do i = 0,nint
   call weights1(xnodes(i),xnodes,nint,nint,1,BFWEIGHTS)
   DiffMat(i,:) = BFWEIGHTS(:,2)
  end do

  !find determine number of quads 
  num_quads = nt*nr

  !create an array of quads and vectors holding polar coordinates
  ALLOCATE(qds(num_quads))
  ALLOCATE(rvec(0:nr))
  ALLOCATE(tvec(0:nt))

  !What we'll do is pick an number of quads to use, then 
  !compute the area of the half annulus. Each quad will contain
  !its endpoints and then we'll give each side a number of quadrature
  !nodes and use the Gordon-Hall mapping.

  stepsize = 0.5_dp/dble(nr)
  do i = 0, nr
    rvec(i) = 0.5_dp + dble(i)*stepsize
  end do

  !precompute PI and true value of area
  pi = ACOS(-1.0_dp)
  true_val = 0.375_dp*pi

  stepsize = pi/dble(nt)
  do i=0,nt
    tvec(i) = dble(i)*stepsize
  end do

  !build each quad with its proper corners
  do i=0,nr-1
    r_start = rvec(i)
    r_end = rvec(i+1)
    do j =0,nt-1
      ind = j + nt*i + 1
      theta_start = tvec(j)
      theta_end = tvec(j+1)
      call allocate_quad(qds(ind),q,nint+1,2)
      !define corners of quad
      qds(ind)%xy(1,:) = (/r_end*COS(theta_start), r_end*SIN(theta_start)/)
      qds(ind)%xy(2,:) = (/r_end*COS(theta_end), r_end*SIN(theta_end)/)
      qds(ind)%xy(3,:) = (/r_start*COS(theta_end), r_start*SIN(theta_end)/)
      qds(ind)%xy(4,:) = (/r_start*COS(theta_start), r_start*SIN(theta_start)/)
      qds(ind)%my_ind = ind

      !Give appropriate boundary types for each element
      !note that this breaks if we only have one subdivision
      !along a direction
      if (i .eq. 0) then
        !these are the innermost elements, here 
        !only side 3 is curved
        qds(ind)%bc_type(1:2) = 10
        qds(ind)%bc_type(4) = 10
        qds(ind)%bc_type(3) = 12
      elseif(i .eq. nr-1) then
        !these are the outermost elements, here 
        !only side 1 is curved 
        qds(ind)%bc_type(2:4) = 10
        qds(ind)%bc_type(1) = 14
      else
        !these are the internal elements, here 
        !every side is straight
        qds(ind)%bc_type(:) = 10
      end if  
      end do
  end do

  !Compute and store the metric on each quad.
  do i = 1,num_quads
    call set_metric(qds(i),xnodes,diffmat,nint)
  end do


  !Now we approximate the area of the annulus
  area_approx = 0.0_dp
  do i = 1,num_quads
    do j = 0,nint
      do k = 0,nint
        area_approx = area_approx + qds(i)%jac(k+1,j+1)*weights(k)*weights(j)
      end do
    end do
  end do

  write(*,*) ABS(area_approx - true_val)

  do i =1,num_quads
    call deallocate_quad(qds(i))
  end do

contains

!!$subroutine set_initial_data(qd,nint)

!!! DEAA, from hwk

!!$end subroutine set_initial_data


!!! DEAA, this is left as is! Rewrite to fit with quadrature points 0:nint....

    subroutine set_metric(qd,xnodes,diffmat,nint)
      use type_defs
      use quad_element
      use problemsetup, only: pis
      implicit none
      type(quad) :: qd
      integer :: nint
      integer :: ix,iy
      real(kind=dp) :: xnodes(0:nint)
      real(kind=dp) :: x_coord_elem(0:nint,0:nint) ,y_coord_elem(0:nint,0:nint),diffmat(0:nint,0:nint)
      real(kind=dp) :: pi1(2),pi2(2),pi3(2),pi4(2),pi2_m(2),pi2_p(2),pi4_m(2),pi4_p(2)
      real(kind=dp) :: xy_loc(2),xy_s(2),xy_e(2),eta,xi

      ! Compute metric
      ! We use a Gordon-Hall mapping
      ! The soubroutine pis must contain the approproate information
      ! for the parametrization of the curvilinear elements
      !
      xy_s = qd%xy(3,1:2)
      xy_e = qd%xy(2,1:2)
      eta = 1.d0
      call pis(pi2_p,eta,xy_s,xy_e,qd%bc_type(2))
      eta = -1.d0
      call pis(pi2_m,eta,xy_s,xy_e,qd%bc_type(2))
      !
      xy_s = qd%xy(4,1:2)
      xy_e = qd%xy(1,1:2)
      eta = 1.d0
      call pis(pi4_p,eta,xy_s,xy_e,qd%bc_type(4))
      eta = -1.d0
      call pis(pi4_m,eta,xy_s,xy_e,qd%bc_type(4))
      !
      do iy = 0,nint
         eta = xnodes(iy)
         !
         xy_s = qd%xy(3,1:2)
         xy_e = qd%xy(2,1:2)
         call pis(pi2,eta,xy_s,xy_e,qd%bc_type(2))
         !
         xy_s = qd%xy(4,1:2)
         xy_e = qd%xy(1,1:2)
         call pis(pi4,eta,xy_s,xy_e,qd%bc_type(4))
         do ix = 0,nint
            xi  = xnodes(ix)
            !
            xy_s = qd%xy(2,1:2)
            xy_e = qd%xy(1,1:2)
            call pis(pi1,xi,xy_s,xy_e,qd%bc_type(1))
            !
            xy_s = qd%xy(3,1:2)
            xy_e = qd%xy(4,1:2)
            call pis(pi3,xi,xy_s,xy_e,qd%bc_type(3))
            xy_loc = (1.d0-eta)/2.d0*pi3+(1.d0+eta)/2.d0*pi1&
                 +(1.d0-xi)/2.d0*(pi2-(1.d0+eta)/2.d0*pi2_p-(1.d0-eta)/2.d0*pi2_m)&
                 +(1.d0+xi)/2.d0*(pi4-(1.d0+eta)/2.d0*pi4_p-(1.d0-eta)/2.d0*pi4_m)
            x_coord_elem(ix,iy) = xy_loc(1)
            y_coord_elem(ix,iy) = xy_loc(2)
         end do
      end do

      qd%x = x_coord_elem
      qd%y = y_coord_elem

      call compute_curve_metric(qd%rx,qd%sx,qd%ry,qd%sy,qd%jac,&
           x_coord_elem,y_coord_elem,Diffmat,nint)
      ! Compute normals and line elements on all sides

    !$
    !$  ! Face 1. corresponds to s = 1 and r \in [-1,1].
    !$  ! Thus the normal is (s_x,s_y) / \sqrt(s_x^2+s_y^2).
    !$  ! The line integral element is dl = \sqrt(x_r^2+y_r^2)| = J * \sqrt(s_x^2+s_y^2).
    !$  ! Compute the norm of the metric.
    !$  qd%dl_face(1:nint,1) = sqrt(qd%sx(1:nint,nint)**2+qd%sy(1:nint,nint)**2)
    !$  ! Compute outward pointing unit normal.
    !$  qd%nx_in(1:nint,1)   = qd%sx(1:nint,nint)/qd%dl_face(1:nint,1)
    !$  qd%ny_in(1:nint,1)   = qd%sy(1:nint,nint)/qd%dl_face(1:nint,1)
    !$  ! Scale by Jacobian to get the metric.
    !$  qd%dl_face(1:nint,1) = qd%dl_face(1:nint,1)*qd%jac(1:nint,nint)
    !$

    end subroutine set_metric

  subroutine compute_curve_metric(r_x,s_x,r_y,s_y,jac,X,Y,D,n)
    use type_defs
    implicit none
    integer :: n
    real(kind=dp), dimension(0:n,0:n) :: r_x,s_x,r_y,s_y,jac,X,Y,D
    real(kind=dp), dimension(0:n,0:n) :: x_r, x_s, y_r, y_s

    integer :: i
    !% Compute the derivatives w.r.t r & s
    do i = 0,n
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

end program annulus_GH

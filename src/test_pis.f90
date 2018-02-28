program test_pis
	use type_defs
	use problemsetup
	use legendre_module, only : lglnodes
	implicit none

	!Here we will test the pis routine to ensure
	!it is correct.

    real(kind=dp) :: xy(2),xy_start(2),xy_end(2),s
    integer, parameter :: n=7
    integer :: curve_type, i
    real(kind=dp) :: w(0:n),x(0:n)
    write(*,*) q

    ! curve_type = 12
    ! xy = 0.0_dp
    ! xy_start(:) = (/-0.5_dp, 0.0_dp/)
    ! xy_end(:) = (/0.5_dp, 0.0_dp/)

    curve_type = 14
    xy = 0.0_dp
    xy_start(:) = (/-1.0_dp, 0.0_dp/)
    xy_end(:) = (/1.0_dp, 0.0_dp/)

    !generate nodes and weights 
	call lglnodes(x,w,n)

    do i = 0, n
    	s = -1.0_dp + dble(i)*0.02_dp
    	call pis(xy,x(i),xy_start,xy_end,curve_type)
    	write(*,'(2(E24.16))') xy(1), xy(2)
    end do


end program test_pis

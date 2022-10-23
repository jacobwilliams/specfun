!
!>
!  Basic tests of all the routines.

    program specfun_tests

    use,intrinsic :: iso_fortran_env, only: wp => real64
    use specfun_module

    implicit none

    integer,parameter :: m  = 2 !! Order of Qmn(x)  ( m = 0,1,2,… )
    integer,parameter :: n  = 2 !! Degree of Qmn(x) ( n = 0,1,2,… )
    integer,parameter :: mm = 2 !! Physical dimension of QM and QD

    complex(wp) :: z, Cdn , Zd , Zf
    real(wp) :: Qm(0:Mm,0:n) !! Qmn(x)
    real(wp) :: Qd(0:Mm,0:n) !! Qmn'(x)
    real(wp) :: x

    interface check
      procedure :: check_scalar, check_real_matrix
    end interface

    !------------

    z = cmplx(1.0_wp, 2.0_wp, wp)
    call cpdsa(-2,z,Cdn)
    call check('cpdsa', 1, Cdn, (-0.38250608436459332_wp,-0.20581319133663531_wp))

    call cpdsa(0,z,Cdn)
    call check('cpdsa', 2, Cdn, (1.1438199904987183_wp,-1.7813940888174007_wp))

    z = cmplx(0.0_wp, 0.0_wp, wp)
    call cpdsa(0,z,Cdn)
    call check('cpdsa', 3, Cdn, (1.0_wp,0.0_wp))

    !------------

    z = cmplx(0.0_wp, 0.0_wp, wp)
    call cfs(z,Zf,Zd)
    call check('cfs', 1, Zf, (0.0_wp,0.0_wp))
    call check('cfs', 1, Zd, (0.0_wp,0.0_wp))

    z = cmplx(1.0_wp, 2.0_wp, wp)
    call cfs(z,Zf,Zd)
    call check('cfs', 2, Zf, (36.725464883991435_wp,15.587751104404592_wp)       )
    call check('cfs', 2, Zd, (267.74676148374817_wp,-4.91839391213067031E-014_wp))

    z = cmplx(0.1_wp, 0.1_wp, wp)
    call cfs(z,Zf,Zd)
    call check('cfs', 3, Zf, (-1.04727137798325495E-003_wp,1.04727137798325495E-003_wp))
    call check('cfs', 3, Zd, (0.0_wp,3.14210945037003767E-002_wp))

    z = cmplx(2.5_wp, 0.1_wp, wp)
    call cfs(z,Zf,Zd)
    call check('cfs', 4, Zf, (0.65749086637407972_wp,-4.17480398822913656E-002_wp))
    call check('cfs', 4, Zd, (-0.48762112846295241_wp,-0.80766983396363090_wp))

    z = cmplx(5.5_wp, 0.1_wp, wp)
    call cfs(z,Zf,Zd)
    call check('cfs', 5, Zf, (0.65579635983731077_wp,-5.93612315727645445E-002_wp))
    call check('cfs', 5, Zd, (-1.0687294488900367_wp,-2.5341170981639647_wp))


    ! real(wp),intent(in) :: x !! Argument of Qmn(x)
    ! integer,parameter :: m  = 2 !! Order of Qmn(x)  ( m = 0,1,2,… )
    ! integer,parameter :: n  = 2 !! Degree of Qmn(x) ( n = 0,1,2,… )
    ! integer,parameter :: mm = 2 !! Physical dimension of QM and QD
    ! real(wp),intent(out) :: Qm(0:Mm,0:n) !! Qmn(x)
    ! real(wp),intent(out) :: Qd(0:Mm,0:n) !! Qmn'(x)

    x = 4.0_wp
    call lqmn(mm,m,n,x,Qm,Qd)
    call check('lqmn', 1, Qm, &
                reshape( [0.25541281188299536_wp,      -0.25819888974716110_wp,       0.53333333333333321_wp,   &
                          2.1651247531981367E-002_wp,  -4.3585992157795711E-002_wp,  0.13333333333333330_wp,    &
                          2.2010792503905276E-003_wp,  -6.6341263992263184E-003_wp,   2.6909868068123705E-002_wp], [Mm+1,n+1]))
    call check('lqmn', 1, Qd, &
                reshape( [-6.6666666666666666E-002_wp,   6.8853037265909620E-002_wp, -0.15111111111111108_wp,   &
                          -1.1253854783671326E-002_wp,   2.2803587390875957E-002_wp,  -7.1111111111111097E-002_wp,&
                          -1.7129240707225676E-003_wp,   5.1789976853051055E-003_wp,  -2.1203625919222910E-002_wp], [Mm+1,n+1]))

    write(*,*) 'all tests passed'

    contains

        subroutine check_scalar(routine, test, z1, z2)
        !! checks if two compex numbers are the same to within a tolerance
        character(len=*),intent(in) :: routine
        integer,intent(in) :: test
        complex(wp),intent(in) :: z1, z2

        character(len=10) :: test_str

        real(wp),parameter :: tol = 1 * epsilon(1.0_wp) !! == tolerance

        if (abs(real(z1,wp)  - real(z2,wp)) > tol .or. &
            abs(cmplx(z1) - cmplx(z2)) > tol) then
            write(test_str,'(i10)') test
            write(*,*) '  z1: ', z1
            write(*,*) '  z2: ', z2
            error stop 'Error in '//trim(routine)//' test '//trim(adjustl(test_str))
        end if

        end subroutine check_scalar

        subroutine check_real_matrix(routine, test, z1, z2)
            !! checks if two real 2d arrays are the same to within a tolerance
            character(len=*),intent(in) :: routine
            integer,intent(in) :: test
            real(wp),dimension(:,:),intent(in) :: z1, z2

            character(len=10) :: test_str
            integer :: i !! counter

            real(wp),parameter :: tol = 1 * epsilon(1.0_wp) !! == tolerance

            do i = 1, size(z1,2)
                if (any(abs(z1(:,i) - z2(:,i)) > tol) .or. &
                    any(abs(z1(:,i) - z2(:,i)) > tol)) then
                    write(test_str,'(i10)') test
                    write(*,*) '  z1: ', z1
                    write(*,*) '  z2: ', z2
                    error stop 'Error in '//trim(routine)//' test '//trim(adjustl(test_str))
                end if
            end do

            end subroutine check_real_matrix

    end program specfun_tests
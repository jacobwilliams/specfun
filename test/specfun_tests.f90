!
!>
!  Basic tests of all the routines.

    program specfun_tests

    use,intrinsic :: iso_fortran_env, only: wp => real64
    use specfun_module

    implicit none

    complex(wp) :: Cdn

    call cpdsa(-2,cmplx(1.0_wp, 2.0_wp, wp),Cdn)
    if (Cdn /= (-0.38250608436459332_wp,-0.20581319133663531_wp)) error stop 'cpdsa error'

    call cpdsa(0,cmplx(1.0_wp, 2.0_wp, wp),Cdn)
    if (Cdn /= (1.1438199904987183_wp,-1.7813940888174007_wp)) error stop 'cpdsa error'

    call cpdsa(0,cmplx(0.0_wp, 0.0_wp, wp),Cdn)
    if (Cdn /= (1.0_wp,0.0_wp)) error stop 'cpdsa error'

    end program specfun_tests
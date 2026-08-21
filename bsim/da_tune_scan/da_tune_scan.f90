!+
! Program da_tune_scan
!
! Program for scanning dynamic aperture over a grid of betatron tunes.
! This combines the functionality of tune_scan and dynamic_aperture:
! - Scans over a grid of betatron frequencies (Qa, Qb, optionally Qz or pz)
! - At each tune point, calculates the dynamic aperture
! - Outputs the DA area (or other metrics) for each tune point
!-

program da_tune_scan_program

use da_ts_mod

implicit none

type (da_ts_params_struct) params
type (da_ts_com_struct) com
type (da_ts_data_struct), allocatable, target :: data(:,:,:)

real(rp) del_time, time0
integer ja, jb, jz

!---------------------------------------------
! Init

call da_ts_init_params (params, com)
allocate (data(0:com%n_a, 0:com%n_b, 0:com%n_z))
call run_timer ('START')
time0 = 0

!---------------------------------------------
! Main loop

do ja = 0, com%n_a
do jb = 0, com%n_b
do jz = 0, com%n_z
  call da_ts_calc_at_point (params, com, [ja, jb, jz], data(ja,jb,jz))
  call run_timer ('READ', del_time)

  if (del_time - time0 > params%timer_print_dtime) then
    print '(a, f10.2, a, 3i5)', 'Time (min): ', del_time/60, '  At point ja, jb, jz = ', ja, jb, jz
    time0 = del_time
  endif
enddo
enddo
enddo

!---------------------------------------------
! Write results

call da_ts_write_results (params, com, data)

end program da_tune_scan_program

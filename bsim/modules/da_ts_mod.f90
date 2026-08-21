module da_ts_mod

use bsim_interface
use mode3_mod
use dynamic_aperture_mod
 
implicit none

integer, parameter :: master_rank$   = 0
integer, parameter :: job_tag$       = 1000
integer, parameter :: have_data_tag$ = 1001
integer, parameter :: results_tag$   = 1002

type da_ts_params_struct
  character(100) :: lat_file = '', dat_out_file = '', quad_mask = ''
  character(40) :: group_knobs(2) = ["", ""]
  character(40) :: test = ""
  real(rp) :: Q_a0 = real_garbage$, Q_a1 = real_garbage$, dQ_a = real_garbage$
  real(rp) :: Q_b0 = real_garbage$, Q_b1 = real_garbage$, dQ_b = real_garbage$
  real(rp) :: Q_z0 = real_garbage$, Q_z1 = real_garbage$, dQ_z = real_garbage$
  real(rp) :: pz0 = real_garbage$, pz1 = real_garbage$, dpz = real_garbage$
  real(rp) :: a_emit = 0, b_emit = 0, sig_pz = 0
  real(rp) :: timer_print_dtime = 120
  integer :: n_turn = 0
  integer :: ix_branch = 0
  logical :: use_phase_trombone = .false.
  logical :: debug = .false.
  logical :: rf_on = .false.
  
  ! Dynamic aperture parameters
  real(rp) :: min_angle = 0
  real(rp) :: max_angle = pi
  integer :: n_angle = 37
  real(rp) :: x_init = 1e-3_rp
  real(rp) :: y_init = 1e-3_rp
  real(rp) :: rel_accuracy = 1e-2_rp
  real(rp) :: abs_accuracy = 1e-5_rp
  character(40) :: da_start_ele = ''
end type

type da_ts_com_struct
  character(100) master_input_file
  type (lat_struct) ring
  type (coord_struct), allocatable :: closed_orb(:)
  real(rp) sig_a, sig_b, sig_pz, a_emit, b_emit
  integer n_a, n_b, n_z
  integer int_Qa, int_Qb
  integer :: mpi_rank = master_rank$
  real(rp) :: time_start = 0
  logical :: using_mpi = .false.
end type

type da_ts_data_struct
  integer :: ix_q(3) = 0
  real(rp) :: tune(3) = 0   ! Either (Qa, Qb, Qz) with RF on or (Qa, Qb, pz) with RF off.
  real(rp) :: da_area = -1  ! Dynamic aperture area
  real(rp) :: da_x_max = -1 ! Maximum x aperture
  real(rp) :: da_y_max = -1 ! Maximum y aperture
  integer :: n_angle_survived = 0  ! Number of angles where particle survived all turns
  logical :: calc_successful = .false.
end type

contains

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine da_ts_init_params (ts, ts_com)

type (da_ts_params_struct), target :: ts
type (da_ts_com_struct), target :: ts_com
type (normal_modes_struct) mode
type (ele_struct), pointer :: ele

integer n_arg
logical err, rf_on

namelist / params / bmad_com, ts

!---------------------------------------------
! Read in the parameters

n_arg = command_argument_count()
if (n_arg > 1) then
  print '(a)', 'Usage: da_tune_scan <input_file>'
  print '(a)', 'Default: <input_file> = da_tune_scan.init'
  stop
endif

ts_com%master_input_file = 'da_tune_scan.init'
if (n_arg == 1) call get_command_argument(1, ts_com%master_input_file)

bmad_com%auto_bookkeeper = .false.

open (unit= 1, file = ts_com%master_input_file, status = 'old', action = 'read')
read(1, nml = params)
close (1)

if (ts%Q_a0 == real_garbage$ .or. ts%Q_a1 == real_garbage$ .or. ts%dQ_a == real_garbage$) then
  print '(a)', 'Error: One of ts%Q_a0, ts%Q_a1, ts%dQ_a is not set! Stopping here.'
  stop
endif

if (ts%Q_b0 == real_garbage$ .or. ts%Q_b1 == real_garbage$ .or. ts%dQ_b == real_garbage$) then
  print '(a)', 'Error: One of ts%Q_b0, ts%Q_b1, ts%dQ_b is not set! Stopping here.'
  stop
endif

if (ts%rf_on) then
  if (ts%Q_z0 == real_garbage$ .or. ts%Q_z1 == real_garbage$ .or. ts%dQ_z == real_garbage$) then
    print '(a)', 'Error: One of ts%Q_z0, ts%Q_z1, ts%dQ_z is not set with RF on! Stopping here.'
    stop
  endif
  ts%Q_z0 = abs(ts%Q_z0)
  ts%Q_z1 = abs(ts%Q_z1)
  ts%dQ_z = abs(ts%dQ_z)

else
  if (ts%pz0 == real_garbage$ .or. ts%pz1 == real_garbage$ .or. ts%dpz == real_garbage$) then
    print '(a)', 'Error: One of ts%pz0, ts%pz1, ts%dpz is not set with RF off! Stopping here.'
    stop
  endif
endif

if (ts%dat_out_file == '') call file_suffixer(ts_com%master_input_file, ts%dat_out_file, 'dat', .true.)

if (ts%group_knobs(1) == '' .neqv. ts%group_knobs(2) == '') then
  print '(a)', 'Error: Both ts%group_knobs(1) and ts%group_knobs(2) strings must be non-blank or both must be non-blank'
  stop
endif

if (ts%use_phase_trombone .and. ts%group_knobs(1) /= '') then
  print '(a)', 'Error: ts%use_phase_trombone and ts%group_knobs cannot both be used at the same time.'
  stop
endif 

!---------------------------------------------
! Calculate number of steps to take from range and step size

ts_com%n_a = 0
ts_com%n_b = 0
ts_com%n_z = 0

if (ts%dQ_a > 0) ts_com%n_a = nint(abs((ts%Q_a1 - ts%Q_a0) / ts%dQ_a))
if (ts%dQ_b > 0) ts_com%n_b = nint(abs((ts%Q_b1 - ts%Q_b0) / ts%dQ_b))

if (ts%rf_on) then
  if (ts%dQ_z > 0) ts_com%n_z = nint(abs((ts%Q_z1 - ts%Q_z0) / ts%dQ_z))
else
  if (ts%dpz > 0)  ts_com%n_z = nint(abs((ts%pz1 - ts%pz0) / ts%dpz))
endif

if (ts_com%mpi_rank == master_rank$) then
  print '(3(a, f10.6), a, i0, a)', 'ts%dQ_a, ts%Q_a0, ts%Q_a1 = [', ts%dQ_a, ', ', ts%Q_a0, ', ', ts%Q_a1, '],   n_a = [0, ', ts_com%n_a, ']'
  print '(3(a, f10.6), a, i0, a)', 'ts%dQ_b, ts%Q_b0, ts%Q_b1 = [', ts%dQ_b, ', ', ts%Q_b0, ', ', ts%Q_b1, '],   n_b = [0, ', ts_com%n_b, ']'
  if (ts%rf_on) then
    print '(3(a, f10.6), a, i0, a)', 'ts%dQ_z, ts%Q_z0, ts%Q_z1 = [', ts%dQ_z, ', ', ts%Q_z0, ', ', ts%Q_z1, '],   n_z = [0, ', ts_com%n_z, ']'
  else
    print '(3(a, f10.6), a, i0, a)', 'ts%dpz, ts%pz0, ts%pz1    = [', ts%dpz,  ', ', ts%pz0,  ', ', ts%pz1,  '],   n_z = [0, ', ts_com%n_z, ']'
  endif
endif

!---------------------------------------------
! Initialize lattice

bmad_com%auto_bookkeeper = .false.
global_com%exit_on_error = .false.

call bmad_parser(ts%lat_file, ts_com%ring, err_flag = err)
if (err) stop

rf_on = rf_is_on(ts_com%ring%branch(0))
if (ts%rf_on .neqv. rf_on) then
  if (ts%rf_on) then
    call set_on_off(rfcavity$, ts_com%ring, on$)
  else
    call set_on_off(rfcavity$, ts_com%ring, off$)
  endif
endif

rf_on = rf_is_on(ts_com%ring%branch(0))
if (ts%rf_on .neqv. rf_on) then
  print '(a)', 'Cannot turn RF on. Will stop here.'
  stop
endif

if (ts%use_phase_trombone) call insert_phase_trombone(ts_com%ring%branch(0))

allocate(ts_com%closed_orb(0:ts_com%ring%n_ele_max))
bmad_com%aperture_limit_on = .true.

if (ts%rf_on) then
  call closed_orbit_calc(ts_com%ring, ts_com%closed_orb, 6)
else
  call closed_orbit_calc(ts_com%ring, ts_com%closed_orb, 4)
endif

call lat_make_mat6(ts_com%ring, -1, ts_com%closed_orb)
call twiss_at_start(ts_com%ring)
call twiss_propagate_all(ts_com%ring)
call radiation_integrals (ts_com%ring, ts_com%closed_orb, mode)
call calc_z_tune (ts_com%ring%branch(0))

!---------------------------------------------

ts_com%int_Qa = int(ts_com%ring%ele(ts_com%ring%n_ele_track)%a%phi / twopi)
ts_com%int_Qb = int(ts_com%ring%ele(ts_com%ring%n_ele_track)%b%phi / twopi)
ele => ts_com%ring%ele(0)

ts_com%a_emit = ts%a_emit
if (ts_com%a_emit == 0) ts_com%a_emit = mode%a%emittance
ts_com%sig_a = sqrt(ts_com%a_emit * ele%a%beta)

ts_com%b_emit = ts%b_emit
if (ts_com%b_emit == 0) ts_com%b_emit = mode%b%emittance
ts_com%sig_b = sqrt(ts_com%b_emit * ele%b%beta)

ts_com%sig_pz = ts%sig_pz
if (ts_com%sig_pz == 0) ts_com%sig_pz = mode%sigE_E

end subroutine da_ts_init_params

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

recursive subroutine da_ts_calc_dynamic_aperture(ts, ts_com, jtune, ts_dat)

type (da_ts_params_struct) ts
type (da_ts_com_struct) ts_com
type (da_ts_data_struct) ts_dat
type (lat_struct), target :: ring
type (coord_struct), allocatable :: closed_orb(:)
type (ele_struct), pointer :: ele

real(rp) pz_val
integer jtune(3), iqm, status
logical ok, err

!
ts_dat = da_ts_data_struct()
ts_dat%ix_q = jtune

ts_dat%tune(1) = ts_com%int_Qa + ts%Q_a0 + jtune(1)*ts%dQ_a
ts_dat%tune(2) = ts_com%int_Qb + ts%Q_b0 + jtune(2)*ts%dQ_b
if (ts%rf_on) then
  ts_dat%tune(3) = -(ts%Q_z0 + jtune(3)*ts%dQ_z)
  iqm = 3
else
  ts_dat%tune(3) = ts%pz0 + jtune(3)*ts%dpz
  iqm = 2
endif

ring = ts_com%ring              ! Use copy in case tune setting fails
closed_orb = ts_com%closed_orb  ! Use copy in case tune setting fails
ele => ring%ele(0)

! Note: Tunes in radians.
ok = set_tune_3d (ring%branch(0), ts_dat%tune, ts%quad_mask, ts%use_phase_trombone, ts%rf_on, ts%group_knobs, .false.)
if (.not. ok) then
  print '(a, 3f9.4)', 'Note: Cannot set tunes (this is normal when close to a resonance) at: ', ts_dat%tune(1:iqm)
  return
endif

if (ts%rf_on) then
  call closed_orbit_calc(ring, closed_orb, 6)
else
  closed_orb(0)%vec(6) = ts_dat%tune(3)
  call closed_orbit_calc(ring, closed_orb, 4)
endif

call lat_make_mat6(ring, -1, closed_orb)
call twiss_at_start(ring, status, ts%ix_branch)
if (status /= ok$) then
  print '(a, 3f10.4)', 'Twiss calc fail at tunes:', ts_dat%tune/twopi
  return
endif

call twiss_propagate_all(ring)
call calc_z_tune (ring%branch(0))

!---------------------------------------------
! Calculate dynamic aperture at this tune point

call da_ts_compute_da_at_point(ring, closed_orb, ts, ts_dat)

end subroutine da_ts_calc_dynamic_aperture

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine da_ts_compute_da_at_point(ring, closed_orb, ts, ts_dat)

type (lat_struct), target :: ring
type (coord_struct), allocatable :: closed_orb(:)
type (da_ts_params_struct) ts
type (da_ts_data_struct) ts_dat

type (aperture_param_struct) ap_param
type (aperture_scan_struct), allocatable :: aperture_scan(:)
type (aperture_point_struct), pointer :: da_point
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele0

real(rp) pz_val(1), angle, x, y, area, x_max, y_max
integer i, n_survived
logical err

!
branch => ring%branch(0)
ele0 => branch%ele(0)

! Setup aperture parameters
ap_param%min_angle = ts%min_angle
ap_param%max_angle = ts%max_angle
ap_param%n_angle = ts%n_angle
ap_param%n_turn = ts%n_turn
ap_param%x_init = ts%x_init
ap_param%y_init = ts%y_init
ap_param%rel_accuracy = ts%rel_accuracy
ap_param%abs_accuracy = ts%abs_accuracy
if (ts%da_start_ele /= '') then
  ap_param%start_ele = ts%da_start_ele
else
  ap_param%start_ele = ''
endif

! Set pz value based on RF state
if (ts%rf_on) then
  pz_val(1) = 0.0_rp
else
  pz_val(1) = ts_dat%tune(3)
endif

! Perform dynamic aperture scan
call dynamic_aperture_scan(aperture_scan, ap_param, pz_val, ring, print_timing=.false.)

if (.not. allocated(aperture_scan) .or. size(aperture_scan) == 0) then
  print '(a)', 'Warning: Dynamic aperture scan returned no data'
  return
endif

! Calculate area and extract max values
area = 0.0_rp
x_max = 0.0_rp
y_max = 0.0_rp
n_survived = 0

do i = 1, ts%n_angle
  da_point => aperture_scan(1)%point(i)
  
  ! Track maximum x and y apertures
  x_max = max(x_max, abs(da_point%x))
  y_max = max(y_max, abs(da_point%y))
  
  ! Check if particle survived all turns
  if (da_point%i_turn >= ts%n_turn) then
    n_survived = n_survived + 1
  endif
  
  ! Calculate area using polar integration
  ! Area contribution from this angular segment
  if (i < ts%n_angle) then
    angle = ts%min_angle + (i-1) * (ts%max_angle - ts%min_angle) / (ts%n_angle - 1)
    area = area + 0.5_rp * (da_point%x**2 + da_point%y**2) * &
                  ((ts%max_angle - ts%min_angle) / (ts%n_angle - 1))
  endif
enddo

! Store results
ts_dat%da_area = area
ts_dat%da_x_max = x_max
ts_dat%da_y_max = y_max
ts_dat%n_angle_survived = n_survived
ts_dat%calc_successful = .true.

deallocate(aperture_scan)

end subroutine da_ts_compute_da_at_point

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine da_ts_write_results (ts, ts_com, ts_dat)

type (da_ts_params_struct) ts
type (da_ts_com_struct) ts_com
type (da_ts_data_struct), target :: ts_dat(0:,0:,0:)
type (da_ts_data_struct), pointer :: t

integer ja, jb, jz

!
open(unit = 23, file = ts%dat_out_file)

write (23, '(a, a)')         '# lat_file                   = ', quote(ts%lat_file)
write (23, '(a, a)')         '# quad_mask                  = ', quote(ts%quad_mask)
write (23, '(a, es12.4)')    '# Q_a0                       = ', ts%Q_a0
write (23, '(a, es12.4)')    '# Q_a1                       = ', ts%Q_a1
write (23, '(a, es12.4)')    '# dQ_a                       = ', ts%dQ_a
write (23, '(a, es12.4)')    '# Q_b0                       = ', ts%Q_b0
write (23, '(a, es12.4)')    '# Q_b1                       = ', ts%Q_b1
write (23, '(a, es12.4)')    '# dQ_b                       = ', ts%dQ_b
if (ts%rf_on) then
  write (23, '(a, es12.4)')  '# Q_z0                       = ', ts%Q_z0
  write (23, '(a, es12.4)')  '# Q_z1                       = ', ts%Q_z1
  write (23, '(a, es12.4)')  '# dQ_z                       = ', ts%dQ_z
else
  write (23, '(a, es12.4)')  '# pz0                        = ', ts%pz0
  write (23, '(a, es12.4)')  '# pz1                        = ', ts%pz1
  write (23, '(a, es12.4)')  '# dpz                        = ', ts%dpz
endif
write (23, '(a, i8)')        '# na_max                     = ', ts_com%n_a
write (23, '(a, i8)')        '# nb_max                     = ', ts_com%n_b
write (23, '(a, i8)')        '# nz_max                     = ', ts_com%n_z
write (23, '(a, i8)')        '# n_turn                     = ', ts%n_turn
write (23, '(a, i8)')        '# n_angle                    = ', ts%n_angle
write (23, '(a, es12.4, a)') '# sigma_a                    = ', ts_com%sig_a,  '  # Used in calculation'
write (23, '(a, es12.4, a)') '# sigma_b                    = ', ts_com%sig_b,  '  # Used in calculation'
write (23, '(a, es12.4, a)') '# sigma_pz                   = ', ts_com%sig_pz, '  # Used in calculation'
write (23, '(a, l4)')        '# radiation_damping_on       = ', bmad_com%radiation_damping_on
write (23, '(a, l4)')        '# radiation_fluctuations_on  = ', bmad_com%radiation_fluctuations_on
write (23, '(a, l4)')        '# rf_on                      = ', ts%rf_on
write (23, '(a, l4)')        '# use_phase_trombone         = ', ts%use_phase_trombone

if (ts%rf_on) then
  write (23, '(a, a4, 2a6, 3a10, a15, 3a15, a15)') '#-', 'ja', 'jb', 'jz', 'Q_a', 'Q_b', 'Q_z', &
                              'DA_area', 'DA_x_max', 'DA_y_max', 'n_angle_survived', 'calc_successful'
else
  write (23, '(a, a4, 2a6, 3a10, a15, 3a15, a15)') '#-', 'ja', 'jb', 'jz', 'Q_a', 'Q_b', 'pz', &
                              'DA_area', 'DA_x_max', 'DA_y_max', 'n_angle_survived', 'calc_successful'
endif

do jz = 0, ts_com%n_z
do jb = 0, ts_com%n_b
do ja = 0, ts_com%n_a
  t => ts_dat(ja, jb, jz)
  if (ts%rf_on) then
    write(23, '(3i6, 3f10.5, es15.6, 2es15.6, i15, a15)') ja, jb, jz, &
                  t%tune(1)-ts_com%int_Qa, t%tune(2)-ts_com%int_Qb, -t%tune(3), &
                  t%da_area, t%da_x_max, t%da_y_max, t%n_angle_survived, &
                  merge('T','F', t%calc_successful)
  else
    write(23, '(3i6, 3f10.5, es15.6, 2es15.6, i15, a15)') ja, jb, jz, &
                  t%tune(1)-ts_com%int_Qa, t%tune(2)-ts_com%int_Qb, t%tune(3), &
                  t%da_area, t%da_x_max, t%da_y_max, t%n_angle_survived, &
                  merge('T','F', t%calc_successful)
  endif
enddo
enddo
enddo

close(23)

end subroutine da_ts_write_results

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine da_ts_print_mpi_info (ts, ts_com, line, do_print)

type (da_ts_params_struct) ts
type (da_ts_com_struct) ts_com

real(rp) time_now
character(*) line
character(20) time_str
logical, optional :: do_print

!
if (.not. logic_option(ts%debug, do_print)) return

call run_timer ('ABS', time_now)
call date_and_time_stamp (time_str)
print '(a, f8.2, 2a, 2x, i0, 2a)', 'dTime:', (time_now-ts_com%time_start)/60, &
                                        ' Now: ', time_str, ts_com%mpi_rank, ': ', trim(line)

end subroutine da_ts_print_mpi_info

end module

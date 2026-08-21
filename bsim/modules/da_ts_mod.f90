!+
! Module da_ts_mod
!
! Module of routines used by the da_tune_scan program.
! This module combines functionality from ts_mod (tune_scan) and 
! dynamic_aperture_mod to scan dynamic aperture over a tune grid.
!-

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
  real(rp) :: da_min_angle = 0
  real(rp) :: da_max_angle = pi
  integer :: da_n_angle = 37
  real(rp) :: da_x_init = 1e-3_rp
  real(rp) :: da_y_init = 1e-3_rp
  real(rp) :: da_rel_accuracy = 1e-2_rp
  real(rp) :: da_abs_accuracy = 1e-5_rp
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
  real(rp) :: da_area = -1      ! Area of dynamic aperture
  real(rp) :: da_x_max = -1     ! Maximum horizontal aperture
  real(rp) :: da_y_max = -1     ! Maximum vertical aperture
  real(rp) :: da_x_min = -1     ! Minimum horizontal aperture  
  real(rp) :: da_y_min = -1     ! Minimum vertical aperture
  integer :: n_angle_survived = 0  ! Number of angles where particle survived all turns
  logical :: calc_successful = .false.
end type

contains

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine da_ts_init_params (params, com)

type (da_ts_params_struct), target :: params
type (da_ts_com_struct), target :: com
type (normal_modes_struct) mode
type (ele_struct), pointer :: ele

integer n_arg
logical err, rf_on

namelist / params / bmad_com, params

!---------------------------------------------
! Read in the parameters

n_arg = command_argument_count()
if (n_arg > 1) then
  print '(a)', 'Usage: da_tune_scan <input_file>'
  print '(a)', 'Default: <input_file> = da_tune_scan.init'
  stop
endif

com%master_input_file = 'da_tune_scan.init'
if (n_arg == 1) call get_command_argument(1, com%master_input_file)

bmad_com%auto_bookkeeper = .false.

open (unit= 1, file = com%master_input_file, status = 'old', action = 'read')
read(1, nml = params)
close (1)

if (params%Q_a0 == real_garbage$ .or. params%Q_a1 == real_garbage$ .or. params%dQ_a == real_garbage$) then
  print '(a)', 'Error: One of params%Q_a0, params%Q_a1, params%dQ_a is not set! Stopping here.'
  stop
endif

if (params%Q_b0 == real_garbage$ .or. params%Q_b1 == real_garbage$ .or. params%dQ_b == real_garbage$) then
  print '(a)', 'Error: One of params%Q_b0, params%Q_b1, params%dQ_b is not set! Stopping here.'
  stop
endif

if (params%rf_on) then
  if (params%Q_z0 == real_garbage$ .or. params%Q_z1 == real_garbage$ .or. params%dQ_z == real_garbage$) then
    print '(a)', 'Error: One of params%Q_z0, params%Q_z1, params%dQ_z is not set with RF on! Stopping here.'
    stop
  endif
  params%Q_z0 = abs(params%Q_z0)
  params%Q_z1 = abs(params%Q_z1)
  params%dQ_z = abs(params%dQ_z)

else
  if (params%pz0 == real_garbage$ .or. params%pz1 == real_garbage$ .or. params%dpz == real_garbage$) then
    print '(a)', 'Error: One of params%pz0, params%pz1, params%dpz is not set with RF off! Stopping here.'
    stop
  endif
endif

if (params%dat_out_file == '') call file_suffixer(com%master_input_file, params%dat_out_file, 'dat', .true.)

if (params%group_knobs(1) == '' .neqv. params%group_knobs(2) == '') then
  print '(a)', 'Error: Both params%group_knobs(1) and params%group_knobs(2) strings must be non-blank or both must be non-blank'
  stop
endif

if (params%use_phase_trombone .and. params%group_knobs(1) /= '') then
  print '(a)', 'Error: params%use_phase_trombone and params%group_knobs cannot both be used at the same time.'
  stop
endif 

!---------------------------------------------
! Calculate number of steps to take from range and step size

com%n_a = 0
com%n_b = 0
com%n_z = 0

if (params%dQ_a > 0) com%n_a = nint(abs((params%Q_a1 - params%Q_a0) / params%dQ_a))
if (params%dQ_b > 0) com%n_b = nint(abs((params%Q_b1 - params%Q_b0) / params%dQ_b))

if (params%rf_on) then
  if (params%dQ_z > 0) com%n_z = nint(abs((params%Q_z1 - params%Q_z0) / params%dQ_z))
else
  if (params%dpz > 0)  com%n_z = nint(abs((params%pz1 - params%pz0) / params%dpz))
endif

if (com%mpi_rank == master_rank$) then
  print '(3(a, f10.6), a, i0, a)', 'params%dQ_a, params%Q_a0, params%Q_a1 = [', params%dQ_a, ', ', params%Q_a0, ', ', params%Q_a1, '],   n_a = [0, ', com%n_a, ']'
  print '(3(a, f10.6), a, i0, a)', 'params%dQ_b, params%Q_b0, params%Q_b1 = [', params%dQ_b, ', ', params%Q_b0, ', ', params%Q_b1, '],   n_b = [0, ', com%n_b, ']'
  if (params%rf_on) then
    print '(3(a, f10.6), a, i0, a)', 'params%dQ_z, params%Q_z0, params%Q_z1 = [', params%dQ_z, ', ', params%Q_z0, ', ', params%Q_z1, '],   n_z = [0, ', com%n_z, ']'
  else
    print '(3(a, f10.6), a, i0, a)', 'params%dpz, params%pz0, params%pz1    = [', params%dpz,  ', ', params%pz0,  ', ', params%pz1,  '],   n_z = [0, ', com%n_z, ']'
  endif
  print '(a, i0, a)', 'Dynamic aperture: n_angle = ', params%da_n_angle
  print '(a, f10.6)', '  min_angle = ', params%da_min_angle
  print '(a, f10.6)', '  max_angle = ', params%da_max_angle
  print '(a, i0)', '  n_turn = ', params%n_turn
endif

!---------------------------------------------
! Initialize lattice

bmad_com%auto_bookkeeper = .false.
global_com%exit_on_error = .false.

call bmad_parser(params%lat_file, com%ring, err_flag = err)
if (err) stop

rf_on = rf_is_on(com%ring%branch(0))
if (params%rf_on .neqv. rf_on) then
  if (params%rf_on) then
    call set_on_off(rfcavity$, com%ring, on$)
  else
    call set_on_off(rfcavity$, com%ring, off$)
  endif
endif

rf_on = rf_is_on(com%ring%branch(0))
if (params%rf_on .neqv. rf_on) then
  print '(a)', 'Cannot turn RF on. Will stop here.'
  stop
endif

if (params%use_phase_trombone) call insert_phase_trombone(com%ring%branch(0))

allocate(com%closed_orb(0:com%ring%n_ele_max))
bmad_com%aperture_limit_on = .true.

if (params%rf_on) then
  call closed_orbit_calc(com%ring, com%closed_orb, 6)
else
  call closed_orbit_calc(com%ring, com%closed_orb, 4)
endif

call lat_make_mat6(com%ring, -1, com%closed_orb)
call twiss_at_start(com%ring)
call twiss_propagate_all(com%ring)
call radiation_integrals (com%ring, com%closed_orb, mode)
call calc_z_tune (com%ring%branch(0))

!---------------------------------------------

com%int_Qa = int(com%ring%ele(com%ring%n_ele_track)%a%phi / twopi)
com%int_Qb = int(com%ring%ele(com%ring%n_ele_track)%b%phi / twopi)
ele => com%ring%ele(0)

com%a_emit = params%a_emit
if (com%a_emit == 0) com%a_emit = mode%a%emittance
com%sig_a = sqrt(com%a_emit * ele%a%beta)

com%b_emit = params%b_emit
if (com%b_emit == 0) com%b_emit = mode%b%emittance
com%sig_b = sqrt(com%b_emit * ele%b%beta)

com%sig_pz = params%sig_pz
if (com%sig_pz == 0) com%sig_pz = mode%sigE_E

end subroutine da_ts_init_params

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine da_ts_calc_at_point(params, com, jtune, da_dat)

type (da_ts_params_struct) params
type (da_ts_com_struct) com
type (da_ts_data_struct) da_dat
type (lat_struct), target :: ring
type (coord_struct), allocatable :: closed_orb(:)
type (ele_struct), pointer :: ele
type (aperture_param_struct) :: ap_param
type (aperture_scan_struct), allocatable :: aperture_scan(:)
type (aperture_point_struct), pointer :: pt

real(rp) tune_vec(3)
integer jtune(3), iqm, i, j
logical ok, err

!

da_dat = da_ts_data_struct()
da_dat%ix_q = jtune

da_dat%tune(1) = com%int_Qa + params%Q_a0 + jtune(1)*params%dQ_a
da_dat%tune(2) = com%int_Qb + params%Q_b0 + jtune(2)*params%dQ_b
if (params%rf_on) then
  da_dat%tune(3) = -(params%Q_z0 + jtune(3)*params%dQ_z)
  iqm = 3
else
  da_dat%tune(3) = params%pz0 + jtune(3)*params%dpz
  iqm = 2
endif

ring = com%ring              ! Use copy in case tune setting fails
closed_orb = com%closed_orb  ! Use copy in case tune setting fails
ele => ring%ele(0)

! Note: Tunes in radians.
ok = set_tune_3d (ring%branch(0), da_dat%tune, params%quad_mask, params%use_phase_trombone, params%rf_on, params%group_knobs, .false.)
if (.not. ok) then
  print '(a, 3f9.4)', 'Note: Cannot set tunes (this is normal when close to a resonance) at: ', da_dat%tune(1:iqm)
  return
endif

! Recalculate closed orbit at new tunes

if (params%rf_on) then
  call closed_orbit_calc(ring, closed_orb, 6)
else
  closed_orb(0)%vec(6) = da_dat%tune(3)
  call closed_orbit_calc(ring, closed_orb, 4)
endif

call lat_make_mat6(ring, -1, closed_orb)
call twiss_at_start(ring, i, params%ix_branch)
if (i /= ok$) then
  print '(a, 3f10.4)', 'Twiss calc fail at tunes:', da_dat%tune/twopi
  return
endif

call twiss_propagate_all(ring)
call calc_z_tune (ring%branch(0))

!---------------------------------------------
! Calculate dynamic aperture at this tune point

! Set up aperture parameters
ap_param%min_angle = params%da_min_angle
ap_param%max_angle = params%da_max_angle
ap_param%n_angle = params%da_n_angle
ap_param%n_turn = params%n_turn
ap_param%x_init = params%da_x_init
ap_param%y_init = params%da_y_init
ap_param%rel_accuracy = params%da_rel_accuracy
ap_param%abs_accuracy = params%da_abs_accuracy
ap_param%start_ele = params%da_start_ele

! Use closed orbit as reference
! For DA calculation, we typically want to scan around the closed orbit at this tune point
! The dynamic_aperture_scan routine expects pz values for each scan point
! Since we're doing DA at fixed energy, we'll use a single pz value (0 for on-energy)

call dynamic_aperture_scan (aperture_scan, ap_param, [0.0_rp], ring, print_timing=params%debug)

if (.not. allocated(aperture_scan)) then
  print '(a)', 'Warning: dynamic_aperture_scan returned unallocated array'
  return
endif

if (size(aperture_scan) < 1) then
  print '(a)', 'Warning: dynamic_aperture_scan returned empty array'
  return
endif

! Calculate DA area and extract metrics
if (allocated(aperture_scan(1)%point)) then
  da_dat%calc_successful = .true.
  
  ! Count surviving angles and find extrema
  da_dat%n_angle_survived = 0
  da_dat%da_x_max = 0
  da_dat%da_y_max = 0
  da_dat%da_x_min = 0
  da_dat%da_y_min = 0
  
  do i = 1, params%da_n_angle
    pt => aperture_scan(1)%point(i)
    if (pt%i_turn >= params%n_turn) then
      da_dat%n_angle_survived = da_dat%n_angle_survived + 1
    endif
    
    ! Track maximum extents in each direction
    da_dat%da_x_max = max(da_dat%da_x_max, pt%x)
    da_dat%da_y_max = max(da_dat%da_y_max, pt%y)
    da_dat%da_x_min = min(da_dat%da_x_min, pt%x)
    da_dat%da_y_min = min(da_dat%da_y_min, pt%y)
  enddo
  
  ! Calculate approximate area using polygon area formula
  ! Area = 0.5 * |sum(x_i * y_{i+1} - x_{i+1} * y_i)|
  da_dat%da_area = 0.0_rp
  do i = 1, params%da_n_angle
    j = mod(i, params%da_n_angle) + 1
    pt => aperture_scan(1)%point(i)
    da_dat%da_area = da_dat%da_area + &
                     aperture_scan(1)%point(j)%x * pt%y - &
                     aperture_scan(1)%point(j)%y * pt%x
  enddo
  da_dat%da_area = 0.5_rp * abs(da_dat%da_area)
  
  if (params%debug) then
    print '(a, f12.6, a, f12.6, a, i4)', '  DA area: ', da_dat%da_area, &
                                          '  X: [', da_dat%da_x_min, ', ', da_dat%da_x_max, ']', &
                                          '  Y: [', da_dat%da_y_min, ', ', da_dat%da_y_max, ']', &
                                          '  Survived: ', da_dat%n_angle_survived
  endif
endif

end subroutine da_ts_calc_at_point

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine da_ts_write_results (params, com, da_dat)

type (da_ts_params_struct) params
type (da_ts_com_struct) com
type (da_ts_data_struct), target :: da_dat(0:,0:,0:)
type (da_ts_data_struct), pointer :: t

integer ja, jb, jz

!

open(unit = 23, file = params%dat_out_file)

write (23, '(a, a)')         '# lat_file                   = ', quote(params%lat_file)
write (23, '(a, a)')         '# quad_mask                  = ', quote(params%quad_mask)
write (23, '(a, es12.4)')    '# Q_a0                       = ', params%Q_a0
write (23, '(a, es12.4)')    '# Q_a1                       = ', params%Q_a1
write (23, '(a, es12.4)')    '# dQ_a                       = ', params%dQ_a
write (23, '(a, es12.4)')    '# Q_b0                       = ', params%Q_b0
write (23, '(a, es12.4)')    '# Q_b1                       = ', params%Q_b1
write (23, '(a, es12.4)')    '# dQ_b                       = ', params%dQ_b
if (params%rf_on) then
  write (23, '(a, es12.4)')  '# Q_z0                       = ', params%Q_z0
  write (23, '(a, es12.4)')  '# Q_z1                       = ', params%Q_z1
  write (23, '(a, es12.4)')  '# dQ_z                       = ', params%dQ_z
else
  write (23, '(a, es12.4)')  '# pz0                        = ', params%pz0
  write (23, '(a, es12.4)')  '# pz1                        = ', params%pz1
  write (23, '(a, es12.4)')  '# dpz                        = ', params%dpz
endif
write (23, '(a, i8)')        '# na_max                     = ', com%n_a
write (23, '(a, i8)')        '# nb_max                     = ', com%n_b
write (23, '(a, i8)')        '# nz_max                     = ', com%n_z
write (23, '(a, i8)')        '# n_turn                     = ', params%n_turn
write (23, '(a, i8)')        '# da_n_angle                 = ', params%da_n_angle
write (23, '(a, es12.4, a)') '# sigma_a                    = ', com%sig_a,  '  # Used in calculation'
write (23, '(a, es12.4, a)') '# sigma_b                    = ', com%sig_b,  '  # Used in calculation'
write (23, '(a, l4)')        '# rf_on                      = ', params%rf_on
write (23, '(a, l4)')        '# use_phase_trombone         = ', params%use_phase_trombone

if (params%rf_on) then
  write (23, '(a, a4, 2a6, 3a10, a12, 5a15, a8)') '#-', 'ja', 'jb', 'jz', 'Q_a', 'Q_b', 'Q_z', 'data_turns', &
                              'DA_area', 'DA_x_max', 'DA_y_max', 'DA_x_min', 'DA_y_min', 'n_survived'
else
  write (23, '(a, a4, 2a6, 3a10, a12, 5a15, a8)') '#-', 'ja', 'jb', 'jz', 'Q_a', 'Q_b', 'pz', 'data_turns', &
                              'DA_area', 'DA_x_max', 'DA_y_max', 'DA_x_min', 'DA_y_min', 'n_survived'
endif

do jz = 0, com%n_z
do jb = 0, com%n_b
do ja = 0, com%n_a
  t => da_dat(ja, jb, jz)
  if (params%rf_on) then
    write(23, '(3i6, 3f10.5, i12, f15.6, 4f15.6, i8)') ja, jb, jz, &
                  t%tune(1)-com%int_Qa, t%tune(2)-com%int_Qa, -t%tune(3), &
                  t%n_angle_survived, t%da_area, t%da_x_max, t%da_y_max, t%da_x_min, t%da_y_min, t%n_angle_survived
  else
    write(23, '(3i6, 3f10.5, i12, f15.6, 4f15.6, i8)') ja, jb, jz, &
                  t%tune(1)-com%int_Qa, t%tune(2)-com%int_Qa, t%tune(3), &
                  t%n_angle_survived, t%da_area, t%da_x_max, t%da_y_max, t%da_x_min, t%da_y_min, t%n_angle_survived
  endif
enddo
enddo
enddo

close(23)

end subroutine da_ts_write_results

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------

subroutine da_ts_print_mpi_info (params, com, line, do_print)

type (da_ts_params_struct) params
type (da_ts_com_struct) com

real(rp) time_now
character(*) line
character(20) time_str
logical, optional :: do_print

!

if (.not. logic_option(params%debug, do_print)) return

call run_timer ('ABS', time_now)
call date_and_time_stamp (time_str)
print '(a, f8.2, 2a, 2x, i0, 2a)', 'dTime:', (time_now-com%time_start)/60, &
                                        ' Now: ', time_str, com%mpi_rank, ': ', trim(line)

end subroutine da_ts_print_mpi_info

end module da_ts_mod

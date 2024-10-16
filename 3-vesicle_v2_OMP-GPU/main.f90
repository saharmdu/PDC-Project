program SLD

  use OMP_LIB

  use mod_param
  use mod_common
  use mod_zero
  use mod_core
  use mod_output

  implicit none

  integer :: istart,iend, t_nondi
  real :: t_simu
  
 
  !t_nondi = 100.  ! effective now, search for iend in zero.f90
  t_nondi = 1
  
  if (restart_from_last) then
  
  	open (unit=18, file=trim(restart_fname), status='old', &
         access='sequential', form='formatted', action='read')
    
    read(18,*)     ! skip the header
    read(18,'(I10)') istart
    
  	close(18)
  
  else
  	istart = 0
  	  
  endif

  !print *, "freq =", freq
  !print *, "shear_rate =", shear_rate
  
  if (oscillatory_shear) then
     t_simu = t_nondi*(2.*pi/freq)
  else
     t_simu = t_nondi/shear_rate
  endif
  !t_simu = 10

  call cpu_time(t_w0)
  t_w0 = omp_get_wtime() 

  call launcher(istart,t_simu,iend,t_nondi)  ! initialization

  call time_marching(istart,iend)    ! the main loop

  call cpu_time(t_w1)
  t_w1 = omp_get_wtime()

  call output_timing
  
  stop
end program SLD

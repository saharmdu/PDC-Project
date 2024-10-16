module mod_param

  implicit none

  !-- Constants --!
  
  real, parameter :: pi = acos(-1.)

  !-- Miscellaneous --!
    
  character(len=2),  parameter :: case_num = "1/"
  character(len=6),  parameter :: datadir = 'data'//trim(case_num)
  character(len=50), parameter :: init_fn = 'D3N500VF0.55Bidi1.4_0.5Cubic_1_.dat' 
  !logical,           parameter :: restart_from_last = .true.
  logical,           parameter :: restart_from_last = .false.
  !logical,           parameter :: start_oscl_from_stdy = .true.
  character(len=50), parameter :: restart_fname = 'restart_500/restart_latest.txt'
  !character(len=50), parameter :: restart_fname = 'data1/restart_sim/restart0020000k.txt'
  !character(len=50), parameter :: restart_fname = 'data5/restart0770000k.txt'
  character(len=1),  parameter :: dsub_num = "0"  ! index of sub folders in /data?, eg. /gap?
  character(len=1),  parameter :: dsub_part = "1"  ! index of real particle data folders in /data: para1?
  
  !-- Physical paramters --!
  
  !integer, parameter ::   np = 500      ! number of particles
  integer, parameter ::   np = 5
  
  real, parameter ::  a_1 = 1.          ! radius of group 1 sphere [length] (the smaller one)
  real, parameter ::  a_2 = 1.4         ! radius of group 2 sphere [length]
  real, parameter ::  rho = 1.          ! particle density [mass/length^3]
  real, parameter ::  k_n = 1e4         ! normal spring constant (collision) [mass/time^2]
  real, parameter :: dp_n = k_n*0.0     ! normal damping constant (collision) [mass/time]
  real, parameter ::  k_t = k_n*2./7.   ! tangential spring constant (collision) [mass/time^2]
  real, parameter :: mu_c = 0.5         ! sliding friction coefficient [1]
                                        
  real, parameter :: mu_f = 1.          ! fluid dynamic viscosity [mass/(time*length)]
  real, parameter :: shear_rate = 1e-2  ! [1/time] (ineffective if oscillatory_shear)

  logical, parameter :: oscillatory_shear = .false.
  real,    parameter :: amp0 = 1     ! oscillation amplitude (at z=Lz) normalized by Lz  [1]
  real,    parameter :: freq = 0.01     ! oscillation frequency [1/time]

  logical, parameter :: rotating_shear = .false.  ! second oscillatory shear with a -pi/2. shift
  
  logical, parameter :: elst_repl = .true.     ! electrostatic repulsion i.e. EDL force
  real,    parameter :: elst_sr   = 10.0        ! electrostatic repulsion scaled shear rate [1]
  real,    parameter :: Debye_len = 0.05*a_1   ! Debye length (appearing in the electrostatic repulsion model) [length]
  logical, parameter :: van_der_W = .false.     ! van der Waals (VdW) attraction
  real,    parameter :: C_Hamaker = 5e-4       ! Hamaker constant [energy]

  real, parameter :: fric_sr = 1e-2                 ! friction-scaled shear rate [1], redundant if fully_frictional

  logical, parameter :: frictional = .true.         ! frictional contact
  logical, parameter :: fully_frictional = .true.  ! always activate friction (i.e. zero critical loading)
  
  !-- Numerical parameters --!
  
  real, parameter :: h_lub_out = min(a_1,a_2)*0.05    ! lubrication cutoff (outer range)
  real, parameter :: h_lub_inn = min(a_1,a_2)*0.001   ! lubrication cutoff (inner range)

  integer, parameter :: integrate_method = 2
  ! 1: Forward Euler
  ! 2: Velocity-Verlet
  
end module mod_param

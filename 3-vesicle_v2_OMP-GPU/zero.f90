module mod_zero

  use mod_param
  use mod_common
  
  implicit none
  
  private
  public launcher

contains
  
  subroutine launcher(istart,t_simu,iend,t_nondi)

    real,    intent(in)  :: t_simu
    integer, intent(in)  :: istart, t_nondi
    integer, intent(out) :: iend

    real :: Stokes_num,stif_sr, T_period

    open( 17,file=trim(datadir)//'log.txt',position='append')
    
    write(17,'(A)') '****************************************************************'
    write(17,'(A)') ''
    write(17,'(A)') '  HLGD, a program for sheared lubrication dynamics simulations.  '
    write(17,'(A)') ''
    write(17,'(A)') '****************************************************************'
    write(17,'(A)') ''

    t9 = 0.
        
    call init_particles(istart)
    call init_parameters(istart)

    if (restart_from_last) then
       write(17,'(A)') 'Restarting the simulation from '//trim(restart_fname)
       write(17,'(A)') 
       call reinit_part_para ! overwriting part of the init_ files above
    endif

    if (oscillatory_shear) then
       Stokes_num = rho*(amp0*freq)*a_1**2/mu_f            ! particle inertia [1] (maximal)
       stif_sr = 2.*(amp0*freq)*a_1/sqrt(k_n/(2.*rho*a_1)) ! stiffness-scaled shear rate [1] (maximal)
    else
       Stokes_num = rho*shear_rate*a_1**2/mu_f            ! particle inertia [1]
       stif_sr = 2.*shear_rate*a_1/sqrt(k_n/(2.*rho*a_1)) ! stiffness-scaled shear rate [1]
    endif
    
    write(17,'(A,I4)')    'Number of particles   = ', np
    write(17,'(A,2F6.2)') 'Radii of particles    = ', a_1, a_2
    write(17,'(A,3F8.4)') 'Box size (Lx, Ly, Lz) = ', lx,ly,lz
    write(17,'(A,1F6.2)') 'Particle vol fraction = ', vol_frac
    write(17,'(A)') ''
    write(17,'(A)') '--- Particle parameters ---'
    write(17,'(A)') ''
    write(17,'(A,F7.0,A)') 'normal spring const     = ',k_n,  '  [mass/time^2]'
    write(17,'(A,F7.0,A)') 'normal damping coef     = ',dp_n, '  [mass/time]'
    if (frictional) then
       write(17,'(A,F7.0,A)') 'tangential spring const = ',k_t,  '  [mass/time^2]'
       write(17,'(A,F8.1,A)') 'sliding coef (Coulomb)  = ',mu_c, ' [1]'
    else
       write(17,'(A)') '(frictionless contact)'
    endif
    if (oscillatory_shear) then
       if (rotating_shear) then          
          write(17,'(A)') ''
          write(17,'(A)') '--- Rotating shear ---'
          write(17,'(A)') ''
          write(17,'(A,F9.2,A)' ) 'Rotation radius (A/Lz)  = ',amp0, '  [1]'
          write(17,'(A,ES9.3,A)') 'Angular frequency       = ',freq, '  [1/time]'
       else
          write(17,'(A)') ''
          write(17,'(A)') '--- Oscillatory shear ---'
          write(17,'(A)') ''
          write(17,'(A,F9.2,A)' ) 'Oscillation amplitude (A/Lz)  = ',amp0,  '  [1]'
          write(17,'(A,ES9.3,A)') 'Oscillation frequency         = ',freq, '  [1/time]'
       endif
    endif
    write(17,'(A)') ''
    write(17,'(A)') '--- Fluid properties ---'
    write(17,'(A)') ''
    write(17,'(A,F6.3,A)') 'dynamic viscosity = ',mu_f,       '  [mass/(time*length)]'
    if (.not. oscillatory_shear) then
       write(17,'(A,F6.3,A)') 'shear rate        = ',shear_rate, '  [1/time]'
    endif
    write(17,'(A)') ''
    write(17,'(A)') '--- Non-dimensional numbers ---'
    write(17,'(A)') ''
    write(17,'(A,ES9.3)') 'Stokes number (St)          = ',Stokes_num
    write(17,'(A,ES9.3)') 'Stiffness-scaled shear rate = ',stif_sr
    if (elst_repl) then
       write(17,'(A,ES9.3)') 'Electrostatic repulsion scaled shear rate = ',elst_sr
    endif
    if (van_der_W) then
       write(17,'(A,ES9.3)') 'Hamaker constant (van der Waals attraction) = ',C_Hamaker
    endif
    if (frictional) then
       if (fully_frictional) then
          write(17,'(A)') 'Friction-scaled shear rate  = infinity (always activate friction)'
       else
          write(17,'(A,ES9.3)') 'Friction-scaled shear rate  = ',fric_sr
       endif
    endif
    write(17,'(A)') ''

    call chkdt(17)
    
    iend = ceiling(t_simu/dt)
    !iend = 20*1000*1000 !900*1000*1000
    write(17,'(A,I10)')    'Total (i)steps = ',iend
    write(17,'(A,I10)') 'Total periods =',t_nondi 
    write(17,'(A)') ''
    close(17)

    ! output something every iout? steps

    if (oscillatory_shear) then
       T_period = 2.*pi/freq
    else
       T_period = 1./shear_rate
    endif
    
    iout0 = floor(T_period/10./dt) !20*1000
    iout1 = floor(T_period/1000./dt)
    iout2 = floor(T_period/10./dt)   ! effective (see data_output)
    iout3 = floor(T_period*5./dt)

    ! CPU timing
    t_nb_tot = 0.; t_comp_tot = 0.

    return
  end subroutine launcher

  !--------------------------------------------------------------

  subroutine init_particles(iter)

    use mod_output

    integer, intent(in) :: iter
    
    integer :: i,j, np1,np2, dummyi
    real :: dummyr
    character(len=1) :: dummyc

    if (.not. restart_from_last) then
       write(17,'(A)') 'Initialize particle positions from file: '//trim(init_fn)
       write(17,'(A)') '' ! unit 17 is opened in launcher
    endif

    open (unit=15, file='seeds/'//trim(init_fn), status='old', &
         access='sequential', form='formatted', action='read')

    ! read the headers of Romain Mari's file
    read(15,*) 
    read(15,*) dummyc,np1,np2,vol_frac,lx,ly,lz,dummyr,dummyr,dummyr

    vol_box = lx*ly*lz  ! box volume
    
    ! read particle positions from Mari's file
    j=1
    do i = 1,np
       
       read(15,*) pp(j),pp(j+1),pp(j+2), radii(i)
      

       if (radii(i) == a_1) then
          p_id(i) = 1
       elseif (radii(i) == a_2) then
          p_id(i) = 2
       else
          write(*,*) 'Error: unrecognized particle raidus. Program terminated.'
          stop
       endif

       ! init velocities and orientations
       pv(j)   = 0.
       pv(j+1) = shear_rate*(pp(j+2)-lz/2.) 
       pv(j+2) = 0.
       po(j)   = 0.
       po(j+1) = 0.
       po(j+2) = 0.
       pw(j)   = 0.
       pw(j+1) = 0.
       pw(j+2) = 0.

       j=j+3
       
    enddo
    
    ppr = pp		!stores real positions

    close(15)

    !call output_vtk(iter)

    return
  end subroutine init_particles

  !--------------------------------------------------------------

  subroutine init_parameters(iter)

    integer, intent(in) :: iter

    integer :: i
    real    :: lambda

    ! velocity shift (Lees-Edwards)
    
    if (oscillatory_shear) then
       y_shift = 0.*amp0*lz*sin(freq*t9)
       v_shift = (amp0*lz*freq)*cos(freq*t9)
       if (rotating_shear) then
          x_shift = amp0*lz*sin(freq*t9 -pi/2.)
          u_shift = (amp0*lz*freq)*cos(freq*t9 -pi/2.)
       endif
    else
       y_shift = 0.
       v_shift = shear_rate*lz

	x_shift = 0.
	u_shift = 0.  
    end if
    
    ! Stokes drag/torque coefficients
    F_stk_factor(1) = 6.*pi*mu_f*a_1; T_stk_factor(1) = 8.*pi*mu_f*a_1**3
    F_stk_factor(2) = 6.*pi*mu_f*a_2; T_stk_factor(2) = 8.*pi*mu_f*a_2**3

    masses(1) = rho*4./3.*pi*a_1**3
    masses(2) = rho*4./3.*pi*a_2**3

    mois(1) = 2./5.*masses(1)*a_1**2
    mois(2) = 2./5.*masses(2)*a_2**2

    ! Lubrication coefficients
    do i = -1,1

       lambda = (a_2/a_1)**i

       R_xiia1(i) = ( 6.*pi*mu_f)*(2.*lambda**2)/(1.+lambda)**3
       R_xiia2(i) = ( 6.*pi*mu_f)*lambda*(1.+7.*lambda+lambda**2)/(5.*(1.+lambda)**3)
       R_yiia(i)  = ( 6.*pi*mu_f)*4.*lambda*(2.+lambda+2.*lambda**2)/(15.*(1.+lambda)**3)
       R_yiib(i)  = (-4.*pi*mu_f)*lambda*(4.+lambda)/(5.*(1.+lambda)**2)
       R_yjib(i)  = (-4.*pi*mu_f)*(1./lambda)*(4.+1./lambda)/(5.*(1.+1./lambda)**2)
       R_yiic(i)  = ( 8.*pi*mu_f)*2.*lambda/(5.*(1.+lambda))
       R_yijc(i)  = ( 8.*pi*mu_f)*lambda**2/(10.*(1.+lambda))
       R_yjjc(i)  = ( 8.*pi*mu_f)*2.*(1./lambda)/(5.*(1.+1./lambda))
       
    enddo

    if (elst_repl) then
       if (oscillatory_shear) then
          F_elst = F_stk_factor(1)*a_1*(amp0*freq)/elst_sr
       else
          F_elst = F_stk_factor(1)*a_1*shear_rate/elst_sr
       endif
    endif


    if (fully_frictional) then
       F_crit = 0.
    else
       F_crit = F_stk_factor(1)*shear_rate*a_1/fric_sr
    endif
    
    if (iter .eq. 0) then
       fr_prev_mx(:,:)   = 0
       fr_curr_mx(:,:)   = 0
       xi_prev_mx(:,:,:) = 0.
       xi_curr_mx(:,:,:) = 0.
    endif
    
    return
  end subroutine init_parameters

  !-------------------------------------------------------------

  subroutine reinit_part_para
    
    integer :: i,j, np1,np2, dummyi

    open (unit=18, file=trim(restart_fname), status='old', &
         access='sequential', form='formatted', action='read')
    
    read(18,*)     ! skip the header
    read(18,*)		!skip the istart value, it has been read in main.f90
    
    !if (start_oscl_from_stdy) then
    !   read(18,*)     ! keep t9=0
    !else
    
    read(18,'(1ES16.8)') t9
    
    !endif
    
    read(18,'(4ES16.8)') vol_frac, lx,ly,lz
    read(18,'(2ES16.8)') y_shift, v_shift
    
    if(t9 .ne. 0 .AND. rotating_shear) then
        read(18,'(2ES16.8)') x_shift, u_shift
    endif       

    ! for t9 = 0: x_shit and u_shift = 0 from init_para

    vol_box = lx*ly*lz  ! box volume

    do i = 1,np
       read(18,'(I6,ES16.8,I3)') dummyi, radii(i), p_id(i)
    enddo

    j=1
    do i = 1,np
       read(18,'(I6,3ES16.8)') dummyi, pp(j),pp(j+1),pp(j+2) 
       j=j+3
    enddo
    
   
 ! add ppr data read
	
	if (t9 .eq. 0) then

		ppr =pp

	else    
		j=1
		do i = 1,np
   			read(18,'(I6,3ES16.8)') dummyi, ppr(j),ppr(j+1),ppr(j+2) 
			j=j+3
		enddo
	endif
	
!---!

    j=1
    do i = 1,np
       read(18,'(I6,3ES16.8)') dummyi, pv(j),pv(j+1),pv(j+2) 
       j=j+3
    enddo

    j=1
    do i = 1,np
       read(18,'(I6,3ES16.8)') dummyi, po(j),po(j+1),po(j+2) 
       j=j+3
    enddo

    j=1
    do i = 1,np
       read(18,'(I6,3ES16.8)') dummyi, pw(j),pw(j+1),pw(j+2) 
       j=j+3
    enddo

    do i = 1,np
       do j = 1,np
          read(18,'(2I6,I3,3ES16.8)') dummyi,dummyi, fr_prev_mx(j,i), xi_prev_mx(1,j,i),xi_prev_mx(2,j,i),xi_prev_mx(3,j,i)
       enddo
    enddo

    close(18)

    return
  end subroutine reinit_part_para

  !--------------------------------------------------------------

  subroutine chkdt(fnum)

    integer, intent(in) :: fnum

    real :: dt_1,dt_2, t_relax, max_F_net

    ! tangential collision relaxation time scale approximated by the lubrication
    t_relax = min(a_1,a_2)*minval(R_yiia)*log(max(a_1,a_2)/h_lub_inn)/k_t
    
    ! choose dt such that it is resolved at least by 10 time steps
    dt_1 = 0.1*t_relax

    ! maximal net forces approximated to be 10 times reference Stokes drag
    max_F_net = F_stk_factor(1)*shear_rate*a_1/masses(1)*10.

    ! choose dt such that the change of velocity is at most 0.01 times reference velocity
    dt_2 = 0.01*shear_rate*a_1/max_F_net

    !write(*,*) dt_1
    !write(*,*) dt_2
    
    dt = 2e-4 !0.5*min(dt_1,dt_2)
    
    write(fnum,'(A,ES9.3)') 'Time step (dt) = ',dt 
    
    return
  end subroutine chkdt

end module mod_zero

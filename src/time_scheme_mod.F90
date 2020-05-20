  module time_scheme_mod
    use stat_mod
    use tend_mod
    use mesh_mod
    use parameters_mod
    use spatial_operators_mod
    use ghost_mod
    use io_mod
    implicit none
    
    contains
    
    subroutine RK4(stat_new,stat_old)
      type(stat_field), intent(inout) :: stat_new
      type(stat_field), intent(in   ) :: stat_old

      integer, parameter :: old = 0
      integer, parameter :: new = 1
      
      integer, parameter :: one   = -1
      integer, parameter :: two   = -2
      integer, parameter :: three = -3
      integer, parameter :: four  = -4
      
      call spatial_operator (stat_old , tend(one))
      call update_stat      (stat(two), stat_old, tend(one), 0.5 * dt)
      !call history_write_tend(tend(one),2)
      !stop
      
      call spatial_operator (stat(two), tend(two))
      call update_stat      (stat(three), stat_old, tend(two), 0.5 * dt)
      
      call spatial_operator (stat(three), tend(three))
      call update_stat      (stat(four), stat_old, tend(three), dt)
      
      call spatial_operator(stat(four), tend(four))
      
      tend(new)%phiG = (tend(one)%phiG + 2. * tend(two)%phiG + 2. * tend(three)%phiG + tend(four)%phiG) / 6.
      tend(new)%u    = (tend(one)%u    + 2. * tend(two)%u    + 2. * tend(three)%u    + tend(four)%u   ) / 6.
      tend(new)%v    = (tend(one)%v    + 2. * tend(two)%v    + 2. * tend(three)%v    + tend(four)%v   ) / 6.
      
      call update_stat      (stat_new, stat_old, tend(new), dt)

    end subroutine RK4
    
    subroutine RK3_TVD(stat_new,stat_old)
      type(stat_field), intent(inout) :: stat_new
      type(stat_field), intent(in   ) :: stat_old

      integer, parameter :: old = 0
      integer, parameter :: new = 1
      
      integer, parameter :: one   = -1
      integer, parameter :: two   = -2
      integer, parameter :: three = -3
      integer, parameter :: four  = -4
      
      call spatial_operator (stat_old, tend(one))
      call update_stat      (stat(two), stat_old, tend(one), dt)
      
      call spatial_operator (stat(two), tend(two))
      call update_stat_RK3_TVD_1(stat(three), stat_old, stat(two), tend(two))
      
      call spatial_operator (stat(three), tend(three))
      call update_stat_RK3_TVD_2(stat_new, stat_old, stat(three), tend(three))

    end subroutine RK3_TVD
    
    subroutine update_stat(stat_new, stat_old, tend, inc_t)
      type(stat_field), intent(inout) :: stat_new
      type(stat_field), intent(in   ) :: stat_old
      type(tend_field), intent(in   ) :: tend
      real            , intent(in   ) :: inc_t
      
      stat_new%phiG = stat_old%phiG + inc_t * tend%phiG
      stat_new%u    = stat_old%u    + inc_t * tend%u   
      stat_new%v    = stat_old%v    + inc_t * tend%v   
      
      call correct_bdy_ghost(stat_new)
      
      stat_new%phi = stat_new%phiG / mesh%sqrtG
      
    end subroutine update_stat
    
    subroutine update_stat_RK3_TVD_1(stat_new, stat_old,stat1, tend)
      type(stat_field), intent(inout) :: stat_new
      type(stat_field), intent(in   ) :: stat_old
      type(stat_field), intent(in   ) :: stat1
      type(tend_field), intent(in   ) :: tend
      
      stat_new%phiG = 0.75 * stat_old%phiG + 0.25 * stat1%phiG + 0.25 * dt * tend%phiG
      stat_new%u    = 0.75 * stat_old%u    + 0.25 * stat1%u    + 0.25 * dt * tend%u   
      stat_new%v    = 0.75 * stat_old%v    + 0.25 * stat1%v    + 0.25 * dt * tend%v   
      
      call correct_bdy_ghost(stat_new)
      
      stat_new%phi = stat_new%phiG / mesh%sqrtG
      
    end subroutine update_stat_RK3_TVD_1
    
    subroutine update_stat_RK3_TVD_2(stat_new, stat_old,stat2, tend)
      type(stat_field), intent(inout) :: stat_new
      type(stat_field), intent(in   ) :: stat_old
      type(stat_field), intent(in   ) :: stat2
      type(tend_field), intent(in   ) :: tend
      
      stat_new%phiG = stat_old%phiG / 3. + 2./3. * stat2%phiG + 2./3. * dt * tend%phiG
      stat_new%u    = stat_old%u    / 3. + 2./3. * stat2%u    + 2./3. * dt * tend%u   
      stat_new%v    = stat_old%v    / 3. + 2./3. * stat2%v    + 2./3. * dt * tend%v  
      
      call correct_bdy_ghost(stat_new)
      
      stat_new%phi = stat_new%phiG / mesh%sqrtG
      
    end subroutine update_stat_RK3_TVD_2
    
    subroutine correct_bdy_ghost(stat)
      type(stat_field), intent(inout) :: stat
      
      call convert_wind_P2SP       (stat)
      call fill_ghost              (stat)
      call convert_wind_SP2P       (stat)
      call convert_wind_cov2contrav(stat)
      
    end subroutine correct_bdy_ghost
    
  end module time_scheme_mod
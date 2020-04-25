MODULE diag_mod
  use constants_mod
  use mesh_mod
  use parameters_mod
  use stat_mod
  
  implicit none
  
    contains
    subroutine calc_total_mass(total_mass,stat)
      type(stat_field), intent(in ) :: stat
      real            , intent(out) :: total_mass
      
      real    massOnCell(Nx,Ny,Nf)
      integer iCell,jCell,iPatch
      
      do iPatch = ifs, ife
        do jCell = 1, Ny
          do iCell = 1, Nx
            massOnCell(iCell,jCell,iPatch) = stat%phiG(iCell,jCell,iPatch)
            !massOnCell(iCell,jCell,iPatch) = mesh%areaCell(iCell,jCell,iPatch) * stat%phi(iCell,jCell,iPatch)
          enddo
        enddo
      enddo
      
      total_mass = sum(massOnCell)
    
    end subroutine calc_total_mass
    
    subroutine calc_total_energy(total_energy,stat)
      type(stat_field), intent(in ) :: stat
      real            , intent(out) :: total_energy
      
      real energyOnCell(Nx,Ny,Nf)
      real KE,PE
      integer iCell,jCell,iPatch
      
      energyOnCell = 0.
      total_energy = 0.
      do iPatch = ifs, ife
        do jCell = 1, Ny
          do iCell = 1, Nx
            KE = 0.5 * stat%phi(iCell,jCell,iPatch) * (stat%u(iCell,jCell,iPatch) * stat%uc(iCell,jCell,iPatch) + stat%v(iCell,jCell,iPatch) * stat%vc(iCell,jCell,iPatch))
            PE = 0.5 * (stat%phi(iCell,jCell,iPatch) + mesh%phis(iCell,jCell,iPatch))**2
            
            energyOnCell(iCell,jCell,iPatch) = mesh%sqrtG(iCell,jCell,iPatch) * (KE + PE)
            !energyOnCell(iCell,jCell,iPatch) = mesh%areaCell(iCell,jCell,iPatch) * (KE + PE)
          enddo
        enddo
      enddo
      
      total_energy = sum(energyOnCell)
    end subroutine calc_total_energy
END MODULE diag_mod


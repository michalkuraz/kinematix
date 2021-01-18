! Copyright 2008 Michal Kuraz, Petr Mayer, Copyright 2016  Michal Kuraz, Petr Mayer, Johanna Bloecher

! This file is part of DRUtES.
! DRUtES is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! DRUtES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with DRUtES. If not, see <http://www.gnu.org/licenses/>.

!> \file manage_pointers.f90
!! \brief Based on problem name definition calls appropriate pointer linker.
!<


module manage_pointers
  public :: set_pointers

  contains


    !> set pointers for the entire problem, except the file read pointers
    subroutine set_pointers()
      use typy
      use global_objs
      use pde_objs
      use globals
      use global4solver
      use core_tools
      use capmat
      use fem_tools
      use feminittools
      use fem
      use femmat
      use schwarz_dd
      use schwarz_dd2subcyc
      use solver_interfaces
      use debug_tools
      use bousspointers
      use re_pointers
      use ADE_pointers
      use Re_dual_pointers
      use heat_pointers
      use drutes_init
      use kinpointer
      use freeze_pointers
      use evappointers

      integer(kind=ikind) :: i, processes
      logical :: symetric
      
      

      select case(cut(drutes_config%name))
        case("REstd")
          write(unit=drutes_config%fullname, fmt=*) "Richards' equation, in", drutes_config%dimen, &
            "D, in pressure head form."
        
        
          call RE_processes(pde_common%processes)
          call pde_constructor(pde_common%processes)
          
          call RE_std(pde(1))
        case("RE")
          write(unit=drutes_config%fullname, fmt=*) "Richards' equation, in", drutes_config%dimen, &
            "D, in total hydraulic head form."
          call RE_processes(pde_common%processes)

          call pde_constructor(pde_common%processes)

          call REstdH(pde(1))
              
              
        case("boussi")   
           write(unit=drutes_config%fullname, fmt=*) " Boussinesq equation for hillslope runoff", &
                 "(1D approximation of groundwater flow)."
            
           call bouss_processes(pde_common%processes)
           call pde_constructor(pde_common%processes)      
           call boussi(pde(1))
                 
              
        case("ADE") 
        
          call ade_processes(pde_common%processes)
          call pde_constructor(pde_common%processes)
          write(unit=drutes_config%fullname, fmt=*) " advection-dispersion-reaction equation in", &
                 drutes_config%dimen, "D"
          call ade()

        case("Re_dual")
            write(unit=drutes_config%fullname, fmt=*) " Richards equation ", &
              "in total hydraulic head form for dual (fracture and matrix) medium"	
             pde_common%processes = 2
             call pde_constructor(pde_common%processes)
             call RE_matrix()
             call RE_fracture()  
      
    
        case("REtest")
          write(unit=drutes_config%fullname, fmt=*) "(debugs) itself"
          pde_common%processes = 3
          call pde_constructor(3_ikind)
          do i=1, 3
            call RE_std(pde(i))
          end do
	  
        case("heat")
        
          call heat_processes(pde_common%processes)
          call pde_constructor(pde_common%processes)
          write(unit=drutes_config%fullname, fmt=*) "heat conduction with convection"
          call heat(pde(:))
          
          
        case("kinwave")
        
          call kinwaveprocs(pde_common%processes)
          call pde_constructor(pde_common%processes)
          write(unit=drutes_config%fullname, fmt=*) "kinematic wave equation for real catchments"
          call kinwavelinker()
          
	  
        case("freeze")
        
          call freeze_processes(pde_common%processes)
          call pde_constructor(pde_common%processes)
          write(unit = drutes_config%fullname, fmt=*) "coupled water and heat flow considering freezing and melting"
          call frz_pointers()
          
       case("LTNE")
        
          call freeze_processes(pde_common%processes)
          call pde_constructor(pde_common%processes)
          write(unit = drutes_config%fullname, fmt=*) "local thermal non equilibrium"
          call frz_pointers()
          
        
           
        case("REevap")
        
            call REevap_proc(pde_common%processes)
            call pde_constructor(pde_common%processes)
            write(unit=drutes_config%fullname, fmt=*) "coupled heat and Richards equation with evaporation"
            call REevap_linker()
        

        case default
          print *, "your new model: ", trim(drutes_config%name), " requires pointer linking"
          print *, "exited from manage_pointers::set_pointers"
          ERROR stop

	  
      end select
      
      
      select case(cut(solver_name))
        case("LDU","LDUbalanced","LDUdefault")
          solve_matrix => LDU_face
        case("BJLDU")
          solve_matrix => blockjacobi_face
        case("PCGdiag","PCGbalanced")
          solve_matrix => CG_normal_face
      end select
!      select case(drutes_config%dimen)
!        case(1)
!          print *, cut(solver_name) ; stop
!!            solve_matrix => LDU_face
!            solve_matrix => blockjacobi_face
!            !solve_matrix => Minres_face
!!             solve_matrix => CG_normal_face
!        case(2)
!          symetric = .true.
!          do i=1, ubound(pde,1)
!            if (.not. pde(i)%symmetric) then
!              symetric = .false.
!              EXIT 
!            end if
!          end do
          
!          if (.not. symetric) then
!            solve_matrix => CG_normal_face
!          else
!            solve_matrix => cg_face
!          end if


!      end select
      
      select case (drutes_config%it_method)
        case(0) 
          pde_common%treat_pde => solve_picard
        case(1)
          pde_common%treat_pde => schwarz_picard
        case(2)
          pde_common%treat_pde => schwarz_subcyc
      end select
          
      select case(pde_common%timeint_method)
        case(0)
          pde_common%time_integ => steady_state_int
          solve_matrix => LDU_face
        case(1)
          pde_common%time_integ => impl_euler_np_diag
        case(2)
          pde_common%time_integ => impl_euler_np_nondiag
      end select


  end subroutine set_pointers


  

end module manage_pointers

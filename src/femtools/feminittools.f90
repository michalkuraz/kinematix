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



!> \file feminittools.f90
!! \brief FEM initialization
!<

!> Prepares integration nodes, reorders the nodes (Dirichlet boundary is excluded from unknowns), prepares surface integration.





module feminittools
  use typy
  use global_objs
  
  public :: feminit
  private :: init_integ, surface_integ, search_bnd_nodes, set_boundary
  private :: reorder

  contains
    
    subroutine feminit()
      use typy
      use integral
      use globals
      use global_objs
      use pde_objs
      use globals2d
      use geom_tools
      use re_constitutive
      use sparsematrix
      use debug_tools
      use printtools
      use core_tools
      use read_inputs

      integer(kind=ikind) :: i,j,k,l, locmatdim, process, limits,  el, nd
      real(kind=rkind), dimension(:,:), allocatable :: a_s
      integer(kind=ikind), dimension(:), allocatable :: itmp_vct
      
      
      locmatdim = (ubound(elements%data,2))*ubound(pde,1)
  
      allocate(cap_mat(locmatdim, locmatdim))
      
      allocate(stiff_mat(locmatdim, locmatdim))
      
      allocate(bside(locmatdim))

      allocate(elnode_prev(locmatdim))

      allocate(a_s(drutes_config%dimen+1, drutes_config%dimen+1))

      call init_integ()


      select case(drutes_config%dimen)
        case(1)
          do i=1, elements%kolik
            elements%areas(i) = abs(nodes%data(elements%data(i,2),1) - nodes%data(elements%data(i,1),1))
            elements%ders(i,1,1) = -1.0_rkind/(nodes%data(elements%data(i,2),1) - nodes%data(elements%data(i,1),1))
          end do

          elements%ders(:,2,1) = -elements%ders(:,1,1)  

        case(2)
          call write_log("calculating elements geometrical properties...")
          do i=1, elements%kolik
            call progressbar(int(100*i/elements%kolik))
            elements%areas(i) = triarea(nodes%data(elements%data(i,1),:), nodes%data(elements%data(i,2),:), &
              nodes%data(elements%data(i,3),:))
            do j=1,3
              a_s(j,1) = nodes%data(elements%data(i,j),1)
              a_s(j,2) = nodes%data(elements%data(i,j),2)
            end do

            a_s(:,3) = (/1.0_rkind, 0.0_rkind, 0.0_rkind/)
           
            
            call plane_derivative(a_s(1,:), a_s(2,:), a_s(3,:) , &
                  elements%ders(i,1,1), elements%ders(i,1,2), i)

            a_s(:,3) = (/0.0_rkind, 1.0_rkind, 0.0_rkind/)

            call plane_derivative(a_s(1,:), a_s(2,:), a_s(3,:) , &
                  elements%ders(i,2,1), elements%ders(i,2,2))

            a_s(:,3) = (/0.0_rkind, 0.0_rkind, 1.0_rkind/)

            call plane_derivative(a_s(1,:), a_s(2,:), a_s(3,:) , &
                  elements%ders(i,3,1), elements%ders(i,3,2))

          end do
      end select

      ! create for each node list of elements where the node belongs     
      do i=1, elements%kolik
        do j=1, ubound(elements%data,2)
          nd = elements%data(i,j)
          call nodes%element(nd)%fill(i)
        end do
      end do
      
      call reorder()

      i = ubound(pde,1)
      
      if (drutes_config%it_method /= 1) then
      	call spmatrix%init(maxval(pde(i)%permut(:)),maxval(pde(i)%permut(:)))
      end if

      allocate(pde_common%bvect(maxval(pde(i)%permut(:)))) 
      allocate(pde_common%xvect(maxval(pde(i)%permut(:)),4))
      

      ! fill nodes%el2integ
      !--------------
      call write_log("analyzing mesh structure...")
      do el=1, elements%kolik
        call progressbar(int(100*el/elements%kolik))
        do j=1, ubound(elements%data,2)
          nd=elements%data(el,j) 
          call nodes%el2integ(nd)%fill(el) 
        end do 
      end do
      !-----------------
      
       call surface_integ()
      
      if (.not. drutes_config%run_from_backup) then   
        do process=1, ubound(pde,1)
          call pde(process)%initcond()
        end do
      else
        call read_scilab(ADJUSTL(trim(backup_file)), 1_ikind)
      end if


      do process=1, ubound(pde,1)
        do j=1, nodes%kolik
          k = pde(process)%permut(j)
          if (k > 0) then
            pde_common%xvect(k,:) = pde(process)%solution(j)
          end if 
        end do
      end do



      deallocate(a_s)

    end subroutine feminit

    !> set ups integration nodes and weights
    subroutine init_integ()
      use typy
      use integral
      use globals2D
      use LinAlg
      use globals
      use geom_tools
      use pde_objs
      
      real(kind=rkind), dimension(:), allocatable :: uzly
      real(kind=rkind), dimension(:), allocatable :: vahy
      real(kind=rkind), dimension(3) :: bod
      real(kind=rkind), dimension(9,2) :: body
      real(kind=rkind), dimension(6) :: rezy
      real(kind=rkind), dimension(9,3,3) :: values
      real(kind=rkind), dimension(3) :: delka
      real(kind=rkind), dimension(3,2) :: ders
      integer(kind=ikind) :: l,j,i, no_points
      real(kind=rkind) :: tmp, tmp2
      
      !number of integration nodes, depends on integration method
      no_points = int(integ_method/10_ikind)

      allocate(gauss_points%point(no_points, drutes_config%dimen))

      allocate(gauss_points%weight(no_points))
      
      allocate(base_fnc(drutes_config%dimen+1, no_points))  

      allocate(integnode_now(no_points, ubound(pde,1)))

      allocate(integnode_prev(no_points, ubound(pde,1)))

      select case(drutes_config%dimen)
  
        case(1)

          if (modulo(no_points,2_ikind) == 0) then
            allocate(uzly(no_points/2))
            allocate(vahy(no_points/2))
          else
            allocate(uzly(no_points/2+1))
            allocate(vahy(no_points/2+1))
          end if

          call getform(no_points, uzly, vahy)


          if (modulo(no_points,2_ikind) == 0) then

            do i = 1, ubound(uzly,1) 
              gauss_points%point(i,1) = 0.5_rkind - uzly(i)/2.0_rkind
              gauss_points%point(i+ubound(uzly,1),1) = 0.5_rkind + uzly(i)/2.0_rkind
              gauss_points%weight(i) = vahy(i)/2.0_rkind
              gauss_points%weight(i+ubound(uzly,1)) = vahy(i)/2.0_rkind
            end do

          else

            gauss_points%point(ubound(uzly,1),1) = 0.5_rkind
            gauss_points%weight(ubound(uzly,1)) = vahy(ubound(vahy,1))/2.0_rkind

            do i = 1, ubound(uzly,1) - 1 
              gauss_points%point(i,1) = 0.5_rkind - uzly(i)/2.0_rkind
              gauss_points%point(i+ubound(uzly,1),1) = 0.5_rkind + uzly(i)/2.0_rkind
              gauss_points%weight(i) = vahy(i)/2.0_rkind
              gauss_points%weight(i+ubound(uzly,1)) = vahy(i)/2.0_rkind
            end do

          end if

          do i=1,no_points
            base_fnc(1,i) = 1 - gauss_points%point(i,1)
            base_fnc(2,i) = gauss_points%point(i,1)
          end do
          
          gauss_points%area = 1.0_rkind
          

        case(2)
          gauss_points%area = 0.5_rkind

          select case(integ_method)
            case(10)
              gauss_points%point(1,:) = (/0.333333333333333_rkind, 0.333333333333333_rkind/)

              gauss_points%weight = 0.5_rkind

            case(20)

              print *, "integration method", integ_method,   "for 2D not allowed, select from the following values"
              print *, "10, 30, 40, 60, 70, 90, 120"

            case(30) 
              gauss_points%point(1,:) = (/0.1666666666666667_rkind, 0.1666666666666667_rkind/)
              gauss_points%point(2,:) = (/0.6666666666666667_rkind, 0.1666666666666667_rkind/)
              gauss_points%point(3,:) = (/0.1666666666666667_rkind, 0.6666666666666667_rkind/)

              gauss_points%weight(1) = 0.1666666666666667_rkind
              gauss_points%weight(2) = 0.1666666666666667_rkind
              gauss_points%weight(3) = 0.1666666666666667_rkind

            case(40) 
              gauss_points%point(1,:) = (/0.2_rkind, 0.2_rkind/)
              gauss_points%point(2,:) = (/0.6_rkind, 0.2_rkind/)
              gauss_points%point(3,:) = (/0.2_rkind, 0.6_rkind/)
              gauss_points%point(4,:) = (/0.333333333333333_rkind,0.333333333333333_rkind/)

              gauss_points%weight(1) = 0.260416666666667_rkind
              gauss_points%weight(2) = 0.260416666666667_rkind
              gauss_points%weight(3) = 0.260416666666667_rkind
              gauss_points%weight(4) = -0.28125_rkind

            case(50)
              print *, "integration method", integ_method,   "for 2D not allowed, define one the following values"
              print *, "10, 30, 40, 60, 70, 90, 120"


            case(60)
              gauss_points%point(1,:) = (/0.816847572980459_rkind, 0.091576213509771_rkind/)
              gauss_points%point(2,:) = (/0.091576213509771_rkind, 0.816847572980459_rkind/)
              gauss_points%point(3,:) = (/0.091576213509771_rkind, 0.091576213509771_rkind/)
              gauss_points%point(4,:) = (/0.108103018168070_rkind, 0.445948490915965_rkind/)
              gauss_points%point(5,:) = (/0.445948490915965_rkind, 0.108103018168070_rkind/)
              gauss_points%point(6,:) = (/0.445948490915965_rkind, 0.445948490915965_rkind/)

              gauss_points%weight(1) = 0.054975871827661_rkind
              gauss_points%weight(2) = 0.054975871827661_rkind
              gauss_points%weight(3) = 0.054975871827661_rkind
              gauss_points%weight(4) = 0.111690794839005_rkind
              gauss_points%weight(5) = 0.111690794839005_rkind
              gauss_points%weight(6) = 0.111690794839005_rkind

            case(70)
              gauss_points%point(1,:) = (/0.333333333333333_rkind, 0.333333333333333_rkind/)
              gauss_points%point(2,:) = (/0.059715871798770_rkind, 0.470142064105115_rkind/)
              gauss_points%point(3,:) = (/0.470142064105115_rkind, 0.059715871798770_rkind/)
              gauss_points%point(4,:) = (/0.470142064105115_rkind, 0.470142064105115_rkind/)
              gauss_points%point(5,:) = (/0.797426985353087_rkind, 0.101286507323456_rkind/)
              gauss_points%point(6,:) = (/0.101286507323456_rkind, 0.797426985353087_rkind/)
              gauss_points%point(7,:) = (/0.101286507323456_rkind, 0.101286507323456_rkind/)

              gauss_points%weight(1) = 0.1125_rkind
              gauss_points%weight(2) = 0.066197076394253_rkind
              gauss_points%weight(3) = 0.066197076394253_rkind
              gauss_points%weight(4) = 0.066197076394253_rkind
              gauss_points%weight(5) = 0.066197076394253_rkind
              gauss_points%weight(6) = 0.066197076394253_rkind
              gauss_points%weight(7) = 0.066197076394253_rkind

            case(80)
                  print *, "integration method", integ_method,   "for 2D not allowed, define one the following values"
              print *, "10, 30, 40, 60, 70, 90, 120"


            case(90)
              !> warning this method differs from the others, it was derived form one dimensional problem transfered to quadrangle, that was transformed into triangle, Petr Mayer's trick :)
              allocate(uzly(2))
              allocate(vahy(2))

              call getform(3_ikind, uzly, vahy)
              bod(1) = (1 - uzly(1))/2
              bod(2) = 0.5_rkind
              bod(3) = (1 + uzly(1))/2

              
              delka(1) = 1 - bod(1)

              delka(2) = 0.5_rkind

              delka(3) = bod(1)  




                !shape of the unit triangle

              ! [0,1]
              ! |\
              ! |  \
              ! |     \
              ! |       \
              ! |         \
              ! |           \
              ! |               \ 
              ! |----7|---8|---9|---\ 
              ! |     |    |    |    \
              ! |----4|---5|---6|---  \
              ! |     |    |    |       \
              ! |----1|---2|---3|-----   -\       
              ! |     |    |    |           \ 
              ! |     |    |    |             \
              ! -----------------------------------------
              ! [0,0]                          [1,0]

              body(1,1) = bod(1)
              body(1,2) = bod(1)*delka(1)
              body(2,1) = bod(2)
              body(2,2) = bod(1)*delka(2)
              body(3,1) = bod(3)
              body(3,2) = bod(1)*delka(3)
              body(4,1) = bod(1)
              body(4,2) = bod(2)*delka(1)
              body(5,1) = bod(2)
              body(5,2) = bod(2)*delka(2)
              body(6,1) = bod(3)
              body(6,2) = bod(2)*delka(3)
              body(7,1) = bod(1)
              body(7,2) = bod(3)*delka(1)
              body(8,1) = bod(2)
              body(8,2) = bod(3)*delka(2)
              body(9,1) = bod(3)
              body(9,2) = bod(3)*delka(3)

            
              
              
              gauss_points%weight(1) = vahy(1)*vahy(1)/4*(1-body(1,1))
              gauss_points%weight(2) = vahy(2)*vahy(1)/4*(1-body(2,1))
              gauss_points%weight(3) = vahy(1)*vahy(1)/4*(1-body(3,1))
              gauss_points%weight(4) = vahy(2)*vahy(1)/4*(1-body(4,1))
              gauss_points%weight(5) = vahy(2)*vahy(2)/4*(1-body(5,1))
              gauss_points%weight(6) = vahy(2)*vahy(1)/4*(1-body(6,1))
              gauss_points%weight(7) = vahy(1)*vahy(1)/4*(1-body(7,1))
              gauss_points%weight(8) = vahy(2)*vahy(1)/4*(1-body(8,1))
              gauss_points%weight(9) = vahy(1)*vahy(1)/4*(1-body(9,1))
              
              gauss_points%point = body

              gauss_points%area = 0.5_rkind	

            case(120)
              gauss_points%point(1,:) = (/0.873821971016996_rkind, 0.063089014491502_rkind/)
              gauss_points%point(2,:) = (/0.063089014491502_rkind, 0.873821971016996_rkind/)
              gauss_points%point(3,:) = (/0.063089014491502_rkind, 0.063089014491502_rkind/)
              gauss_points%point(4,:) = (/0.501426509658179_rkind, 0.249286745170910_rkind/)
              gauss_points%point(5,:) = (/0.249286745170910_rkind, 0.501426509658179_rkind/)
              gauss_points%point(6,:) = (/0.249286745170910_rkind, 0.249286745170910_rkind/)
              gauss_points%point(7,:) = (/0.636502499121399_rkind, 0.310352451033785_rkind/)
              gauss_points%point(8,:) = (/0.310352451033785_rkind, 0.053145049844816_rkind/)
              gauss_points%point(9,:) = (/0.636502499121399_rkind, 0.053145049844816_rkind/)
              gauss_points%point(10,:) = (/0.636502499121399_rkind, 0.053145049844816_rkind/)
              gauss_points%point(11,:) = (/0.310352451033785_rkind, 0.636502499121399_rkind/)
              gauss_points%point(12,:) = (/0.053145049844816_rkind, 0.310352451033785_rkind/)

              gauss_points%weight(1) = 0.025422453185103_rkind
              gauss_points%weight(2) = 0.025422453185103_rkind
              gauss_points%weight(3) = 0.025422453185103_rkind
              gauss_points%weight(4) = 0.058393137863189_rkind
              gauss_points%weight(5) = 0.058393137863189_rkind
              gauss_points%weight(6) = 0.058393137863189_rkind
              gauss_points%weight(7) = 0.041425537809187_rkind
              gauss_points%weight(8) = 0.041425537809187_rkind
              gauss_points%weight(9) = 0.041425537809187_rkind
              gauss_points%weight(10) = 0.041425537809187_rkind
              gauss_points%weight(11) = 0.041425537809187_rkind
              gauss_points%weight(12) = 0.041425537809187_rkind


          end select

          call plane_derivative((/0.0_rkind, 0.0_rkind, 1.0_rkind/), (/1.0_rkind, 0.0_rkind, 0.0_rkind/), (/0.0_rkind, &
            1.0_rkind, 0.0_rkind/), ders(1,1), ders(1,2))

          call plane_derivative((/0.0_rkind, 0.0_rkind, 0.0_rkind/), (/1.0_rkind, 0.0_rkind, 1.0_rkind/), (/0.0_rkind, &
            1.0_rkind, 0.0_rkind/), ders(2,1), ders(2,2))
          
          call plane_derivative((/0.0_rkind, 0.0_rkind, 0.0_rkind/), (/1.0_rkind, 0.0_rkind, 0.0_rkind/), (/0.0_rkind, &
            1.0_rkind, 1.0_rkind/), ders(3,1), ders(3,2))
          do i=1, ubound(gauss_points%point,1)
            base_fnc(1,i) = 1 + gauss_points%point(i,1)*ders(1,1) + gauss_points%point(i,2)*ders(1,2)
          end do

          do j=2,3
            do i=1,ubound(gauss_points%point,1)
               base_fnc(j,i) = gauss_points%point(i,1)*ders(j,1) + gauss_points%point(i,2)*ders(j,2)
            end do
          end do

	      case(3)

          !to be implemented

	 
	    end select

      if (allocated(uzly)) then
        deallocate(uzly)
      end if

      if (allocated(vahy)) then
        deallocate(vahy)
      end if
      


    end subroutine init_integ


    subroutine reorder()
      use typy
      use globals
      use global_objs
      use pde_objs
      use core_tools
      use debug_tools

      integer(kind=ikind) :: i, counter, bc, j, last, proc, proc_start, el_id, nd_id

      counter = 1
      proc_start = 0
  
      do proc=1, ubound(pde,1)
        pde(proc)%procbase_fnc(1) = counter

        do i=1, nodes%kolik
          pde(proc)%permut(i) = i
          if (nodes%edge(i) /= 0) then
            el_id = nodes%element(i)%data(1)
            do j=1, ubound(elements%data,2)
              if (i == elements%data(el_id, j)) then
                nd_id = j
                EXIT
              end if
            end do
            call pde(proc)%bc(nodes%edge(i))%value_fnc(pde(proc), el_id,nd_id, code=bc)
          else
            bc = 0
          end if
          
          select case(bc)
            case(1)
              pde(proc)%permut(i) = 0
            case default
              pde(proc)%permut(i) = counter 
              counter  = counter + 1
          end select  
        end do
        proc_start = counter-1
        pde(proc)%procbase_fnc(2) = counter - 1
      end do

      

      allocate(pde_common%invpermut(counter))
      
      do proc=1, ubound(pde,1)
        do i=1, ubound(pde(proc)%permut,1)
          j = pde(proc)%permut(i)
          if (j > 0) then
            pde_common%invpermut(j) = i
          end if
        end do
      end do



    end subroutine reorder 


      !> evaluates an array that contains length of the element edge, if the element edge is boundary edge, if the element 
      !! edge is apart the boundary, zero value is supplied
      !!
      !!         3
      !!       /|  
      !!      / |                   
      !!  l3 /  |  l2
      !!    /   | 
      !!   /____|
      !!  1      2
      !!     l1
      !!
      !! this figure explains the boundary length order in the lenght array, 
      !! thus the line between the nodes 1-2 has order length(1), the line between the nodes 3-1 (or 1-3) has order lenght(3)
      !< 
      subroutine surface_integ()
        use typy
        use globals
        use globals1D
        use global_objs
        use geom_tools
        use core_tools
        use debug_tools

        integer(kind=ikind) :: i,j,k,l, found


        select case(drutes_config%dimen)
          case(1)
          
            call find_neighbours(elements, nodes)
  
            nodes%boundary(1) = .true.
            nodes%boundary(nodes%kolik) = .true.
            nodes%boundary(2:nodes%kolik-1) = .false.

            elements%length = 1.0_rkind

            elements%nvect_z = 0.0_rkind

            elements%nvect_z(1,1) = cos(angle_1D)

            elements%nvect_z(elements%kolik,2) = -cos(angle_1D)

          case(2)
          
!            if (drutes_config%check4mass) then
            
              call write_log("creating graph of the discretization mesh, and searching for boundary nodes...")
              
              call find_neighbours(elements, nodes)
              
              call set_boundary()
              
!            end if

            elements%length = 0.0_rkind

            elements%nvect_z = 0.0_rkind
            
      

            do i=1, elements%kolik
              found = -1
              do j=1,3
                if (nodes%boundary(elements%data(i,j))) then
                  do k=1,3
                    if (k/=j .and. nodes%edge(elements%data(i,j)) == nodes%edge(elements%data(i,k)) & 
                          .and. found /= elements%data(i,j) ) then

                      elements%length(i,j) = elements%length(i,j) + &
                                             dist(nodes%data(elements%data(i,j),:), nodes%data(elements%data(i,k),:))/2.0
                                             
                      elements%nvect_z(i,j) = get_nz(i,j,k)

                      elements%length(i,k) = elements%length(i,k) + &
                                              dist(nodes%data(elements%data(i,j),:), nodes%data(elements%data(i,k),:))/2.0
                      elements%nvect_z(i,k) = get_nz(i, j,k)

                      found = elements%data(i,k)

                    end if
                  end do
                end if
              end do
            end do
          case default
            print *, "not implemented, called from feminittools::surface_integ"
            ERROR STOP
        end select
            
            

	
      end subroutine surface_integ


      !> this procedure identifies boundary nodes on 2D mesh only. The method searches for all surrounding nodes and evaluates angles. If the angle is lower then 360dg, then this particular node must be the boundary node. 
      subroutine search_bnd_nodes(id, found)
        use typy
        use globals
        use global_objs
        use geom_tools
        use debug_tools

        integer(kind=ikind), intent(in) :: id
        logical, intent(out) :: found
        integer(kind=ikind) :: i,j, k, l, counter, m
        real(kind=rkind), dimension(2) :: A,B,C
        real(kind=rkind) :: ang

        ang = 0.0_rkind
        do i=1, nodes%el2integ(id)%pos
        
          j=nodes%el2integ(id)%data(i)
          do k=1,3
            if (elements%data(j,k) == id) then
              C = nodes%data(elements%data(j,k),:)
              select case(k)
                case(1)
                      A = nodes%data(elements%data(j,2),:)
                      B = nodes%data(elements%data(j,3),:)
                case(2)
                      A = nodes%data(elements%data(j,1),:)
                      B = nodes%data(elements%data(j,3),:)
                case(3)
                      A = nodes%data(elements%data(j,1),:)
                      B = nodes%data(elements%data(j,2),:)
              end select
              ang = ang + angle(A,B,C)*180.0_rkind/(4*datan(1.0d0))
            end if
          end do
          if (ang >= 360 - 100*epsilon(ang) ) then
            EXIT
          end if
        end do


        if (ang < 360 - 100*epsilon(ang)) then
          nodes%boundary(id) = .true.
          found = .true.
        else
          nodes%boundary(id) = .false.
          found = .false.
        end if


      end subroutine search_bnd_nodes
      
      subroutine set_boundary()
        use typy
        use globals
        use global_objs
        use debug_tools
        use core_tools
        use printtools
        use simplelinalg
        
        integer(kind=ikind) :: i, j, k, el, nd, elfirst, counter, nodes_at_bc, elprev, n, ii, jj, el1, el2
        logical :: found
        logical, dimension(:), allocatable, save :: nd_processed
        integer(kind=ikind), dimension(:,:), allocatable :: el_combi

        nodes%boundary = .false.
              
        if (.not. allocated(nd_processed)) then
          allocate(nd_processed(nodes%kolik))
          nd_processed = .false.
        end if
        

        do i=1, elements%kolik
          outer: do j=1, ubound(elements%data,2)
            if (elements%neighbours(i,j) == 0) then
              inner: do k = 1, ubound(elements%data,2)
                      nd = elements%data(i,k)
                      if (.not. nd_processed(nd)) then
                        call search_bnd_nodes(elements%data(i,k), found)
                        nd_processed(nd) = .true.
                      end if
              end do inner
            end if
          end do outer   
        end do	
        

        
        nodes%boundary_order=0
        
        allocate(el_combi(ubound(elements%data,2),2))
        

        do el=1, elements%kolik
          j=0
          do i=1, ubound(elements%data,2)
            nd=elements%data(el,i)
            if (nodes%boundary(nd)) j=j+1
          end do
          
          found=.false.

           if (j>1 .and. j<=ubound(elements%data,2)) then
            el_combi=0
            do i=1,ubound(elements%data,2)
              if (i<ubound(elements%data,2)) then
                el_combi(i,:) = (/elements%data(el,i), elements%data(el, i+1)/)
              else
                el_combi(i,:) =  (/elements%data(el,i), elements%data(el, 1)/)
              end if
            end do
            
            found = .true.
            elcomb: do i=1, ubound(el_combi,1)
              if (nodes%boundary(el_combi(i,1)) .and. nodes%boundary(el_combi(i,2))) then
                el_test: do ii=1, nodes%element(el_combi(i,1))%pos
                  el1 = nodes%element(el_combi(i,1))%data(ii)
                  do jj=1, nodes%element(el_combi(i,2))%pos
                    el2 = nodes%element(el_combi(i,2))%data(jj)
                    if (el1==el2 .and. el1 /= el .and. el2 /= el) then
                      found=.false.
                      exit el_test
                    end if
                  end do 
                end do el_test
                if (found) then
                  call elements%border(el)%nrfill(el_combi(i,1))
                  call elements%border(el)%nrfill(el_combi(i,2))
                end if
              end if
            end do elcomb
          
          end if
            
	 
        end do
        
        do i=1, elements%kolik
          if (elements%border(i)%pos > 0) then
            call elements%bcel%nrfill(i)
          end if
        end do


      end subroutine set_boundary
      

end module feminittools

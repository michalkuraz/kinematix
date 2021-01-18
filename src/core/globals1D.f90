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

!> \file globals1D.f90
!! \brief contains specific variables for 1D problem
!<

module globals1D
  use typy
  !> 1D mesh data
  real(kind=rkind), public ::	length_1D
  !> 1D domain discretization array
  real(kind=rkind), dimension(:,:), allocatable,  public :: deltax_1D
  !> 1D materials
  real(kind=rkind), dimension(:,:), allocatable, public :: materials_1D
  !> permutation array - how-to permutate, for 1D
  integer(kind=ikind), dimension(:), allocatable, public :: permut_vector
  integer(kind=ikind), dimension(:), allocatable, public :: permut_inverse
  !> environment aniso angle (specific for 1D -- angle of flow direction
  real(kind=rkind), public :: angle_1D

end module globals1D

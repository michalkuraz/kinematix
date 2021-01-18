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



module kinglobs
  use typy
  use global_objs
  
  type, public :: surfacend_str
    real(kind=rkind), dimension(3) :: xyz
    logical :: boundary
  end type surfacend_str
  
  type, public :: surface_el_str
    real(kind=rkind) :: sx, sy
    integer(kind=ikind) :: cover
  end type surface_el_str
  
  type, public :: inf_model_str
    character(len=7) :: name
    real(kind=rkind) :: Ks, S, A
  end type inf_model_str
  
  type, public :: solute_str
    real(kind=rkind) :: horb, csinit, rhos, lambda
  end type solute_str

  type(surfacend_str), dimension(:), allocatable, public :: watershed_nd
  type(surface_el_str), dimension(:), allocatable, public :: watershed_el
  
  real(kind=rkind), dimension(:), allocatable, public :: manning
  
  real(kind=rkind), dimension(:), allocatable, public :: oneDslopes
  
  
  type, public :: raindata_str
    real(kind=rkind), dimension(:), allocatable, public :: xy
    real(kind=rkind), dimension(:,:), allocatable, public :: series
  end type raindata_str
  
  type(raindata_str), dimension(:), allocatable, public :: raindata
  integer(kind=ikind), dimension(:), allocatable, public :: el2pt
  
  logical, public :: backwater
  
  type(inf_model_str), dimension(:), allocatable :: inf_model
  
  logical, public :: with_solutes
  
  type(solute_str), dimension(:), allocatable, public :: kinsols
  
  
  
  

end module kinglobs



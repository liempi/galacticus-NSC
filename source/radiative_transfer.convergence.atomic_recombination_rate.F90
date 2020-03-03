!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
!!    Andrew Benson <abenson@carnegiescience.edu>
!!
!! This file is part of Galacticus.
!!
!!    Galacticus is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version.
!!
!!    Galacticus is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with Galacticus.  If not, see <http://www.gnu.org/licenses/>.

  !# <radiativeTransferConvergence name="radiativeTransferConvergenceHydrogenRecombinationRate">
  !#  <description>A task which performs radiative transfer.</description>
  !# </radiativeTransferConvergence>
  type, extends(radiativeTransferConvergenceClass) :: radiativeTransferConvergenceHydrogenRecombinationRate
     !% Implementation of a radiative transfer convergence class based on the recombination rate of hydrogren.
     private
     double precision :: toleranceRelative
     double precision :: recombinationRateTotal, recombinationRateTotalPrevious
   contains
     procedure :: testConvergence => hydrogenRecombinationRateTestConvergence
  end type radiativeTransferConvergenceHydrogenRecombinationRate
  
  interface radiativeTransferConvergenceHydrogenRecombinationRate
     !% Constructors for the {\normalfont \ttfamily hydrogenRecombinationRate} radiative transfer matter class.
     module procedure hydrogenRecombinationRateConstructorParameters
     module procedure hydrogenRecombinationRateConstructorInternal
  end interface radiativeTransferConvergenceHydrogenRecombinationRate
  
contains

  function hydrogenRecombinationRateConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily hydrogenRecombinationRate} radiative transfer matter class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameters, inputParameter
    implicit none
    type            (radiativeTransferConvergenceHydrogenRecombinationRate)                :: self
    type            (inputParameters                                      ), intent(inout) :: parameters
    double precision                                                                       :: toleranceRelative
    
    !# <inputParameter>
    !#   <name>toleranceRelative</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>1.0d-3</defaultValue>
    !#   <description>The relative tolerance in total recombination rate required to declare convergence.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    self=radiativeTransferConvergenceHydrogenRecombinationRate(toleranceRelative)
    return
  end function hydrogenRecombinationRateConstructorParameters

  function hydrogenRecombinationRateConstructorInternal(toleranceRelative) result(self)
    !% Internal constructor for the {\normalfont \ttfamily hydrogenRecombinationRate} radiative transfer matter class.
    implicit none
    type            (radiativeTransferConvergenceHydrogenRecombinationRate)                :: self
    double precision                                                       , intent(in   ) :: toleranceRelative
    !# <constructorAssign variables="toleranceRelative"/>

    self%recombinationRateTotal        =-huge(0.0d0)
    self%recombinationRateTotalPrevious=-huge(0.0d0)
    return
  end function hydrogenRecombinationRateConstructorInternal
  
  subroutine hydrogenRecombinationRateTestConvergence(self,radiativeTransferMatter_,properties,statusCell,converged)
    !% Test convergence in the computational domain cell.
    use :: Galacticus_Display        , only : Galacticus_Display_Message   , verbosityStandard
    use :: MPI_Utilities             , only : mpiSelf
    use :: Radiative_Transfer_Matters, only : radiativeTransferMatterAtomic, radiativeTransferPropertiesMatterAtomic
    implicit none
    class    (radiativeTransferConvergenceHydrogenRecombinationRate), intent(inout) :: self
    class    (radiativeTransferMatterClass                         ), intent(inout) :: radiativeTransferMatter_
    class    (radiativeTransferPropertiesMatter                    ), intent(inout) :: properties
    integer                                                         , intent(in   ) :: statusCell
    logical                                                         , intent(  out) :: converged
    character(len=128                                              )                :: message
    
    ! Reset accumulated recombination rate for the first cell.
    if (statusCell == statusCellFirst) then
       self%recombinationRateTotalPrevious=self%recombinationRateTotal
       self%recombinationRateTotal        =0.0d0
    end if
    ! Accumulate recombination rate.
    select type (radiativeTransferMatter_)
    type is (radiativeTransferMatterAtomic)
       select type (properties)
       type is (radiativeTransferPropertiesMatterAtomic)
          self%recombinationRateTotal=+self                    %recombinationRateTotal                &
               &                      +radiativeTransferMatter_%recombinationRateHydrogen(properties)                   
       end select
    end select
    ! Test convergence for the last cell.
    if (statusCell == statusCellLast) then
       converged=                              abs(self%recombinationRateTotal-self%recombinationRateTotalPrevious) &
            &    <                                                                                                  &
            &     self%toleranceRelative*0.5d0*   (self%recombinationRateTotal+self%recombinationRateTotalPrevious)
       if (mpiSelf%isMaster()) then
          write (message,'(a,e12.6,a)') 'total atomic recombination rate = ',self%recombinationRateTotal,' s⁻¹ ('
          if (.not.converged) then
             message=trim(message)//'not converged)'
          else
             message=trim(message)//    'converged)'
          end if
          call Galacticus_Display_Message(trim(message),verbosityStandard)
       end if
    else
       converged=.true.
    end if
    return
  end subroutine hydrogenRecombinationRateTestConvergence
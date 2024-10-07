!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024
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

  !!{
  Implements a node operator class that performs star formation in nuclear star cluster.
  !!}

  use :: NSC_Timescales , only : NSCTimescaleClass
  !![
  <nodeOperator name="nodeOperatordarkCoreMassEvolution">
   <description>A node operator class that performs black hole formation.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatordarkCoreMassEvolution
     !!{
     A node operator class that performs the black hole evolution in dark cores.
     !!}
     private
     class(NSCTimescaleClass), pointer :: timescaleNSC_ => null()

   contains
     final     ::                                        darkCoreMassEvolutionDestructor
     procedure :: differentialEvolution               => darkCoreMassEvolutionDifferentialEvolution
  end type nodeOperatordarkCoreMassEvolution
  
  interface nodeOperatordarkCoreMassEvolution
     !!{
     Constructors for the {\normalfont \ttfamily darkCoreMassEvolutionDarkCore} node operator class.
     !!}
     module procedure darkCoreMassEvolutionConstructorParameters
     module procedure darkCoreMassEvolutionConstructorInternal
  end interface nodeOperatordarkCoreMassEvolution
  
contains

  function darkCoreMassEvolutionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily starFormation} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatordarkCoreMassEvolution)                :: self
    type (inputParameters                  ), intent(inout) :: parameters
    class(NSCTimescaleClass                ), pointer       :: timescaleNSC_
    !![
    <objectBuilder class="timescaleNSC"  name="timescaleNSC_"  source="parameters"/>
    !!]
    self=nodeOperatordarkCoreMassEvolution(timescaleNSC_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="timescaleNSC_"/>
    !!]
    return
  end function darkCoreMassEvolutionConstructorParameters

  function darkCoreMassEvolutionConstructorInternal(timescaleNSC_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily darkCoreMassEvolutionDarkCore} node operator class.
    !!}
    implicit none
    type   (nodeOperatordarkCoreMassEvolution)                        :: self
    class  (NSCTimescaleClass                ), intent(in   ), target :: timescaleNSC_
    !![
    <constructorAssign variables="*timescaleNSC_"/>
    !!]
    return
  end function darkCoreMassEvolutionConstructorInternal

  subroutine darkCoreMassEvolutionDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily darkCoreMassEvolutionDarkCore} node operator class.
    !!}
    implicit none
    type (nodeOperatordarkCoreMassEvolution), intent(inout) :: self

    !![
    <objectDestructor name="self%timescaleNSC_"/>
    !!]
    return
  end subroutine darkCoreMassEvolutionDestructor
  
  subroutine darkCoreMassEvolutionDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Perform mass of the dark core.
    !!}
    use :: Galacticus_Nodes              , only : interruptTask, nodeComponentNSC, nodeComponentDarkCore, nodeComponentDarkCoreStandard, &
          &                                       propertyInactive, treeNode
    implicit none
    class           (nodeOperatordarkCoreMassEvolution), intent(inout), target  :: self
    type            (treeNode                         ), intent(inout), target  :: node
    logical                                            , intent(inout)          :: interrupt
    procedure       (interruptTask                    ), intent(inout), pointer :: functionInterrupt
    integer                                            , intent(in   )          :: propertyType
    class           (nodeComponentNSC                 )               , pointer :: NSC
    class           (nodeComponentDarkCore            )               , pointer :: darkCore
    double precision                                                            :: dynFrictionTimescale, massDarkCore
    
    if (propertyInactive(propertyType)) return

    ! Check for a realistic dark core, return immediately if nuclear star cluster is unphysical.
    darkCore => node%darkCore()
    if (ratedarkCoreMassEvolution <= 0.0d0) return

    select type (darkCore)
    class is (nodeComponentDarkCoreStandard)
      NSC => node%NSC     ()
      !Return inmediatly if the mass of BHs in NSC is zero.
      if (NSC%massBHs() == 0.0d0) return
      dynFrictionTimescale = self%timescaleNSC_%timescale(node) 
      if (dynFrictionTimescale > 0.0d0) then 
        massDarkCoreRate   = NSC% massBHs()/dynFrictionTimescale
        call NSC     % massBHsRate    (-ratedarkCoreMassEvolution)
        call darkCore% massStellarRate(+ratedarkCoreMassEvolution)
      else
        return
      end if 
    end select
    return
  end subroutine darkCoreMassEvolutionDifferentialEvolution

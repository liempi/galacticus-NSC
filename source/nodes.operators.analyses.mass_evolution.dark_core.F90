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
  use :: Dark_Core_Growth_Rates, only: darkCoreGrowthRatesClass

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
     class (darkCoreGrowthRatesClass), pointer :: darkCoreGrowthRates_ => null()
   contains
     final     ::                          darkCoreMassEvolutionDestructor
     procedure :: differentialEvolution => darkCoreMassEvolutionDifferentialEvolution
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
    class(darkCoreGrowthRatesClass         ), pointer       :: darkCoreGrowthRates_
    !![
    <objectBuilder class="darkCoreGrowthRates" name="darkCoreGrowthRates_" source="parameters"/>
    !!]

    self=nodeOperatordarkCoreMassEvolution(darkCoreGrowthRates_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkCoreGrowthRates_"/>
    !!]
    return
  end function darkCoreMassEvolutionConstructorParameters

  function darkCoreMassEvolutionConstructorInternal(darkCoreGrowthRates_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily darkCoreMassEvolutionDarkCore} node operator class.
    !!}
    implicit none
    type (nodeOperatordarkCoreMassEvolution)                        :: self
    class(darkCoreGrowthRatesClass         ), intent(in   ), target :: darkCoreGrowthRates_
    !![
    <constructorAssign variables="*darkCoreGrowthRates_"/>
    !!]
    return
  end function darkCoreMassEvolutionConstructorInternal

  subroutine darkCoreMassEvolutionDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily NuclearStarClusterGrowth} node operator class.
    !!}
    implicit none
    type(nodeOperatordarkCoreMassEvolution), intent(inout) :: self
    
    !![
    <objectDestructor name="self%darkCoreGrowthRates_"/>
    !!]
    return
  end subroutine darkCoreMassEvolutionDestructor
  
  subroutine darkCoreMassEvolutionDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Perform mass of the dark core.
    !!}
    use :: Galacticus_Nodes    , only : interruptTask, nodeComponentNSC, propertyInactive, nodeComponentNSCStandard, &
           &                            treeNode
    use :: Histories           , only : operator(*)  , history
    use :: Abundances_Structure, only : operator(*)
    implicit none
    class           (nodeOperatordarkCoreMassEvolution), intent(inout), target  :: self
    type            (treeNode                         ), intent(inout), target  :: node
    logical                                            , intent(inout)          :: interrupt
    procedure       (interruptTask                    ), intent(inout), pointer :: functionInterrupt
    integer                                            , intent(in   )          :: propertyType
    class           (nodeComponentNSC                 )               , pointer :: nuclearStarCluster
    double precision                                                            :: rateMassDarkCore
    type            (history                          )                         :: historyTransferRate

    if (propertyInactive(propertyType)) return
    
    nuclearStarCluster => node%NSC     ()
    select type (nuclearStarCluster)
      class default
        ! Generic type, do nothing
        return
      class is (nodeComponentNSCStandard)
        rateMassDarkCore = self%darkCoreGrowthRates_%rate(node)
        if (rateMassDarkCore>0.0d0.and.nuclearStarCluster%massStellar()>0.0d0) then
          call nuclearStarCluster%massDarkCoreRate     (+rateMassDarkCore)
          call nuclearStarCluster% massStellarRate     (-rateMassDarkCore)
          call nuclearStarCluster%abundancesStellarRate(                                        &
        &                                               -rateMassDarkCore                       & 
        &                                               *nuclearStarCluster%abundancesStellar() &
        &                                               /nuclearStarCluster%massStellar()       &
        &                                              )
          historyTransferRate=nuclearStarCluster%stellarPropertiesHistory()
          if (historyTransferRate%exists()) &
              & call nuclearStarCluster%stellarPropertiesHistoryRate(-rateMassDarkCore     &
                   &                                                 *historyTransferRate  &
                   &                                                )
          call historyTransferRate%destroy()
           !! Star formation history.
           historyTransferRate=nuclearStarCluster%starFormationHistory()
          if (historyTransferRate%exists()) &
              & call nuclearStarCluster%starFormationHistoryRate    (-rateMassDarkCore    &
                   &                                                 *historyTransferRate &
             &                                                      )
          call historyTransferRate%destroy()
        end if 
      end select
      return
  end subroutine darkCoreMassEvolutionDifferentialEvolution

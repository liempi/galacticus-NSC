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
  Implements a node operator class that performs black hole formation in nuclear star cluster.
  !!}

  use :: Star_Formation_Rates_Nuclear_Star_Clusters, only : starFormationRateNuclearStarClustersClass
  !![
  <nodeOperator name="nodeOperatorBlackHoleFormationNSC">
   <description>A node operator class that performs black hole formation.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorBlackHoleFormationNSC
     !!{
     A node operator class that performs the black hole evolution in dark cores.
     !!}
     private
     class  (starFormationRateNuclearStarClustersClass), pointer :: starFormationRateNuclearStarClusters_ => null()
     double precision                                            :: efficiency

   contains
     final     ::                                        blackHoleFormationDestructor
     procedure :: differentialEvolution               => blackHoleFormationDifferentialEvolution
  end type nodeOperatorBlackHoleFormationNSC
  
  interface nodeOperatorBlackHoleFormationNSC
     !!{
     Constructors for the {\normalfont \ttfamily blackHoleFormationDarkCore} node operator class.
     !!}
     module procedure blackHoleFormationConstructorParameters
     module procedure blackHoleFormationConstructorInternal
  end interface nodeOperatorBlackHoleFormationNSC
  
contains

  function blackHoleFormationConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily starFormation} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorBlackHoleFormationNSC          )                :: self
    type (inputParameters                            ), intent(inout) :: parameters
    class  (starFormationRateNuclearStarClustersClass), pointer       :: starFormationRateNuclearStarClusters_
    double precision                                                  :: efficiency
    !![
    <inputParameter>
    <name>efficiency</name>
    <defaultValue>0.01d0</defaultValue>
    <description> Free parameter regulating the black hole formation rate</description>
    <source>parameters</source>
    </inputParameter>
    <objectBuilder class="starFormationRateNuclearStarClusters" name="starFormationRateNuclearStarClusters_" source="parameters"/>
    !!]
    self=nodeOperatorBlackHoleFormationNSC(efficiency, starFormationRateNuclearStarClusters_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="starFormationRateNuclearStarClusters_"/>
    !!]
    return
  end function blackHoleFormationConstructorParameters

  function blackHoleFormationConstructorInternal(efficiency,starFormationRateNuclearStarClusters_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily blackHoleFormationDarkCore} node operator class.
    !!}
    implicit none
    type   (nodeOperatorBlackHoleFormationNSC        )                        :: self
    class  (starFormationRateNuclearStarClustersClass), intent(in   ), target :: starFormationRateNuclearStarClusters_
    double precision                                  , intent(in   ), target :: efficiency
    !![
    <constructorAssign variables="efficiency"/>
    <constructorAssign variables="*starFormationRateNuclearStarClusters_"/>
    !!]
    return
  end function blackHoleFormationConstructorInternal

  subroutine blackHoleFormationDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily blackHoleFormationDarkCore} node operator class.
    !!}
    implicit none
    type (nodeOperatorBlackHoleFormationNSC), intent(inout) :: self

    !![
    <objectDestructor name="self%starFormationRateNuclearStarClusters_"/>
    !!]
    return
  end subroutine blackHoleFormationDestructor
  
  subroutine blackHoleFormationDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Perform black Hole formation in a nuclear star cluster.
    !!}
    use :: Galacticus_Nodes              , only : interruptTask, nodeComponentNSC, nodeComponentDarkCore, nodeComponentDarkCoreStandard, &
          &                                       propertyInactive, treeNode
    use :: Histories                     , only :  operator(*), history
    implicit none
    class           (nodeOperatorBlackHoleFormationNSC), intent(inout), target  :: self
    type            (treeNode                         ), intent(inout), target  :: node
    logical                                            , intent(inout)          :: interrupt
    procedure       (interruptTask                    ), intent(inout), pointer :: functionInterrupt
    integer                                            , intent(in   )          :: propertyType
    class           (nodeComponentNSC                 )               , pointer :: NSC
    class           (nodeComponentDarkCore            )               , pointer :: darkCore
    double precision                                                            :: rateBlackHoleFormation
    type            (history                          )                         :: historyTransferRate
    
    if (propertyInactive(propertyType)) return

    ! Check for a realistic dark core, return immediately if nuclear star cluster is unphysical.
    darkCore => node%darkCore()
    NSC      => node%NSC     ()

    rateBlackHoleFormation = self%efficiency*self%starFormationRateNuclearStarClusters_%rate(node) 

    if (rateBlackHoleFormation <= 0.0d0) return

    select type (darkCore)
    class is (nodeComponentDarkCoreStandard)
      call NSC% massStellarRate(-rateBlackHoleFormation)
      call NSC%     massBHsRate(+rateBlackHoleFormation)

      historyTransferRate=NSC%stellarPropertiesHistory()
      if (historyTransferRate%exists()) &
            & call NSC%stellarPropertiesHistoryRate (-rateBlackHoleFormation*historyTransferRate)
      call historyTransferRate%destroy()
      !! Star formation history.
      historyTransferRate=NSC%starFormationHistory()
      if (historyTransferRate%exists()) &
            & call NSC%starFormationHistoryRate (-rateBlackHoleFormation*historyTransferRate)
    end select
    return
  end subroutine blackHoleFormationDifferentialEvolution
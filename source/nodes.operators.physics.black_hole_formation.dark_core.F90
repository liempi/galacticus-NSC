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

  use :: Star_Formation_Rates_NSC , only : starFormationRateNSCClass
  !![
  <nodeOperator name="nodeOperatorBlackHoleFormation">
   <description>A node operator class that performs black hole formation.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorBlackHoleFormation
     !!{
     A node operator class that performs the black hole evolution in dark cores.
     !!}
     private
     class(starFormationRateNSCClass), pointer :: starFormationRateNSC_ => null()
     double precision                          :: efficiency

   contains
     final     ::                                        blackHoleFormationDestructor
     procedure :: differentialEvolution               => blackHoleFormationDifferentialEvolution
  end type nodeOperatorBlackHoleFormation
  
  interface nodeOperatorBlackHoleFormation
     !!{
     Constructors for the {\normalfont \ttfamily blackHoleFormationDarkCore} node operator class.
     !!}
     module procedure blackHoleFormationConstructorParameters
     module procedure blackHoleFormationConstructorInternal
  end interface nodeOperatorBlackHoleFormation
  
contains

  function blackHoleFormationConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily starFormation} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorBlackHoleFormation)                :: self
    type (inputParameters               ), intent(inout) :: parameters
    class(starFormationRateNSCClass     ), pointer       :: starFormationRateNSC_
    double precision                                     :: efficiency
    !![
    <inputParameter>
    <name>efficiency</name>
    <defaultValue>0.01d0</defaultValue>
    <description> Free parameter regulating the black hole formation rate</description>
    <source>parameters</source>
    </inputParameter>
    <objectBuilder class="starFormationRateNSC"      name="starFormationRateNSC_"  source="parameters"/>
    !!]
    self=nodeOperatorBlackHoleFormation(efficiency, starFormationRateNSC_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="starFormationRateNSC_"/>
    !!]
    return
  end function blackHoleFormationConstructorParameters

  function blackHoleFormationConstructorInternal(efficiency,starFormationRateNSC_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily blackHoleFormationDarkCore} node operator class.
    !!}
    implicit none
    type   (nodeOperatorBlackHoleFormation)                        :: self
    class  (starFormationRateNSCClass     ), intent(in   ), target :: starFormationRateNSC_
    double precision                       , intent(in   ), target :: efficiency
    !![
    <constructorAssign variables="efficiency"/>
    <constructorAssign variables="*starFormationRateNSC_"/>
    !!]
    return
  end function blackHoleFormationConstructorInternal

  subroutine blackHoleFormationDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily blackHoleFormationDarkCore} node operator class.
    !!}
    implicit none
    type (nodeOperatorBlackHoleFormation), intent(inout) :: self

    !![
    <objectDestructor name="self%starFormationRateNSC_"/>
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
    class           (nodeOperatorBlackHoleFormation), intent(inout), target  :: self
    type            (treeNode                      ), intent(inout), target  :: node
    logical                                         , intent(inout)          :: interrupt
    procedure       (interruptTask                 ), intent(inout), pointer :: functionInterrupt
    integer                                         , intent(in   )          :: propertyType
    class           (nodeComponentNSC              )               , pointer :: NSC
    class           (nodeComponentDarkCore         )               , pointer :: darkCore
    double precision                                                         :: rateBlackHoleFormation
    type            (history                       )                         :: historyTransferRate
    
    if (propertyInactive(propertyType)) return

    ! Check for a realistic dark core, return immediately if nuclear star cluster is unphysical.
    darkCore => node%darkCore()
    NSC      => node%NSC     ()

    rateBlackHoleFormation = self%efficiency*self%starFormationRateNSC_%rate(node) 

    if (rateBlackHoleFormation <= 0.0d0) return

    select type (darkCore)
    type is (nodeComponentDarkCore)
      if (rateBlackHoleFormation > 0.0d0 .and. NSC%massStellar() > 0.0d0) then
        interrupt=.true.
        functionInterrupt => DarkCoreCreate
      end if
      return
    class is (nodeComponentDarkCoreStandard)
      call NSC     % massStellarRate(-rateBlackHoleFormation)
      call darkCore% massStellarRate(+rateBlackHoleFormation)

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

    subroutine DarkCoreCreate(node,timeEnd)
  !!{
    Creates the dark ocre via interrupt.
  !!}
    use :: Galacticus_Nodes, only : interruptTask   , nodeComponentDarkCore, nodeComponentDarkCoreStandard, &
          &                         propertyInactive, treeNode
    implicit none
    type (treeNode             ), intent(inout), target  :: node
    double precision            , intent(in   ), optional:: timeEnd
    class(nodeComponentDarkCore),                pointer :: darkCore
    !$GLC attributes unused :: timeEnd
    darkCore => node%darkCore(autoCreate=.true.)
    return 
  end subroutine DarkCoreCreate
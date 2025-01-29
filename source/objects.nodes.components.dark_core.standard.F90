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
Contains a module which implements the standard dark core node component.
!!}

module Node_Component_Dark_Core_Standard
  !!{
  Implements the standard dark core node component.
  !!}
  use :: Dark_Matter_Halo_Scales         , only : darkMatterHaloScaleClass
  use :: Satellite_Merging_Mass_Movements, only : mergerMassMovementsClass
  use :: Star_Formation_Histories        , only : starFormationHistory            , starFormationHistoryClass
  use :: Stellar_Population_Properties   , only : stellarPopulationPropertiesClass
  implicit none
  private
  public :: Node_Component_Dark_Core_Standard_Scale_Set        , Node_Component_Dark_Core_Standard_Pre_Evolve          , &
       &    Node_Component_Dark_Core_Standard_Post_Step        , Node_Component_Dark_Core_Standard_Thread_Uninitialize , &
       &    Node_Component_Dark_Core_Standard_Initialize       , Node_Component_Dark_Core_Standard_Inactive            , &
       &    Node_Component_Dark_Core_Standard_State_Store      , Node_Component_Dark_Core_Standard_State_Retrieve      , &
       &    Node_Component_Dark_Core_Standard_Thread_Initialize                   

  !![
  <component>
   <class>darkCore</class>
   <name>standard</name>
   <isDefault>true</isDefault>
   <createFunction isDeferred="true" />
   <properties>
    <property>
      <name>isInitialized</name>
      <type>logical</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
    </property>
    <property>
      <name>massStellar</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
      <output unitsInSI="massSolar" comment="Mass of black holes in the standard dark core."/>
    </property>
    <property>
      <name>fractionMassRetained</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
    </property>
    <property>
      <name>radius</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
      <output unitsInSI="megaParsec" comment="Radial scale length in the standard dark core."/>
    </property>
    <property>
      <name>halfMassRadius</name>
      <attributes isSettable="false" isGettable="true" isEvolvable="false" isVirtual="true" />
      <type>double</type>
      <rank>0</rank>
      <output unitsInSI="megaParsec" comment="Radial scale length in the standard dark core."/>
      <getFunction>Node_Component_Dark_Core_Standard_Half_Mass_Radius</getFunction>
    </property>
    <property>
      <name>velocity</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
      <getFunction>Node_Component_Dark_Core_Standard_Velocity</getFunction>
      <output unitsInSI="kilo" comment="Circular velocity of the standard dark core at scale length."/>    
    </property>
    <property>
      <name>stellarPropertiesHistory</name>
      <type>history</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
    </property>
    <property>
      <name>starFormationHistory</name>
      <type>history</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
    </property>
   </properties>
    <bindings>
    <binding method="massDistribution" function="Node_Component_Dark_Core_Standard_Mass_Distribution" bindsTo="component"/>
    <binding method="massBaryonic"     function="Node_Component_Dark_Core_Standard_Mass_Baryonic"     bindsTo="component"/>
   </bindings>
     <functions>objects.nodes.components.dark_core.standard.bound_functions.inc</functions>
  </component>
  !!]

  ! Objects used by this component.
  class(darkMatterHaloScaleClass        ), pointer :: darkMatterHaloScale_
  class(stellarPopulationPropertiesClass), pointer :: stellarPopulationProperties_
  class(starFormationHistoryClass       ), pointer :: starFormationHistory_
  class(mergerMassMovementsClass        ), pointer :: mergerMassMovements_
  !$omp threadprivate(darkMatterHaloScale_,stellarPopulationProperties_,starFormationHistory_,mergerMassMovements_)

  ! Parameters controlling the physical implementation.
  double precision                            :: toleranceAbsoluteMass      , toleranceRelativeMetallicity          

  ! A threadprivate object used to track to which thread events are attached.
  integer :: thread
  !$omp threadprivate(thread)

contains

  !![
  <nodeComponentInitializationTask>
   <unitName>Node_Component_Dark_Core_Standard_Initialize</unitName>
  </nodeComponentInitializationTask>
  !!]
  subroutine Node_Component_Dark_Core_Standard_Initialize(parameters)
    !!{
    Initializes the tree node standard dark core methods module.
    !!}
    use :: Error                                 , only : Error_Report
    use :: Galacticus_Nodes                      , only : defaultDarkCoreComponent , nodeComponentDarkCoreStandard
    use :: Input_Parameters                      , only : inputParameter           , inputParameters
    implicit none
    type(inputParameters              ), intent(inout) :: parameters
    type(inputParameters              )                :: subParameters

    if (defaultDarkCoreComponent%standardIsActive()) then
       ! Find our parameters.
       subParameters=parameters%subParameters('componentDarkCore')
       
       ! Read parameters controlling the physical implementation.
       !![
       <inputParameter>
         <name>toleranceAbsoluteMass</name>
         <defaultValue>1.0d-6</defaultValue>
         <description>The mass tolerance used to judge whether the dark core is physically plausible.</description>
         <source>subParameters</source>
       </inputParameter>
       <inputParameter>
         <name>toleranceRelativeMetallicity</name>
         <defaultValue>1.0d-4</defaultValue>
         <description>The metallicity tolerance for ODE solution.</description>
         <source>subParameters</source>
       </inputParameter>
       !!]
    end if
    return
  end subroutine Node_Component_Dark_Core_Standard_Initialize

  !![
  <nodeComponentThreadInitializationTask>
   <unitName>Node_Component_Dark_Core_Standard_Thread_Initialize</unitName>
  </nodeComponentThreadInitializationTask>
  !!]
  subroutine Node_Component_Dark_Core_Standard_Thread_Initialize(parameters)
    !!{
      Initializes the standard dark core component module for each thread.
    !!}
    use :: Events_Hooks                          , only : dependencyDirectionAfter , dependencyRegEx     , openMPThreadBindingAtLevel, &
          &                                               satelliteMergerEvent
    use :: Error                                 , only : Error_Report
    use :: Galacticus_Nodes                      , only : defaultDarkCoreComponent
    use :: Input_Parameters                      , only : inputParameter           , inputParameters
    use :: Mass_Distributions                    , only : massDistributionSpherical, kinematicsDistributionLocal
    use :: Node_Component_Dark_Core_Standard_Data, only : massDistributionStellar_ , kinematicDistribution_  
    use :: Galactic_Structure_Options            , only : componentTypeDarkCore    , massTypeStellar            
                  
    implicit none
    type            (inputParameters), intent(inout) :: parameters
    type            (dependencyRegEx), dimension(1)  :: dependencies
    type            (inputParameters)                :: subParameters

    ! Check if this implementation is selected. If so, initialize the mass distribution.
    if (defaultDarkCoreComponent%standardIsActive()) then
       dependencies(1)=dependencyRegEx(dependencyDirectionAfter,'^remnantStructure:')
       call satelliteMergerEvent      %attach(thread,satelliteMerger      ,openMPThreadBindingAtLevel,label='nodeComponentDarkCoreStandard',dependencies=dependencies)
       ! Find our parameters.
       subParameters=parameters%subParameters('componentDarkCore')
       !![
       <objectBuilder class="darkMatterHaloScale"                                            name="darkMatterHaloScale_"         source="subParameters"                    />
       <objectBuilder class="stellarPopulationProperties"                                    name="stellarPopulationProperties_" source="subParameters"                    />
       <objectBuilder class="starFormationHistory"                                           name="starFormationHistory_"        source="subParameters"                    />
       <objectBuilder class="mergerMassMovements"                                            name="mergerMassMovements_"         source="subParameters"                    />
       <objectBuilder class="massDistribution"      parameterName="massDistributionDarkCore" name="massDistributionStellar_"     source="subParameters" threadPrivate="yes" >
        <default>
         <massDistributionDarkCore value="betaProfile">
          <beta value="5.0d0/3.0d0"/>
          <dimensionless value="true"/>
         </massDistributionDarkCore>
        </default>
       </objectBuilder>
       !!]
       ! Validate the Dark Core mass distribution
       select type(massDistributionStellar_)
       class is (massDistributionSpherical)
        ! The dark core distribution must have spherical symmetry.
        class default
          call Error_Report('only spherically symmetric mass distributions are allowed'//{introspection:location})
        end select
        if (.not.massDistributionStellar_%isDimensionless()) call Error_Report('dark core mass distribution must be dimensionless'//{introspection:location})
        call massDistributionStellar_%setTypes(componentTypeDarkCore,massTypeStellar)
        ! Construct the kinematic distribution
        allocate(kinematicDistribution_)
        !![
        <referenceConstruct object="kinematicDistribution_" constructor="kinematicsDistributionLocal(alpha=1.0d0/sqrt(2.0d0))"/>
        !!]
    end if
    return
  end subroutine Node_Component_Dark_Core_Standard_Thread_Initialize

  !![
  <nodeComponentThreadUninitializationTask>
   <unitName>Node_Component_Dark_Core_Standard_Thread_Uninitialize</unitName>
  </nodeComponentThreadUninitializationTask>
  !!]
  subroutine Node_Component_Dark_Core_Standard_Thread_Uninitialize()
    !!{
    Uninitializes the standard dark core component module for each thread.
    !!}
    use :: Events_Hooks                          , only : satelliteMergerEvent
    use :: Galacticus_Nodes                      , only : defaultDarkCoreComponent
    use :: Node_Component_Dark_Core_Standard_Data, only : massDistributionStellar_, kinematicDistribution_
    implicit none

    if (defaultDarkCoreComponent%standardIsActive()) then
        if (satelliteMergerEvent      %isAttached(thread,satelliteMerger      )) call satelliteMergerEvent      %detach(thread,satelliteMerger      )
       !![
       <objectDestructor name="darkMatterHaloScale_"         />
       <objectDestructor name="stellarPopulationProperties_" />
       <objectDestructor name="starFormationHistory_"        />
       <objectDestructor name="mergerMassMovements_"         />
       <objectDestructor name="massDistributionStellar_"     />
       <objectDestructor name="kinematicDistribution_"       />
      !!]
    end if
    return
  end subroutine Node_Component_Dark_Core_Standard_Thread_Uninitialize

  !![
  <preEvolveTask>
  <unitName>Node_Component_Dark_Core_Standard_Pre_Evolve</unitName>
  </preEvolveTask>
  !!]
  subroutine Node_Component_Dark_Core_Standard_Pre_Evolve(node)
    !!{
    Ensure the dark core has been initialized.
    !!}
    use :: Galacticus_Nodes, only : defaultDarkCoreComponent, nodeComponentDarkCore, nodeComponentDarkCoreStandard, treeNode
    implicit none
    type (treeNode             ), intent(inout), pointer :: node
    class(nodeComponentDarkCore)               , pointer :: darkCore

    ! Check if we are the default method.
    if (.not.defaultDarkCoreComponent%standardIsActive()) return

    ! Get the dark core component.
    darkCore => node%darkCore()

    ! Check if an standard dark core component exists.
    select type (darkCore)
    class is (nodeComponentDarkCoreStandard)
       ! Initialize the dark Core
       call Node_Component_Dark_Core_Standard_Create(node)
    end select
    return
  end subroutine Node_Component_Dark_Core_Standard_Pre_Evolve

  !![
  <postStepTask>
    <unitName>Node_Component_Dark_Core_Standard_Post_Step</unitName>
  </postStepTask>
  !!]
  subroutine Node_Component_Dark_Core_Standard_Post_Step(node,status)
    !!{
    Do processing of the node required after evolution.
    !!}
    use :: Error                         , only : Error_Report
    use :: Galacticus_Nodes              , only : defaultDarkCoreComponent, nodeComponentDarkCore , nodeComponentDarkCoreStandard, &
          &                                       treeNode
    use :: Interface_GSL                 , only : GSL_Success             , GSL_Continue
    implicit none
    type            (treeNode             ), intent(inout), pointer :: node
    integer                                , intent(inout)          :: status
    class           (nodeComponentDarkCore)               , pointer :: darkCore

    ! Return immediately if this class is not in use.
    if (.not.defaultDarkCoreComponent%standardIsActive()) return
    ! Get the dark core component.
    darkCore => node%darkCore()
    ! Check if an standard dark core component exists.
    select type (darkCore)
    class is (nodeComponentDarkCoreStandard)
       ! Note that "status" is not set to failure as these changes in state of the dark core should not change any calculation of
       ! differential evolution rates as a negative stellar mass was unphysical anyway.
       ! Trap negative stellar masses.
       if (darkCore%massStellar() < 0.0d0) then
          call darkCore%massStellarSet(0.0d0)
          call darkCore%radiusSet     (0.0d0)  
          call darkCore%velocitySet   (0.0d0)    
          ! Indicate that ODE evolution should continue after this state change.
          if (status == GSL_Success) status=GSL_Continue
       end if
    end select
    return
  end subroutine Node_Component_Dark_Core_Standard_Post_Step

  subroutine Node_Component_Dark_Core_Standard_Create(node)
    !!{
    Create properties in an standard dark core component.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkCore, treeNode
    use :: Histories       , only : history
    implicit none
    type   (treeNode             ), intent(inout), target  :: node
    class  (nodeComponentDarkCore)               , pointer :: darkCore

    ! Get the dark core component.
    darkCore => node%darkCore()

    ! Exit if already initialized.
    if (darkCore%isInitialized()) return

    call darkCore%fractionMassRetainedSet( 1.0d0)
    call darkCore%       isInitializedSet(.true.)
    return
  end subroutine Node_Component_Dark_Core_Standard_Create

  !![
  <scaleSetTask>
   <unitName>Node_Component_Dark_Core_Standard_Scale_Set</unitName>
  </scaleSetTask>
  !!]
  subroutine Node_Component_Dark_Core_Standard_Scale_Set(node)
    !!{
    Set scales for properties of {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes              , only : defaultDarkCoreComponent, nodeComponentDarkCore, nodeComponentDarkCoreStandard, nodeComponentNSC, &
          &                                       treeNode
    implicit none
    type            (treeNode              ), intent(inout), pointer :: node
    class           (nodeComponentDarkCore )               , pointer :: darkCore
    class           (nodeComponentNSC      )               , pointer :: NSC 
    double precision                        , parameter              :: massMinimum                   =1.0d0
    double precision                        , parameter              :: fractionTolerance             =1.0d-4

    ! Check if we are the default method.
    if (.not.defaultDarkCoreComponent%standardIsActive()) return
    ! Get the dark core component.
    darkCore => node%darkCore()
    ! Check if an standard dark core component exists.
    select type (darkCore)
    class is (nodeComponentDarkCoreStandard)
       NSC => node%NSC()

       call darkCore%massStellarScale      (max(1.0d-3*NSC%massStellar(), massMinimum))

       ! Set the scale for the retained stellar mass fraction.
       call darkCore%fractionMassRetainedScale(fractionTolerance*darkCore%fractionMassRetained())
    end select
    return
  end subroutine Node_Component_Dark_Core_Standard_Scale_Set

  !![
  <inactiveSetTask>
   <unitName>Node_Component_Dark_Core_Standard_Inactive</unitName>
  </inactiveSetTask>
  !!]
  subroutine Node_Component_Dark_Core_Standard_Inactive(node)
    !!{
    Set Jacobian zero status for properties of {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDarkCore, nodeComponentDarkCoreStandard, treeNode
    implicit none
    type (treeNode             ), intent(inout), pointer :: node
    class(nodeComponentDarkCore)               , pointer :: darkCore

    ! Get the dark core component.
    darkCore => node%darkCore()
    ! Check if an standard dark core component exists.
    select type (darkCore)
      class is (nodeComponentDarkCoreStandard)
      !Do nothing
    end select
    return
  end subroutine Node_Component_Dark_Core_Standard_Inactive

  subroutine satelliteMerger(self,node)
    !!{
    Transfer any standard dark core associated with {\normalfont \ttfamily node} to its host halo.
    !!}
    use :: Error                           , only : Error_Report
    use :: Galacticus_Nodes                , only : nodeComponentDarkCore    , nodeComponentDarkCoreStandard   , nodeComponentSpheroid           , &
         &                                          nodeComponentDisk        , treeNode
    use :: Satellite_Merging_Mass_Movements, only : destinationMergerSpheroid, destinationMergerDisk           , enumerationDestinationMergerType
    implicit none
    class           (*                               ), intent(inout) :: self
    type            (treeNode                        ), intent(inout) :: node
    class           (nodeComponentDarkCore           ), pointer       :: darkCore
    class           (nodeComponentSpheroid           ), pointer       :: spheroidHost           
    class           (nodeComponentDisk               ), pointer       :: diskHost               
    type            (treeNode                        ), pointer       :: nodeHost
    type            (enumerationDestinationMergerType)                :: destinationGasSatellite, destinationGasHost       , &
         &                                                               destinationStarsHost   , destinationStarsSatellite
    logical                                                           :: mergerIsMajor
    !$GLC attributes unused :: self

    ! Check that the dark core is of the standard class.
    darkCore => node%darkCore()
    select type (darkCore)
    class is (nodeComponentDarkCoreStandard)
       ! Find the node to merge with.
       nodeHost    => node%mergesWith  (                 )
       spheroidHost=> nodeHost%spheroid(autoCreate=.true.)
       diskHost    => nodeHost%disk    (autoCreate=.true.)

       ! Get mass movement descriptors.
       call mergerMassMovements_%get(node,destinationGasSatellite,destinationStarsSatellite,destinationGasHost,destinationStarsHost,mergerIsMajor)
       ! Move the gas component of the standard dark core to the host.
       select case (destinationGasSatellite%ID)
       case (destinationMergerDisk%ID)
          call diskHost%massGasSet            (diskHost    %massGas            ())
          call diskHost%abundancesGasSet      (diskHost    %abundancesGas      ())
          call diskHost%angularMomentumSet    (diskHost    %angularMomentum    ())

       case (destinationMergerSpheroid%ID)
          call spheroidHost%massGasSet        (spheroidHost%massGas            ())
          call spheroidHost%abundancesGasSet  (spheroidHost%abundancesGas      ())
          call spheroidHost%angularMomentumSet(spheroidHost%angularMomentum    ())
       case default
          call Error_Report('unrecognized movesTo descriptor'//{introspection:location})
       end select

       ! Move the stellar (black hole) property of the standard dark core to the host.  Here we assume that the Dark Core is dissolved,
       ! without adding the mass to the spheroid/disk stellar mass component.
       ! This must be fixed together with the star (black hole) formation history and stellar properties! 
       select case (destinationStarsSatellite%ID)
       case (destinationMergerDisk%ID)
          call diskHost   %massStellarSet        (diskHost %massStellar        ())
          call diskHost   %abundancesStellarSet  (diskHost %abundancesStellar  ())
          call diskHost   %luminositiesStellarSet(diskHost %luminositiesStellar())
          call diskHost   %angularMomentumSet    (diskHost %angularMomentum    ())
       case (destinationMergerSpheroid%ID)
          call spheroidHost%massStellarSet        (spheroidHost%massStellar        ())
          call spheroidHost%abundancesStellarSet  (spheroidHost%abundancesStellar  ())
          call spheroidHost%luminositiesStellarSet(spheroidHost%luminositiesStellar())
       case default
          call Error_Report('unrecognized movesTo descriptor'//{introspection:location})
       end select
       call darkCore%        massStellarSet(                  0.0d0)
       call darkCore%             radiusSet(                  0.0d0)
    end select
    return
  end subroutine satelliteMerger

  !![
  <stateStoreTask>
   <unitName>Node_Component_Dark_Core_Standard_State_Store</unitName>
  </stateStoreTask>
  !!]
  subroutine Node_Component_Dark_Core_Standard_State_Store(stateFile,gslStateFile,stateOperationID)
    !!{
    Write the tablulation state to file.
    !!}
    use            :: Display                               , only : displayMessage          , verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding                         , only : c_ptr                   , c_size_t
    use            :: Node_Component_Dark_Core_Standard_Data, only : massDistributionStellar_, kinematicDistribution_
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Storing state for: componentDarkCore -> standard',verbosity=verbosityLevelInfo)
    !![
    <stateStore variables="massDistributionStellar_ kinematicDistribution_ stellarPopulationProperties_ darkMatterHaloScale_ starFormationHistory_ mergerMassMovements_"/>
    !!]
    return
  end subroutine Node_Component_Dark_Core_Standard_State_Store

  !![
  <stateRetrieveTask>
   <unitName>Node_Component_Dark_Core_Standard_State_Retrieve</unitName>
  </stateRetrieveTask>
  !!]
  subroutine Node_Component_Dark_Core_Standard_State_Retrieve(stateFile,gslStateFile,stateOperationID)
    !!{
    Retrieve the tabulation state from the file.
    !!}
    use            :: Display                               , only : displayMessage          , verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding                         , only : c_ptr                   , c_size_t
    use            :: Node_Component_Dark_Core_Standard_Data, only : massDistributionStellar_, kinematicDistribution_
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Retrieving state for: componentDarkCore -> standard',verbosity=verbosityLevelInfo)
   !![
    <stateRestore variables="massDistributionStellar_ kinematicDistribution_ stellarPopulationProperties_ darkMatterHaloScale_ starFormationHistory_ mergerMassMovements_"/>
    !!]
    return
  end subroutine Node_Component_Dark_Core_Standard_State_Retrieve

end module Node_Component_Dark_Core_Standard

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
  use :: Galactic_Structure              , only : galacticStructureClass
  implicit none
  private
  public :: Node_Component_Dark_Core_Standard_Scale_Set        , Node_Component_Dark_Core_Standard_Pre_Evolve                  , &
       &    Node_Component_Dark_Core_Standard_Post_Step        , Node_Component_Dark_Core_Standard_Thread_Uninitialize         , &
       &    Node_Component_Dark_Core_Standard_Initialize       , Node_Component_Dark_Core_Standard_Calculation_Reset           , &
       &    Node_Component_Dark_Core_Standard_State_Store      , Node_Component_Dark_Core_Standard_State_Retrieve              , &
       &    Node_Component_Dark_Core_Standard_Thread_Initialize, Node_Component_Dark_Core_Standard_Inactive                  

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
      <name>massGas</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" makeGeneric="false" />
      <output unitsInSI="massSolar" comment="Mass of gas in the standard dark core."/>
    </property>
    <property>
      <name>abundancesGas</name>
      <type>abundances</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" makeGeneric="false" />
      <output unitsInSI="massSolar" comment="Mass of metals in the gas phase of the standard dark core."/>
    </property>
    <property>
      <name>angularMomentum</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" makeGeneric="false" />
      <output unitsInSI="massSolar*megaParsec*kilo" comment="Angular momentum of the standard dark core."/>
      <getFunction>Node_Component_Dark_Core_Standard_Angular_Momentum</getFunction>
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
    <binding method="enclosedMass"              function="Node_Component_Dark_Core_Standard_Enclosed_Mass"             bindsTo="component" />
    <binding method="acceleration"              function="Node_Component_Dark_Core_Standard_Acceleration"              bindsTo="component" />
    <binding method="tidalTensor"               function="Node_Component_Dark_Core_Standard_Tidal_Tensor"              bindsTo="component" />
    <binding method="density"                   function="Node_Component_Dark_Core_Standard_Density"                   bindsTo="component" />
    <binding method="densitySphericalAverage"   function="Node_Component_Dark_Core_Standard_Density_Spherical_Average" bindsTo="component" />
    <binding method="potential"                 function="Node_Component_Dark_Core_Standard_Potential"                 bindsTo="component" />
    <binding method="rotationCurve"             function="Node_Component_Dark_Core_Standard_Rotation_Curve"            bindsTo="component" />
    <binding method="rotationCurveGradient"     function="Node_Component_Dark_Core_Standard_Rotation_Curve_Gradient"   bindsTo="component" />
    <binding method="chandrasekharIntegral"     function="Node_Component_Dark_Core_Standard_Chandrasekhar_Integral"    bindsTo="component" />
   </bindings>
   <functions>objects.nodes.components.dark_core.standard.bound_functions.inc</functions>
  </component>
  !!]

  ! Objects used by this component.
  class(darkMatterHaloScaleClass        ), pointer :: darkMatterHaloScale_
  class(stellarPopulationPropertiesClass), pointer :: stellarPopulationProperties_
  class(starFormationHistoryClass       ), pointer :: starFormationHistory_
  class(mergerMassMovementsClass        ), pointer :: mergerMassMovements_
  class(galacticStructureClass          ), pointer :: galacticStructure_
  !$omp threadprivate(darkMatterHaloScale_,stellarPopulationProperties_,starFormationHistory_,mergerMassMovements_,galacticStructure_)

  ! Internal count of abundances.
  integer                                     :: abundancesCount

  ! Parameters controlling the physical implementation.
  double precision                            :: toleranceAbsoluteMass      , toleranceRelativeMetallicity          
  logical                                     :: DarkCoreNegativeAngularMomentumAllowed

  ! The largest and smallest angular momentum, in units of that of a circular orbit at the virial radius, considered to be physically plausible for a dark core.
  double precision, parameter                 :: angularMomentumMaximum                    =1.0d+1
  double precision, parameter                 :: angularMomentumMinimum                    =1.0d-6

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
    use :: Abundances_Structure                  , only : Abundances_Property_Count
    use :: Error                                 , only : Error_Report
    use :: Galacticus_Nodes                      , only : defaultDarkCoreComponent , nodeComponentDarkCoreStandard
    use :: Input_Parameters                      , only : inputParameter           , inputParameters
    implicit none
    type(inputParameters              ), intent(inout) :: parameters
    type(inputParameters              )                :: subParameters

    if (defaultDarkCoreComponent%standardIsActive()) then
       ! Get number of abundance properties.
       abundancesCount  =Abundances_Property_Count            ()

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
       <inputParameter>
         <name>DarkCoreNegativeAngularMomentumAllowed</name>
         <defaultValue>.true.</defaultValue>
         <description>Specifies whether or not negative angular momentum is allowed for the dark core.</description>
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
    use :: Events_Hooks                          , only : dependencyDirectionAfter         , dependencyRegEx     , openMPThreadBindingAtLevel, &
          &                                               satelliteMergerEvent
    use :: Error                                 , only : Error_Report
    use :: Galacticus_Nodes                      , only : defaultDarkCoreComponent
    use :: Input_Parameters                      , only : inputParameter                   , inputParameters
    use :: Mass_Distributions                    , only : massDistributionSymmetrySpherical
    use :: Node_Component_Dark_Core_Standard_Data, only : massDistributionDarkCore                    
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
       <objectBuilder class="galacticStructure"                                              name="galacticStructure_"           source="subParameters"                    />
       <objectBuilder class="massDistribution"      parameterName="massDistributionDarkCore" name="massDistributionDarkCore"     source="subParameters" threadPrivate="yes" >
        <default>
         <massDistributionDarkCore value="hernquist">
          <dimensionless value="true"/>
         </massDistributionDarkCore>
        </default>
       </objectBuilder>
       !!]
       if (.not.massDistributionDarkCore%isDimensionless()) call Error_Report('dark core mass distribution must be dimensionless'//{introspection:location})
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
    use :: Node_Component_Dark_Core_Standard_Data, only : massDistributionDarkCore
    implicit none

    if (defaultDarkCoreComponent%standardIsActive()) then
        if (satelliteMergerEvent      %isAttached(thread,satelliteMerger      )) call satelliteMergerEvent      %detach(thread,satelliteMerger      )
       !![
       <objectDestructor name="darkMatterHaloScale_"         />
       <objectDestructor name="stellarPopulationProperties_" />
       <objectDestructor name="starFormationHistory_"        />
       <objectDestructor name="mergerMassMovements_"         />
       <objectDestructor name="galacticStructure_"           />
       <objectDestructor name="massDistributionDarkCore"     />
       !!]
    end if
    return
  end subroutine Node_Component_Dark_Core_Standard_Thread_Uninitialize

  !![
  <calculationResetTask>
    <unitName>Node_Component_Dark_Core_Standard_Calculation_Reset</unitName>
  </calculationResetTask>
  !!]
  subroutine Node_Component_Dark_Core_Standard_Calculation_Reset(node,uniqueID)
    !!{
    Reset standard dark core structure calculations.
    !!}
    use :: Galacticus_Nodes                      , only : treeNode
    use :: Kind_Numbers                          , only : kind_int8
    use :: Node_Component_Dark_Core_Standard_Data, only : Node_Component_Dark_Core_Standard_Reset
    implicit none
    type   (treeNode ), intent(inout) :: node
    integer(kind_int8), intent(in   ) :: uniqueID
    !$GLC attributes unused :: node

    call Node_Component_Dark_Core_Standard_Reset(uniqueID)
    return
  end subroutine Node_Component_Dark_Core_Standard_Calculation_Reset

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
    Trim histories attached to the dark core.
    !!}
    use :: Abundances_Structure          , only : abs                     , zeroAbundances
    use :: Display                       , only : displayMessage          , verbosityLevelWarn
    use :: Error                         , only : Error_Report
    use :: Galacticus_Nodes              , only : defaultDarkCoreComponent, nodeComponentDarkCore , nodeComponentDarkCoreStandard, &
          &                                       nodeComponentSpin       , nodeComponentBasic    , treeNode
    use :: Interface_GSL                 , only : GSL_Success             , GSL_Continue
    use :: ISO_Varying_String            , only : assignment(=)           , operator(//)          , varying_string
    use :: String_Handling               , only : operator(//)
    implicit none
    type            (treeNode             ), intent(inout), pointer :: node
    integer                                , intent(inout)          :: status
    class           (nodeComponentDarkCore)               , pointer :: darkCore
    class           (nodeComponentBasic   )               , pointer :: basic
    class           (nodeComponentSpin    )               , pointer :: spin
    double precision                       , parameter              :: angularMomentumTolerance=1.0d-2
    double precision                       , save                   :: fractionalErrorMaximum  =0.0d+0
    double precision                                                :: massDarkCore                   , fractionalError, &
         &                                                             specificAngularMomentum
    character       (len=20               )                         :: valueString
    type            (varying_string       ), save                   :: message
    !$omp threadprivate(message)

    ! Return immediately if this class is not in use.
    if (.not.defaultDarkCoreComponent%standardIsActive()) return
    ! Get the dark core component.
    darkCore => node%darkCore()
    ! Check if an standard dark core component exists.
    select type (darkCore)
    class is (nodeComponentDarkCoreStandard)
       ! Note that "status" is not set to failure as these changes in state of the dark core should not change any calculation of
       ! differential evolution rates as a negative gas/stellar mass was unphysical anyway.
       !
       ! Trap negative gas masses.
       if (darkCore%massGas() < 0.0d0) then
          ! Check if this exceeds the maximum previously recorded error.
          fractionalError=   abs(darkCore%massGas       ()) &
               &          /(                                &
               &             abs(darkCore%massGas       ()) &
               &            +abs(darkCore%massStellar()) &
               &           )
          !$omp critical (Standard_Dark_Core_Post_Evolve_Check)
          if (fractionalError > fractionalErrorMaximum) then
             ! Report a warning.
             message='Warning: dark core has negative gas mass (fractional error exceeds any previously reported):'//char(10)
             message=message//'  Node index        = '//node%index() //char(10)
             write (valueString,'(e12.6)') darkCore%massGas       ()
             message=message//'  dark core gas mass     = '//trim(valueString)//char(10)
             write (valueString,'(e12.6)') darkCore%massStellar()
             message=message//'  dark core stellar mass = '//trim(valueString)//char(10)
             write (valueString,'(e12.6)') fractionalError
             message=message//'  Error measure     = '//trim(valueString)//char(10)
             if (fractionalErrorMaximum == 0.0d0) then
                ! This is the first time this warning has been issued, so give some extra information.
                message=message//'  Gas mass will be reset to zero (in future cases also).'//char(10)
                message=message//'  Future cases will be reported only when they exceed the previous maximum error measure.'//char(10)
                message=message//'  Negative masses are due to numerically inaccuracy in the ODE solutions.'//char(10)
                message=message//'  If significant, consider using a higher tolerance in the ODE solver.'
             end if
             call displayMessage(message,verbosityLevelWarn)
             ! Store the new maximum fractional error.
             fractionalErrorMaximum=fractionalError
          end if
          !$omp end critical (Standard_Dark_Core_Post_Evolve_Check)
          ! Get the specific angular momentum of the dark core material
          massDarkCore= darkCore%massGas       () &
               &       +darkCore%massStellar()
          if (massDarkCore == 0.0d0) then
             specificAngularMomentum=0.0d0
             call darkCore%        massStellarSet(                  0.0d0)
          else
             specificAngularMomentum=darkCore%angularMomentum()/massDarkCore
             if (specificAngularMomentum < 0.0d0) specificAngularMomentum=darkCore%radius()*darkCore%velocity()
          end if
          ! Reset the gas, abundances and angular momentum of the dark core.
          call darkCore%        massGasSet(                                    0.0d0)
          call darkCore%  abundancesGasSet(                           zeroAbundances)
          call darkCore%angularMomentumSet(specificAngularMomentum*darkCore%massStellar())
          ! Indicate that ODE evolution should continue after this state change.
          if (status == GSL_Success) status=GSL_Continue
       end if
       ! Trap negative stellar masses.
       if (darkCore%massStellar() < 0.0d0) then
          ! Check if this exceeds the maximum previously recorded error.
          fractionalError=   abs(darkCore%massStellar()) &
               &          /(                                &
               &             abs(darkCore%massGas       ()) &
               &            +abs(darkCore%massStellar()) &
               &           )
          !$omp critical (Standard_Dark_Core_Post_Evolve_Check)
          if (fractionalError > fractionalErrorMaximum) then
             ! Report a warning.
             message='Warning: dark core has negative stellar mass (fractional error exceeds any previously reported):'//char(10)
             message=message//'  Node index        = '//node%index() //char(10)
             write (valueString,'(e12.6)') darkCore%massGas       ()
             message=message//'  dark core gas mass     = '//trim(valueString)//char(10)
             write (valueString,'(e12.6)') darkCore%massStellar()
             message=message//'  dark core black holes mass = '//trim(valueString)//char(10)
             write (valueString,'(e12.6)') fractionalError
             message=message//'  Error measure     = '//trim(valueString)//char(10)
             if (fractionalErrorMaximum == 0.0d0) then
                ! This is the first time this warning has been issued, so give some extra information.
                message=message//'  Stellar mass will be reset to zero (in future cases also).'//char(10)
                message=message//'  Future cases will be reported only when they exceed the previous maximum error measure.'//char(10)
                message=message//'  Negative masses are due to numerically inaccuracy in the ODE solutions.'//char(10)
                message=message//'  If significant, consider using a higher tolerance in the ODE solver.'
             end if
             call displayMessage(message,verbosityLevelWarn)
             ! Store the new maximum fractional error.
             fractionalErrorMaximum=fractionalError
          end if
          !$omp end critical (Standard_Dark_Core_Post_Evolve_Check)
          ! Get the specific angular momentum of the dark core material
          massDarkCore= darkCore%massGas       () &
                  &    +darkCore%massStellar()
          if (massDarkCore == 0.0d0) then
             specificAngularMomentum=0.0d0
             call darkCore%      massGasSet(         0.0d0)
             call darkCore%abundancesGasSet(zeroAbundances)
          else
             specificAngularMomentum=darkCore%angularMomentum()/massDarkCore
             if (specificAngularMomentum < 0.0d0) specificAngularMomentum=darkCore%radius()*darkCore%velocity()
          end if
          ! Reset the stellar, abundances and angular momentum of the dark core.
          call darkCore%   massStellarSet(                                     0.0d0)
          call darkCore%  angularMomentumSet(specificAngularMomentum*darkCore%massGas())
          ! Indicate that ODE evolution should continue after this state change.
          if (status == GSL_Success) status=GSL_Continue
       end if
       ! Trap negative angular momentum.
       if (darkCore%angularMomentum() < 0.0d0) then
          spin  => node%spin ()
          basic => node%basic()
          if (darkCore%massStellar()+darkCore%massGas () <=0.0d0) then
             call darkCore%angularMomentumSet(0.0d0)

          else if (.not.darkCoreNegativeAngularMomentumAllowed) then
             if  (                                     &
                  &    abs(darkCore%angularMomentum()) &
                  &   /(                               &
                  &        darkCore%massStellar ()  &
                  &     +  darkCore%massGas        ()  &
                  &    )                          &
                  &   <                           &
                  &    angularMomentumTolerance   &
                  &   *spin    %angularMomentum() &
                  &   /basic   %mass           () &
                  & ) then
                call darkCore%angularMomentumSet(0.0d0)
             else
                message='negative angular momentum in dark core with positive mass'
                write (valueString,'(e12.6)') darkCore  %angularMomentum()
                message=message//char(10)//' -> angular momentum       = '//trim(valueString)
                write (valueString,'(e12.6)') darkCore  %massStellar ()
                message=message//char(10)//' -> stellar mass           = '//trim(valueString)
                write (valueString,'(e12.6)') darkCore  %massGas        ()
                message=message//char(10)//' -> gas mass               = '//trim(valueString)
                write (valueString,'(e12.6)') +spin %angularMomentum() &
                     &                        /basic%mass           ()
                message=message//char(10)//' -> angular momentum scale = '//trim(valueString)
                call Error_Report(message//{introspection:location})
             end if
          end if
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
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDarkCore, treeNode
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
    use :: Abundances_Structure          , only : abs                     , abundances           , max                          , operator(*), &
          &                                       unitAbundances
    use :: Galacticus_Nodes              , only : defaultDarkCoreComponent, nodeComponentDarkCore, nodeComponentDarkCoreStandard, treeNode
    implicit none
    type            (treeNode              ), intent(inout), pointer :: node
    class           (nodeComponentDarkCore )               , pointer :: darkCore
    double precision                        , parameter              :: massMinimum                   =1.0d0
    double precision                        , parameter              :: angularMomentumMinimum        =1.0d-1
    double precision                        , parameter              :: fractionTolerance             =1.0d-4
    type            (abundances            )                         :: abundancesScale

    ! Check if we are the default method.
    if (.not.defaultDarkCoreComponent%standardIsActive()) return
    ! Get the dark core component.
    darkCore => node%darkCore()
    ! Check if an standard dark core component exists.
    select type (darkCore)
    class is (nodeComponentDarkCoreStandard)

       call darkCore%angularMomentumScale  (angularMomentumMinimum)                                        
       call darkCore%massGasScale          (massMinimum           )
       call darkCore%massStellarScale      (massMinimum           )

       ! Set the scale for the retained stellar mass fraction.
       call darkCore%fractionMassRetainedScale(fractionTolerance*darkCore%fractionMassRetained())
       ! Set scales for abundances if necessary.
       if (abundancesCount > 0) then
          ! Set scale for abundances.
          abundancesScale= massMinimum *unitAbundances
          ! Set scale for gas abundances.
          call darkCore%abundancesGasScale    (abundancesScale)
          ! Set scale for stellar abundances.
       end if
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
    end select
    return
  end subroutine Node_Component_Dark_Core_Standard_Inactive

  subroutine satelliteMerger(self,node)
    !!{
    Transfer any standard dark core associated with {\normalfont \ttfamily node} to its host halo.
    !!}
    use :: Abundances_Structure            , only : zeroAbundances
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
       call darkCore%      massGasSet(         0.0d0)
       call darkCore%abundancesGasSet(zeroAbundances)


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
       call darkCore%    angularMomentumSet(                  0.0d0)

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
    use            :: Display                               , only : displayMessage                    , verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding                         , only : c_ptr                             , c_size_t
    use            :: Node_Component_Dark_Core_Standard_Data, only : massDistributionDarkCore
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Storing state for: componentDarkCore -> standard',verbosity=verbosityLevelInfo)
    !![
    <stateStore variables="massDistributionDarkCore darkMatterHaloScale_ mergerMassMovements_ "/>
    <workaround type="gfortran" PR="92836" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=92836">
     <description>Internal file I/O in gfortran can be non-thread safe.</description>
    </workaround>
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
    use            :: Node_Component_Dark_Core_Standard_Data, only : massDistributionDarkCore
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Retrieving state for: componentDarkCore -> standard',verbosity=verbosityLevelInfo)
    !![
    <stateRestore variables="massDistributionDarkCore darkMatterHaloScale_ mergerMassMovements_"/>
    <workaround type="gfortran" PR="92836" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=92836">
     <description>Internal file I/O in gfortran can be non-thread safe.</description>
    </workaround>
    !!]
    return
  end subroutine Node_Component_Dark_Core_Standard_State_Retrieve

end module Node_Component_Dark_Core_Standard

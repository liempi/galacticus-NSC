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
     integer          :: darkCoreRadiusID            , darkCoreGasMassID  , &
       &                 darkCoreVelocityDispersionID, darkCoreTimeScaleID
     double precision :: temperature
   contains
     procedure :: differentialEvolution      => darkCoreMassEvolutionDifferentialEvolution
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
    double precision                                        :: temperature

    !![
    <inputParameter>
      <name>temperature</name>
      <defaultValue>100d0</defaultValue>
      <description>The temperature of the gas of the nuclear star cluster.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=nodeOperatordarkCoreMassEvolution(temperature)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function darkCoreMassEvolutionConstructorParameters

  function darkCoreMassEvolutionConstructorInternal(temperature) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily darkCoreMassEvolutionDarkCore} node operator class.
    !!}
    implicit none
    type   (nodeOperatordarkCoreMassEvolution)                        :: self
    double precision                          , intent(in   )         :: temperature       
    !![
    <constructorAssign variables="temperature"/>
    <addMetaProperty component="NSC" name="darkCoreRadius"              id="self%darkCoreRadiusID"                  isEvolvable="no"  isCreator="no"  />
    <addMetaProperty component="NSC" name="darkCoreGasMass"             id="self%darkCoreGasMassID"                 isEvolvable="no"  isCreator="yes" />
    <addMetaProperty component="NSC" name="darkCoreVelocityDispersion"  id="self%darkCoreVelocityDispersionID"      isEvolvable="no"  isCreator="yes" />
    <addMetaProperty component="NSC" name="darkCoreTimeScale"           id="self%darkCoreTimeScaleID"               isEvolvable="no"  isCreator="yes" />
    !!]
    return
  end function darkCoreMassEvolutionConstructorInternal
  
  subroutine darkCoreMassEvolutionDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Perform mass of the dark core.
    !!}
    use :: Galacticus_Nodes                , only : interruptTask                  , nodeComponentNSC,  propertyInactive, nodeComponentNSCStandard, &
        &                                           treeNode
    use :: Galactic_Structure_Options      , only : componentTypeNuclearStarCluster, massTypeGaseous
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Prefixes    , only : kilo
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal , megaParsec, gigayear
    use :: Coordinates                     , only : coordinateSpherical            , assignment(=)
    use :: Ideal_Gases_Thermodynamics      , only : Ideal_Gas_Sound_Speed
    use :: Mass_Distributions              , only : massDistributionClass
    !use :: Stellar_Populations_Initial_Mass_Functions, only : initialMassFunctionChabrier2001

    implicit none
    class           (nodeOperatordarkCoreMassEvolution), intent(inout), target  :: self
    type            (treeNode                         ), intent(inout), target  :: node
    logical                                            , intent(inout)          :: interrupt
    procedure       (interruptTask                    ), intent(inout), pointer :: functionInterrupt
    integer                                            , intent(in   )          :: propertyType
    class           (nodeComponentNSC                 )               , pointer :: nuclearStarCluster
    class           (massDistributionClass            ), pointer                :: massDistributionNuclearStarCluster_ 
    !type           (initialMassFunctionChabrier2001  ),                        :: initialMassFunction
    type            (coordinateSpherical              )                         :: coordinates
    double precision                                   , parameter              :: unnormalizedNumberOfStars = 1.11d0               , nuclearStarClusterMeanMass=0.9d0     
    double precision                                                            :: gasFraction                                      , soundSpeed                      , &
      &                                                                            frictionForceStar                                , darkCoreRadius                  , &
      &                                                                            nuclearStarClusterVelocity                       , massDarkCoreRates               , &
      &                                                                            nuclearStarClusterMassStellar                    , effectiveVelocity               , &
      &                                                                            nuclearStarClusterNumberOfStars                  , nuclearStarClusterRadius        , &  
      &                                                                            nuclearStarClusterGasMassEnclosed                , nuclearStarClusterMassGas       , &
      &                                                                            dynamicalFrictionTimescaleGasNuclearStarCluster  , nuclearStarClusterDensityGas    , &
      &                                                                            dynamicalFrictionTimescaleStarsNuclearStarCluster, nuclearStarClusterTemperatureGas, &
      &                                                                            dynamicalFrictionTimescale                       , darkCoreVelocityDispersion      , &
      &                                                                            darkCoreMass
    if (propertyInactive(propertyType)) return
    
    nuclearStarCluster             => node%NSC     ()

    select type (nuclearStarCluster)
      class default
        ! Generic type, do nothing
        return
      class is (nodeComponentNSCStandard)
        nuclearStarClusterRadius       =  nuclearStarCluster%      radius() !Mpc
        nuclearStarClusterMassGas      =  nuclearStarCluster%     massGas() 
        nuclearStarClusterMassStellar  =  nuclearStarCluster% massStellar()
        darkCoreMass                   =  nuclearStarCluster%massDarkCore()
        darkCoreRadius                 =  nuclearStarCluster%floatRank0MetaPropertyGet(self%darkCoreRadiusID)
        nuclearStarClusterNumberOfStars=  nuclearStarClusterMassStellar*unnormalizedNumberOfStars

        if (nuclearStarClusterMassStellar <= 0.0d0.or.nuclearStarClusterRadius <= 0.0d0.or.nuclearStarClusterMassGas<=0.0d0) return

        gasFraction = nuclearStarClusterMassGas/nuclearStarClusterMassStellar

        effectiveVelocity          = sqrt(gravitationalConstant_internal*(nuclearStarClusterMassGas+nuclearStarClusterMassStellar)/nuclearStarClusterRadius) ! km / s
        nuclearStarClusterVelocity = sqrt(gravitationalConstant_internal*                           nuclearStarClusterMassStellar /nuclearStarClusterRadius)*(1+gasFraction) !km / s
    
        coordinates                = [nuclearStarClusterRadius,0.0d0,0.0d0]
    
        massDistributionNuclearStarCluster_ => node                               %massDistribution      (componentTypeNuclearStarCluster,massTypeGaseous)
        nuclearStarClusterDensityGas        =  massDistributionNuclearStarCluster_%density               (coordinates                                    )
        nuclearStarClusterGasMassEnclosed   =  massDistributionNuclearStarCluster_%massEnclosedBySphere  (darkCoreRadius                                 ) 
        nuclearStarClusterTemperatureGas    =  self%temperature
        darkCoreVelocityDispersion          = sqrt((gravitationalConstant_internal*darkCoreMass*(darkCoreMass*nuclearStarClusterGasMassEnclosed))/(0.5d0*nuclearStarClusterGasMassEnclosed*darkCoreRadius))

        call nuclearStarCluster%floatRank0MetaPropertySet(self%darkCoreVelocityDispersionID,darkCoreVelocityDispersion)
        call nuclearStarCluster%floatRank0MetaPropertySet(self%darkCoreGasMassID           ,nuclearStarClusterGasMassEnclosed)

        !![
          <objectDestructor name="massDistributionNuclearStarCluster_"/>
        !!]    
        soundSpeed = Ideal_Gas_Sound_Speed(nuclearStarClusterTemperatureGas)
        frictionForceStar = 0.11d0*Pi*(gravitationalConstant_internal**2)*nuclearStarClusterDensityGas*(nuclearStarClusterMeanMass**2)/(soundSpeed**2.0d0)

        dynamicalFrictionTimescaleGasNuclearStarCluster   = (0.5d0*nuclearStarClusterMassGas*effectiveVelocity)/(nuclearStarClusterNumberOfStars*frictionForceStar*nuclearStarClusterVelocity)*(megaParsec/kilo/gigayear)
        dynamicalFrictionTimescaleStarsNuclearStarCluster =  nuclearStarClusterVelocity / (frictionForceStar/nuclearStarClusterMeanMass)*(megaParsec/kilo/gigayear)
        dynamicalFrictionTimescale                        =  (dynamicalFrictionTimescaleGasNuclearStarCluster**(-1.0d0) + dynamicalFrictionTimescaleStarsNuclearStarCluster**(-1.0d0))**(-1.0d0)
        call nuclearStarCluster%floatRank0MetaPropertySet(self%darkCoreTimeScaleID,dynamicalFrictionTimescale)

    ! Trap cases where the stellar mass is to low to be considered as a nuclear star cluster.
    !Equations to be used: eq 16., 17

    ! t_dyn_friction_gas = 1/2 Mg v_eff^2 / (N F V)
    ! F∼ 0.11* 4π G^2⟨M⟩^2 ρ /c_s2

    !initialMassFunction = initialMassFunctionChabrier2001(                                                 &
    !    &                                                     massLower         =self%massLower          , &
    !    &                                                     massTransition    =self%massTransition     , &
    !    &                                                     massUpper         =self%massUpper          , &
    !    &                                                     exponent          =self%exponent           , &
    !    &                                                     massCharacteristic=self%massCharacteristic , &
    !    &                                                     sigma             =self%sigma                &
    !    &                                                )
    ! int dN/dM MdM = 1 M☉ as the mass function is normalized.
    !massInInitialMassFunction =  1.0d0 
    ! Determinates the constant to match the stellar mass of the mass function and the stellar mass of the NSC
    ! We know that M_stellar = int dN/dM MdM = C * 1M☉
    !C = massStellar/massInInitialMassFunction
    !unnormalizedNumber = initialMassFunction%numberCumulative(massLower=self%massLower, massUpper=self%massUpper)

    !meanMass = massStellar/N
    !For Chabrier the mean mass is 0.90 M☉

    !if (velocity <= 0.0d0) return
    ! Get the Crossing time in Gyr.
    !crossTimescale=+MpcPerKmPerSToGyr &
    !      &        *radius                  &
    !      &        /velocity
    ! Let's compute the relaxing time
    !relaxTimescale        = 0.138d0*(((1+q)**4)/log(N*gamma))*crossTimescale 
    !dynFrictionTimescale  = 3.330d0*(meanMass/125.0d0)*relaxTimescale

        massDarkCoreRates = nuclearStarCluster%massStellarBlackHoles()/dynamicalFrictionTimescale

        if (nuclearStarClusterMassStellar > 0.0d0) then
          call                 nuclearStarCluster%massStellarBlackHolesRate(-massDarkCoreRates)
          call                 nuclearStarCluster%massDarkCoreRate         (+massDarkCoreRates)
        end if 
      end select
      return
  end subroutine darkCoreMassEvolutionDifferentialEvolution

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
  Implements a black hole seed based on collapse of nuclear star clusters due to runaway stellar collisions.
  !!}
 
  !![
  <blackHoleSeeds name="blackHoleSeedsDarkCores">
    <description>
      A model of black hole seeds in which seeds are formed due to the collapse of nuclear star clusters into a black hole,
      based on the model of \cite{vergara_global_2023} and \cite{escala_observational_2021}.
    </description>
  </blackHoleSeeds>
  !!]

  type, extends(blackHoleSeedsClass) :: blackHoleSeedsDarkCores
     !!{
     A black hole seeds class in which seeds are formed due to the collapse of nuclear star clusters into a black hole,
     based on the model of \cite{vergara_global_2023} and \cite{escala_observational_2021}.
     !!}
     private
     double precision                                   :: massEfficiency              , velocityThreshold, &
      &                                                    massDarkCoreThreshold
     integer                                            :: darkCoreVelocityDispersionID 
    contains
     procedure :: timescale        => darkCoreTimescale
     procedure :: mass             => darkCoresMass
     procedure :: spin             => darkCoresSpin
     procedure :: formationChannel => darkCoresFormationChannel
  end type blackHoleSeedsDarkCores
  
  interface blackHoleSeedsDarkCores
     !!{
     Constructors for the {\normalfont \ttfamily DarkCores} black hole seeds class.
     !!}
     module procedure darkCoresConstructorParameters
     module procedure darkCoresConstructorInternal
  end interface blackHoleSeedsDarkCores

contains

  function darkCoresConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily DarkCores} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (blackHoleSeedsDarkCores)                :: self
    type            (inputParameters        ), intent(inout) :: parameters
    double precision                                         :: massEfficiency       , velocityThreshold, &
      &                                                         massDarkCoreThreshold
    !![
    <inputParameter>
      <name>massEfficiency</name>
      <defaultValue>1.0d-1</defaultValue>
      <description>Specifies the efficiency of the mass converted into a black hole seed.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>velocityThreshold</name>
      <defaultValue>1.0d3</defaultValue>
      <description>Specifies the velocity dispersion of the dark core to apply the seeding prescription.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massDarkCoreThreshold</name>
      <defaultValue>1.0d2</defaultValue>
      <description>Specifies the minimum mass of the dark core to apply the seeding prescription.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=blackHoleSeedsDarkCores(massEfficiency, velocityThreshold, massDarkCoreThreshold)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function darkCoresConstructorParameters
  
  function darkCoresConstructorInternal(massEfficiency,velocityThreshold,massDarkCoreThreshold) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily DarkCores} node operator class.
    !!}
    implicit none
    type            (blackHoleSeedsDarkCores)                :: self
    double precision                         , intent(in   ) :: massEfficiency
    double precision                         , intent(in   ) :: velocityThreshold
    double precision                         , intent(in   ) :: massDarkCoreThreshold
    !![
    <constructorAssign variables="massEfficiency, velocityThreshold, massDarkCoreThreshold"/>
    !!]
    !![
    <addMetaProperty component="NSC" name="darkCoreVelocityDispersion" id="self%darkCoreVelocityDispersionID" isEvolvable="no" isCreator="no" />
    !!]
    return
  end function darkCoresConstructorInternal

  double precision function darkCoreTimescale(self, node)
    !!{
      Returns the timescale associated to the seeding mechanism.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentNSC              , nodeComponentNSCStandard, treeNode
    use :: Numerical_Constants_Physical    , only : speedLight
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Numerical_Constants_Prefixes    , only : giga
    implicit none
    class           (blackHoleSeedsDarkCores), intent(inout)          :: self
    type            (treeNode               ), intent(inout)          :: node
    class           (nodeComponentNSC       )               , pointer :: nuclearStarCluster
    double precision                                                  :: nuclearStarClusterStellarMass, darkCoreVelocityDispersion, &
      &                                                                  nuclearStarClusterRadius
    double precision                         , parameter              :: massMinimumStellarMassBH=8.0d0   ! M⊙
    double precision                         , parameter              :: velocityNormalization   =1.0d1   ! km s
    double precision                         , parameter              :: timeNormalization       =5.0d-12 ! Gyr
    double precision                         , parameter              :: fractionBinaryEnergy    =1.0d0
    double precision                         , parameter              :: meanExcentricity        =1.0d0/sqrt(2.0d0)

    nuclearStarCluster           => node              %                      NSC(                                 )
    
    !Initialize the value
    darkCoreTimescale=0.0d0

    select type (nuclearStarCluster)
      class default
      ! Generic type, do nothing.
      return
      class is (nodeComponentNSCStandard)
        nuclearStarClusterRadius     =  nuclearStarCluster%                   radius(                                 ) 
        nuclearStarClusterStellarMass=  nuclearStarCluster%              massStellar(                                 )
        darkCoreVelocityDispersion   =  nuclearStarCluster%floatRank0MetaPropertyGet(self%darkCoreVelocityDispersionID)


        if (nuclearStarClusterRadius<0.0d0.or.nuclearStarClusterStellarMass<0.0d0.or.darkCoreVelocityDispersion<=0.0d0) return 

        ! We use the gravitational wave timescale as the binding energy stored in the binaries is lost via GW emission.
        ! The binaries cease to be a source of heating for the cluster and core collapse takes place. 
        ! The decay time of a BH binary with an initial separation, a, and eccentricity, e, is (Peters 1964)
        ! from https://arxiv.org/pdf/2406.13072 Eq. (3)
        darkCoreTimescale=+timeNormalization                    &
          &               *speedLight**5                        & ! The timeNormalization corrects the units
          &               *massMinimumStellarMassBH             & ! The timescale is in given in yr so we convert to Gyr.
          &               *gravitationalConstant_internal       &
          &               *fractionBinaryEnergy**(-4)           & 
          &               *(+1.0d0                              &
          &                 -meanExcentricity**2                &
          &                )**(7.0d0/2.0d0)                     &
          &               /(                                    & 
          &                 +darkCoreVelocityDispersion         & 
          &                 /velocityNormalization              &
          &                )**8                                 &
          &               /giga
    end select
    return
  end function darkCoreTimescale

  double precision function darkCoresMass(self,node) result(mass)
      !!{
        Compute the nuclear star cluster collapse condition.
      !!}
    use :: Galacticus_Nodes                , only : nodeComponentNSC               , nodeComponentNSCStandard, treeNode  
    use :: Galactic_Structure_Options      , only : componentTypenuclearStarCluster, massTypeStellar
    implicit none
    class           (blackHoleSeedsDarkCores), intent(inout)          :: self
    type            (treeNode               ), intent(inout)          :: node
    class           (nodeComponentNSC       )               , pointer :: nuclearStarCluster
    double precision                                                  :: velocityDispersionDarkCore
    
    ! Get the nuclear star cluster component.
    nuclearStarCluster => node%NSC()
    ! Assume zero mass.
    mass=+0.0d0
    ! Detect the type of the nuclear star cluster component.
    select type (nuclearStarCluster)
      class default
          ! Generic type, do nothing.
          return
      class is (nodeComponentNSCStandard)
          ! Standard class, get the properties of the nuclear star cluster component.
          ! Unphysical nuclear star cluster, do nothing.
          if   (                                           &
             &   nuclearStarCluster%massStellar() <= 0.0d0 &
             &  .or.                                       &
             &   nuclearStarCluster%radius     () <= 0.0d0 &
             & ) return

          ! Get the velocity dispersion of the dark core.
          velocityDispersionDarkCore = nuclearStarCluster%floatRank0MetaPropertyGet(self%darkCoreVelocityDispersionID)

          ! Generic type - interrupt and create a black hole if velocity dispersion is greater than the threshold.
          if ( self%velocityThreshold <= velocityDispersionDarkCore .and. nuclearStarCluster%massDarkCore()>= self%massDarkCoreThreshold) then
            mass= max(                                                        &
              &       +self%massEfficiency*nuclearStarCluster%massDarkCore(), & 
              &       +8.0d0                                                  & 
              &      )
            ! Adjust black hole stellar mass of the nuclear star cluster
            call nuclearStarCluster%massDarkCoreSet(                                   &
              &                                     +nuclearStarCluster%massDarkCore() &
              &                                     -mass                              &
              &                                    )
          else
            mass=+0.0d0
          end if 
    end select
    return
  end function darkCoresMass

  double precision function darkCoresSpin(self,node) result(spin)
    !!{
    Compute the spin of the seed black hole.
    !!}
    implicit none
    class(blackHoleSeedsDarkCores), intent(inout) :: self
    type (treeNode               ), intent(inout) :: node
    !$GLC attributes unused :: self, node

    ! Assume zero spin.
    spin=0.0d0
    return
  end function darkCoresSpin

  function darkCoresFormationChannel (self,node) result(channel)
    !!{
    Compute the spin of the seed black hole.
    !!}
    implicit none
    type (enumerationBlackHoleFormationChannelType)                :: channel
    class(blackHoleSeedsDarkCores                 ), intent(inout) :: self
    type (treeNode                                ), intent(inout) :: node
    !$GLC attributes unused :: self, node

    channel=blackHoleFormationChannelDarkCoreCollapse
    return
  end function darkCoresFormationChannel

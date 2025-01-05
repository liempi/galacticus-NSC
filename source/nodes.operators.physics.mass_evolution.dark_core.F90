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

  use :: Mass_Distributions, only : massDistributionClass, kinematicsDistributionClass

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
     double precision                       :: efficiency         , massLower         , &
       &                                       massTransition     , massUpper         , &
       &                                       exponent           , massCharacteristic, &
       &                                       sigma
   contains
     final     ::                               darkCoreMassEvolutionDestructor
     procedure :: differentialEvolution      => darkCoreMassEvolutionDifferentialEvolution
  end type nodeOperatordarkCoreMassEvolution
  
  interface nodeOperatordarkCoreMassEvolution
     !!{
     Constructors for the {\normalfont \ttfamily darkCoreMassEvolutionDarkCore} node operator class.
     !!}
     module procedure darkCoreMassEvolutionConstructorParameters
     module procedure darkCoreMassEvolutionConstructorInternal
  end interface nodeOperatordarkCoreMassEvolution

  class (massDistributionClass ), pointer :: massDistribution_, massDistributionStellarNSC_ 
  !$omp threadprivate(massDistribution_,massDistributionStellarNSC_)

  
contains

  function darkCoreMassEvolutionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily starFormation} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatordarkCoreMassEvolution)                :: self
    type (inputParameters                  ), intent(inout) :: parameters
    double precision                                        :: efficiency         , massLower         , &
       &                                                       massTransition     , massUpper         , &
       &                                                       exponent           , massCharacteristic, &
       &                                                       sigma  
    !![
    <inputParameter>
      <name>efficiency</name>
      <defaultValue>0.01d0</defaultValue>
      <description>The efficiency of star formation for the Crossing time method.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massLower</name>
      <description>The lower mass limit for the Chabrier 2001 IMF.</description>
      <defaultValue>0.10d0</defaultValue>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massTransition</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The transition limit for the Chabrier 2001 IMF.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
    <name>massUpper</name>
    <description>The upper mass limit for the Chabrier 2001 IMF.</description>
    <defaultValue>1.25d2</defaultValue>
    <source>parameters</source>
  </inputParameter>
  <inputParameter>
    <name>exponent</name>
    <defaultValue>-2.30d0</defaultValue>
    <description>The exponent apperaing in Chabrier initial mass function .</description>
    <source>parameters</source>
  </inputParameter>
  <inputParameter>
    <name>massCharacteristic</name>
    <defaultValue>0.08d0</defaultValue>
    <description>} Characteristic mass of the lognormal part of the Chabrier 2001 IMF.</description>
    <source>parameters</source>
  </inputParameter>
  <inputParameter>
    <name>sigma</name>
    <defaultValue>0.69d0</defaultValue>
    <description>The exponent of the power law part of the Chabrier 2001 IMF.</description>
    <source>parameters</source>
  </inputParameter>
    !!]
    self=nodeOperatordarkCoreMassEvolution(efficiency,massLower,massTransition,massUpper,exponent,massCharacteristic,sigma)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function darkCoreMassEvolutionConstructorParameters

  function darkCoreMassEvolutionConstructorInternal(efficiency,massLower,massTransition,massUpper,exponent,massCharacteristic,sigma) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily darkCoreMassEvolutionDarkCore} node operator class.
    !!}
    implicit none
    type   (nodeOperatordarkCoreMassEvolution)                        :: self
    double precision                          , intent(in   )         :: efficiency        , massLower         , &
          &                                                              massTransition    , massUpper         , &
          &                                                              exponent          , massCharacteristic, &
          &                                                              sigma  
    !![
    <constructorAssign variables="efficiency,massLower,massTransition,massUpper,exponent,massCharacteristic,sigma"/>
    !!]
    return
  end function darkCoreMassEvolutionConstructorInternal

  subroutine darkCoreMassEvolutionDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily darkCoreMassEvolutionDarkCore} node operator class.
    !!}
    implicit none
    type (nodeOperatordarkCoreMassEvolution), intent(inout) :: self
    return
  end subroutine darkCoreMassEvolutionDestructor
  
  subroutine darkCoreMassEvolutionDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Perform mass of the dark core.
    !!}
    use :: Galacticus_Nodes                          , only : interruptTask                  , nodeComponentNSC, nodeComponentDarkCore, nodeComponentDarkCoreStandard, &
          &                                                   propertyInactive               , treeNode
    use :: Galactic_Structure_Options                , only : componentTypeNSC               , massTypeStellar
    use :: Numerical_Constants_Astronomical          , only : Mpc_per_km_per_s_To_Gyr
    use :: Stellar_Populations_Initial_Mass_Functions, only : initialMassFunctionChabrier2001

    implicit none
    class           (nodeOperatordarkCoreMassEvolution), intent(inout), target  :: self
    type            (treeNode                         ), intent(inout), target  :: node
    logical                                            , intent(inout)          :: interrupt
    procedure       (interruptTask                    ), intent(inout), pointer :: functionInterrupt
    integer                                            , intent(in   )          :: propertyType
    class           (nodeComponentNSC                 )               , pointer :: NSC
    class           (nodeComponentDarkCore            )               , pointer :: darkCore
    type            (initialMassFunctionChabrier2001  ),                        :: initialMassFunction
    double precision                                                            :: velocity                 , radius              , &
      &                                                                            massStellar              , massGas             , &
      &                                                                            q                        , N_un                , &
      &                                                                            N                        , C                   , &
      &                                                                            gamma                    , meanMass            , &
      &                                                                            crossTimescale           , relaxTimescale      , &
      &                                                                            massInInitialMassFunction, dynFrictionTimescale, &
      &                                                                            massDarkCoreRate 
    
    if (propertyInactive(propertyType)) return

    ! Check for a realistic dark core, return immediately if nuclear star cluster is unphysical.
    darkCore => node%darkCore()

    select type (darkCore)
    class is (nodeComponentDarkCoreStandard)
      NSC        => node%NSC     ()
      radius     =  self%efficiency*NSC%     radius() !Mpc
      massGas    =                  NSC%    massGas() 
      massStellar=                  NSC%massStellar()
      gamma      = 0.4
    
      !Return inmediatly if the mass of BHs in NSC is zero.
      if (NSC%massBHs() == 0.0d0) return

      ! Trap cases where there is no stellar component and return 0.0.
      if (massStellar > 0.0d0) then
        q = massGas/massStellar
      else
        return
      end if

      massDistributionStellarNSC_ => node%massDistribution(componentType=componentTypeNSC, massType=massTypeStellar)
      velocity = massDistributionStellarNSC_%rotationCurve(radius*1.0e-6)*(1+q)
    
      initialMassFunction =initialMassFunctionChabrier2001(                                                 &
          &                                                     massLower         =self%massLower         , &
          &                                                     massTransition    =self%massTransition    , &
          &                                                     massUpper         =self%massUpper         , &
          &                                                     exponent          =self%exponent          , &
          &                                                     massCharacteristic=self%massCharacteristic, &
          &                                                     sigma             =self%sigma               &
          &                                                )
      massInInitialMassFunction =  1
      ! Determinates the constant to match the stellar mass of the mass function and the stellar
      ! mass of the NSC
      C = massStellar/massInInitialMassFunction

      N_un = initialMassFunction%numberCumulative(                                   &
          &                                       massLower         =self%massLower, &
          &                                       massUpper         =self%massUpper   )
      N        = C*N_un
      meanMass = massStellar/N_un
    
      if (velocity <= 0.0d0) then
        dynFrictionTimescale =0.0d0
      else if (self%efficiency == 0.0d0) then
        dynFrictionTimescale =0.0d0
      else
        ! Get the Crossing time in Gyr.
        crossTimescale=+Mpc_per_km_per_s_To_Gyr &
              &        *radius                  &
              &        /velocity
      ! Let's compute the relaxing time
       relaxTimescale        = 0.138*(((1+q)**4)/(log(N*gamma)))*crossTimescale 
       dynFrictionTimescale  = 0.333*(meanMass/self%massLower)  *relaxTimescale
       end if

      if (dynFrictionTimescale > 0.0d0) then 
        massDarkCoreRate   = NSC% massBHs()/dynFrictionTimescale
        call NSC     % massBHsRate    (-massDarkCoreRate)
        call darkCore% massStellarRate(+massDarkCoreRate)
      end if 
    end select
    return
  end subroutine darkCoreMassEvolutionDifferentialEvolution

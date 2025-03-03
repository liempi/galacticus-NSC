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
     double precision                       :: efficiency
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
    double precision                                        :: efficiency

    !![
    <inputParameter>
      <name>efficiency</name>
      <defaultValue>0.01d0</defaultValue>
      <description>The efficiency of star formation for the Crossing time method.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=nodeOperatordarkCoreMassEvolution(efficiency)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function darkCoreMassEvolutionConstructorParameters

  function darkCoreMassEvolutionConstructorInternal(efficiency) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily darkCoreMassEvolutionDarkCore} node operator class.
    !!}
    implicit none
    type   (nodeOperatordarkCoreMassEvolution)                        :: self
    double precision                          , intent(in   )         :: efficiency    
    !![
    <constructorAssign variables="efficiency"/>
    !!]
    return
  end function darkCoreMassEvolutionConstructorInternal
  
  subroutine darkCoreMassEvolutionDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Perform mass of the dark core.
    !!}
    use :: Galacticus_Nodes                          , only : interruptTask                  , nodeComponentDarkCore  , nodeComponentNSC, nodeComponentDarkCoreStandard, &
          &                                                   propertyInactive               , treeNode
    use :: Galactic_Structure_Options                , only : componentTypeNuclearStarCluster, massTypeStellar
    use :: Numerical_Constants_Astronomical          , only : gravitationalConstant_internal , MpcPerKmPerSToGyr
    !use :: Stellar_Populations_Initial_Mass_Functions, only : initialMassFunctionChabrier2001

    implicit none
    class           (nodeOperatordarkCoreMassEvolution), intent(inout), target  :: self
    type            (treeNode                         ), intent(inout), target  :: node
    logical                                            , intent(inout)          :: interrupt
    procedure       (interruptTask                    ), intent(inout), pointer :: functionInterrupt
    integer                                            , intent(in   )          :: propertyType
    class           (nodeComponentNSC                 )               , pointer :: NSC
    class           (nodeComponentDarkCore            )               , pointer :: darkCore
    !type            (initialMassFunctionChabrier2001  ),                        :: initialMassFunction
    double precision                                                            :: velocity                 , radius              , &
      &                                                                            massStellar              , massGas             , &
      &                                                                            q                        , unnormalizedNumber  , &
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
      
      ! Trap cases where the stellar mass is to low to be considered as a nuclear star cluster.
      if (massStellar < 1.0d3 .or. radius <= 0.0d0) return

      q = massGas/massStellar
  
      !massDistributionStellarNSC_ => node%massDistribution(componentType=componentTypeNSC, massType=massTypeStellar)
      !velocity = massDistributionStellarNSC_%rotationCurve(1.0e6*radius)*(1+q)
      velocity = sqrt(gravitationalConstant_internal*massStellar/radius)*(1+q)

      !initialMassFunction = initialMassFunctionChabrier2001(                                                 &
      !    &                                                     massLower         =self%massLower          , &
      !    &                                                     massTransition    =self%massTransition     , &
      !    &                                                     massUpper         =self%massUpper          , &
      !    &                                                     exponent          =self%exponent           , &
      !    &                                                     massCharacteristic=self%massCharacteristic , &
      !    &                                                     sigma             =self%sigma                &
      !    &                                                )
      ! int dN/dM MdM = 1 M☉ as the mass function is normalized.
      massInInitialMassFunction =  1.0d0 
      ! Determinates the constant to match the stellar mass of the mass function and the stellar mass of the NSC
      ! We know that M_stellar = int dN/dM MdM = C * 1M☉
      !C = massStellar/massInInitialMassFunction
      !unnormalizedNumber = initialMassFunction%numberCumulative(massLower=self%massLower, massUpper=self%massUpper)
      unnormalizedNumber = 1.11d0     
      N                  = massStellar*unnormalizedNumber
      !meanMass = massStellar/N
      !For Chabrier the mean mass is 0.90 M☉
      meanMass = 0.90d0

      if (velocity <= 0.0d0) return

      ! Get the Crossing time in Gyr.
      crossTimescale=+MpcPerKmPerSToGyr &
            &        *radius                  &
            &        /velocity
      ! Let's compute the relaxing time
      relaxTimescale        = 0.138d0*(((1+q)**4)/log(N*gamma))*crossTimescale 
      dynFrictionTimescale  = 3.330d0*(meanMass/125.0d0)*relaxTimescale

      massDarkCoreRate = NSC% massBHs()/dynFrictionTimescale
      call NSC     % massBHsRate    (-massDarkCoreRate)
      call darkCore% massStellarRate(+massDarkCoreRate) 
    end select
    return
  end subroutine darkCoreMassEvolutionDifferentialEvolution

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
  Implementation of a timescale for star formation which scales with the component Crossing time.
  !!}

  use :: Galactic_Structure                        , only : galacticStructureClass


  !![
  <NSCTimescale name="NSCTimescaleDynamicalFrictionTime">
   <description>
    A timescale class in which the crossing timescale is computed as
    time. Specifically:
    \begin{equation}
     t_{\rm cross} = \epsilon_r \left( {R \over V(\epsilon_r)},
    \end{equation}
    where $\epsilon_r=${\normalfont \ttfamily [efficiency]}  and $r$ and $V$
    are the characteristic radius and velocity respectively of the component.
   </description>
  </NSCTimescale>
  !!]
  type, extends(NSCTimescaleClass) :: NSCTimescaleDynamicalFrictionTime
     !!{
     Implementation of a timescale for star formation which scales with the Crossing time.
     !!}
     private
     class  (galacticStructureClass         ), pointer :: galacticStructure_         => null()
     double precision                                  :: efficiency
     double precision                                  :: massLower
     double precision                                  :: massTransition
     double precision                                  :: massUpper
     double precision                                  :: exponent
     double precision                                  :: massCharacteristic
     double precision                                  :: sigma
   contains
  
     final     ::              DynamicalFrictionTimeTimescaleDestructor 
     procedure :: timescale => DynamicalFrictionTimeTimescale
  end type NSCTimescaleDynamicalFrictionTime

  interface NSCTimescaleDynamicalFrictionTime
     !!{
     Constructors for the {\normalfont \ttfamily DynamicalFrictionTime} timescale for star formation class.
     !!}
     module procedure DynamicalFrictionTimeConstructorParameters
     module procedure DynamicalFrictionTimeConstructorInternal
  end interface NSCTimescaleDynamicalFrictionTime

contains

  function DynamicalFrictionTimeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily DynamicalFrictionTime} timescale for star formation class which takes a parameter set as
    input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (NSCTimescaleDynamicalFrictionTime)                :: self
    type            (inputParameters                  ), intent(inout) :: parameters
    class           (galacticStructureClass           ), pointer       :: galacticStructure_
    double precision                                                   :: efficiency
    double precision                                                   :: massLower
    double precision                                                   :: massTransition
    double precision                                                   :: massUpper
    double precision                                                   :: exponent
    double precision                                                   :: massCharacteristic
    double precision                                                   :: sigma  

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
    <objectBuilder class="galacticStructure"         name="galacticStructure_"         source="parameters"/>
    !!]
    self=NSCTimescaleDynamicalFrictionTime(efficiency,massLower,massTransition,massUpper,exponent,massCharacteristic,sigma,galacticStructure_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="galacticStructure_" />
    !!]
    return
  end function DynamicalFrictionTimeConstructorParameters

  function DynamicalFrictionTimeConstructorInternal(efficiency,massLower,massTransition,massUpper,exponent,massCharacteristic,sigma,galacticStructure_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily DynamicalFrictionTime} timescale for star formation class.
    !!}
    implicit none
    type            (NSCTimescaleDynamicalFrictionTime)                        :: self
    class           (galacticStructureClass           ), intent(in   ), target :: galacticStructure_
    double precision                                   , intent(in   )         :: efficiency
    double precision                                   , intent(in   )         :: massLower
    double precision                                   , intent(in   )         :: massTransition
    double precision                                   , intent(in   )         :: massUpper
    double precision                                   , intent(in   )         :: exponent
    double precision                                   , intent(in   )         :: massCharacteristic
    double precision                                   , intent(in   )         :: sigma  
 
    
    !![
    <constructorAssign variables="efficiency,massLower,massTransition,massUpper,exponent,massCharacteristic,sigma,*galacticStructure_"/>
    !!]

    return
  end function DynamicalFrictionTimeConstructorInternal

  subroutine DynamicalFrictionTimeTimescaleDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily NSCTimescaleDynamicalFrictionTime}  class
    !!} 
    implicit none
    type(NSCTimescaleDynamicalFrictionTime), intent(inout) :: self
    !![
    <objectDestructor name="self%galacticStructure_"/>
    !!]
    return
  end subroutine DynamicalFrictionTimeTimescaleDestructor
  
  double precision function DynamicalFrictionTimeTimescale(self,node)
    !!{
    Returns the crossing timescale (in Gyr) for star formation in the given {\normalfont \ttfamily component}. The timescale is given by

    !!}
    use :: Galacticus_Nodes                          , only : nodeComponentNSC                , treeNode
    use :: Galactic_Structure_Options                , only : componentTypeNSC                , massTypeStellar
    use :: Stellar_Populations_Initial_Mass_Functions, only : initialMassFunctionChabrier2001
    use :: Numerical_Integration2                    , only : integratorCompositeTrapezoidal1D
    implicit none
    class           (NSCTimescaleDynamicalFrictionTime), intent(inout) :: self
    type            (treeNode                         ), intent(inout) :: treeNode
    type            (initialMassFunctionChabrier2001  ),               :: initialMassFunction
    type            (integratorCompositeTrapezoidal1D )                :: integrator_
    class           (nodeComponentNSC                 ), pointer       :: NSC
  
    double precision                                                   :: velocity             , radius               , &
      &                                                                   massStellar          , massGas              , &
      &                                                                   q                    , N_un                 , &
      &                                                                   gamma                , meanMass             , &
      &                                                                   CrossingTimeTimescale, RelaxingTimeTimescale         
    ! Check for zero velocity.
    NSC        => node%           NSC            ()
    radius     =  self%efficiency*NSC%     radius() !Mpc
    massGas    =                  NSC%    gasMass() 
    massStellar=                  NSC%stellarMass()
    gamma      = 0.4
    
    ! Trap cases where there is no stellar component and return 0.0.
    if (massStellar > 0.0d0) then
      q = massGas/massStellar
    else
      DynamicalFrictionTimeTimescale = +0.0d0
      return
    end if 
    
    velocity =  self%galacticStructure_%velocityRotation   (                                    &
         &                                                     node                           , &
         &                                                     radiusNSC                      , &
         &                                                     componentType=componentTypeNSC , &
         &                                                     massType     =massTypeStellar    &
         &                                                  )                                   &
         &     *(1+q)

    initialMassFunction =initialMassFunctionChabrier2001  (                                                                                                                                                &
         &                                                     massLower         =self%massLower                                                                                                    , &
         &                                                     massTransition    =self%massTransition                                                                                                    , &
         &                                                     massUpper         =self%massUpper                                                                                                    , &
         &                                                     exponent          =self%exponent                                                                                                    , &
         &                                                     massCharacteristic=self%massCharacteristic                                                                                                    , &
         &                                                     sigma             =self%sigma                                                                                                      &
         &                                                )
    ! first we need to call the integrator
    call integrator_%initialize  (24           )
    call integrator_%toleranceSet(1.0d-7,1.0d-7)
    call integrator_%integrandSet(initialMassFunctionIntegrand)

    massInInitialMassFunction =  integrator_%evaluate(                                     &
         &                                              initialMassFunction%massMinimum(), &
         &                                              initialMassFunction%massMaximum()  &
         &                                           )
    ! Determinates the constant to match the stellar mass of the mass function and the stellar
    ! mass of the NSC
    C = massStellar/massInInitialMassFunction

    call integrator_%integrandSet(initialMassFunctionIntegrandNumber)

    N_un = integrator_%evaluate(                                    &
         &                       initialMassFunction%massMinimum(), &
         &                       initialMassFunction%massMaximum()  &
         &                     )
    N        = C*N_un
    meanMass = massInInitialMassFunction/N_un
    
    if (velocity <= 0.0d0) then
       DynamicalFrictionTimeTimescale=0.0d0
    else if (self%efficiency == 0.0d0)
       DynamicalFrictionTimeTimescale=0.0d0
    else
       ! Get the Crossing time in Gyr.
       CrossingTimeTimescale=+Mpc_per_km_per_s_To_Gyr &
            &        *radius                  &
            &        /velocity
       ! Let's compute the relaxing time

      RelaxingTimeTimescale           = 0.138*(((1+q)**4)/(log(N*gamma)))*DynamicalFrictionTimeTimescale 
      DynamicalFrictionTimeTimescale  = 0.33*(meanMass/self%massLower)*RelaxingTimeTimescale
    end if
    return
  end function DynamicalFrictionTimeTimescale


  double precision function initialMassFunctionIntegrand(mass)
    !!{
    Integrand used to find the total mass in the initial mass function.
    !!}
    implicit none
    double precision, intent(in   ) :: mass

    initialMassFunctionIntegrand=+       mass  &
         &                       *initialMassFunction%phi(mass)
    return
  end function initialMassFunctionIntegrand

  double precision function initialMassFunctionIntegrandNumber(mass)
  !!{
  Integrand used to find the total number of stars in the initial mass function.
  !!}
  implicit none
  double precision, intent(in   ) :: mass

  initialMassFunctionIntegrand= initialMassFunction%phi(mass)
  return
end function initialMassFunctionIntegrandNumber
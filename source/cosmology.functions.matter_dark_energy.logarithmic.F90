!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
  An implementation of the cosmological functions class for cosmologies consisting of collisionless
  matter and dark energy with an equation of state of the form: $P=\rho^w$ with $w(a)=w_0+w_1 a (1-a)$.
  !!}

  use :: Root_Finder, only : rootFinder

  integer         , parameter :: ageTableNPointsPerDecade     =300
  double precision, parameter :: ageTableNPointsPerOctave     =dble(ageTableNPointsPerDecade)*log(2.0d0)/log(10.0d0)
  double precision, parameter :: ageTableIncrementFactor      =exp(int(ageTableNPointsPerOctave+1.0d0)*log(10.0d0)/dble(ageTableNPointsPerDecade))
  integer         , parameter :: distanceTableNPointsPerDecade=100

  ! Factor by which one component of Universe must dominate others such that we can ignore the others.
  double precision, parameter :: factorDominate               =100.0d0

  ! Variables used in root finding.
  double precision            :: factorDominateCurrent
  !$omp threadprivate(factorDominateCurrent)

  !![
  <cosmologyFunctions name="cosmologyFunctionsMatterDarkEnergyLogarithmic">
   <description>
    Cosmological relations are computed assuming a universe that contains only collisionless matter and dark energy with a logarithmic
    equation of state $w(a)=w_0+w_1\ln{a}$, with $w_0=${\normalfont \ttfamily
    [darkEnergyEquationOfStateW0]}, and $w_1=${\normalfont \ttfamily [darkEnergyEquationOfStateW1]}.
   </description>
  </cosmologyFunctions>
  !!]
  type, extends(cosmologyFunctionsMatterLambda) :: cosmologyFunctionsMatterDarkEnergyLogarithmic
     !!{
     A cosmological functions class for cosmologies consisting of matter plus dark energy with equation of state $w(a)=w_0+w_1\ln{a}$.
     !!}
     private
     double precision             :: darkEnergyEquationOfStateW0, darkEnergyEquationOfStateW1
     type            (rootFinder) :: finderDomination           , finderEquality
   contains
     !![
     <methods>
       <method description="Set a module-scope pointer to {\normalfont \ttfamily self}."                         method="targetSelf"                  />
       <method description="Return the derivative of the dark energy exponent with respect to expansion factor." method="exponentDarkEnergyDerivative"/>
     </methods>
     !!]
     procedure :: cosmicTime                    => matterDarkEnergyLogarithmicCosmicTime
     procedure :: omegaDarkEnergyEpochal        => matterDarkEnergyLogarithmicOmegaDarkEnergyEpochal
     procedure :: hubbleParameterEpochal        => matterDarkEnergyLogarithmicHubbleParameterEpochal
     procedure :: hubbleParameterRateOfChange   => matterDarkEnergyLogarithmicHubbleParameterRateOfChange
     procedure :: equationOfStateDarkEnergy     => matterDarkEnergyLogarithmicEquationOfStateDarkEnergy
     procedure :: exponentDarkEnergy            => matterDarkEnergyLogarithmicExponentDarkEnergy
     procedure :: exponentDarkEnergyDerivative  => matterDarkEnergyLogarithmicExponentDarkEnergyDerivative
     procedure :: equalityEpochMatterDarkEnergy => matterDarkEnergyLogarithmicEqualityEpochMatterDarkEnergy
     procedure :: dominationEpochMatter         => matterDarkEnergyLogarithmicDominationEpochMatter
     procedure :: distanceComoving              => matterDarkEnergyLogarithmicDistanceComoving
     procedure :: timeAtDistanceComoving        => matterDarkEnergyLogarithmicTimeAtDistanceComoving
     procedure :: distanceComovingConvert       => matterDarkEnergyLogarithmicDistanceComovingConvert
     procedure :: expansionFactorTabulate       => matterDarkEnergyLogarithmicMakeExpansionFactorTable
     procedure :: targetSelf                    => matterDarkEnergyLogarithmicTargetSelf
  end type cosmologyFunctionsMatterDarkEnergyLogarithmic

  ! Module scope pointer to the current object.
  class(cosmologyFunctionsMatterDarkEnergyLogarithmic), pointer :: self_ => null()
  !$omp threadprivate(self_)

  interface cosmologyFunctionsMatterDarkEnergyLogarithmic
     !!{
     Constructors for the matter plus dark energy cosmological functions class.
     !!}
     module procedure matterDarkEnergyLogarithmicConstructorParameters
     module procedure matterDarkEnergyLogarithmicConstructorInternal
  end interface cosmologyFunctionsMatterDarkEnergyLogarithmic

contains

  function matterDarkEnergyLogarithmicConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the matter plus dark energy cosmological functions class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (cosmologyFunctionsMatterDarkEnergyLogarithmic)                :: self
    type            (inputParameters                              ), intent(inout) :: parameters
    class           (cosmologyParametersClass                     ), pointer       :: cosmologyParameters_
    double precision                                                               :: darkEnergyEquationOfStateW0, darkEnergyEquationOfStateW1

    !![
    <inputParameter>
      <name>darkEnergyEquationOfStateW0</name>
      <source>parameters</source>
      <defaultValue>-1.0d0</defaultValue>
      <description>The equation of state parameter for dark energy, $w_0$, defined such that $P=\rho^w$ with $w(a)=w_0+w_1\ln{a}$.</description>
    </inputParameter>
    <inputParameter>
      <name>darkEnergyEquationOfStateW1</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>The equation of state parameter for dark energy, $w_1$, defined such that $P=\rho^w$ with $w(a)=w_0+w_1\ln{a}$.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !!]
    ! Use it to construct a matter plus dark energy cosmological functions class.
    self=cosmologyFunctionsMatterDarkEnergyLogarithmic(                             &
         &                                             cosmologyParameters_       , &
         &                                             darkEnergyEquationOfStateW0, &
         &                                             darkEnergyEquationOfStateW1  &
         &                                            )
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    !!]
    return
  end function matterDarkEnergyLogarithmicConstructorParameters

  function matterDarkEnergyLogarithmicConstructorInternal(cosmologyParameters_,darkEnergyEquationOfStateW0,darkEnergyEquationOfStateW1) result(self)
    !!{
    Constructor for the matter plus dark energy cosmological functions class.
    !!}
    use :: Cosmology_Parameters, only : cosmologyParametersClass
    use :: Root_Finder         , only : rangeExpandMultiplicative
    implicit none
    type            (cosmologyFunctionsMatterDarkEnergyLogarithmic)                        :: self
    class           (cosmologyParametersClass                     ), intent(in   ), target :: cosmologyParameters_
    double precision                                               , intent(in   )         :: darkEnergyEquationOfStateW0, darkEnergyEquationOfStateW1
    double precision                                                                       :: rangeExpandDownward        , rangeExpandUpward          , &
         &                                                                                    darkEnergyExponentCurrent
    !![
    <constructorAssign variables="*cosmologyParameters_, darkEnergyEquationOfStateW0, darkEnergyEquationOfStateW1"/>
    !!]
    
    self%collapsingUniverse                  =.false.
    self%enableRangeChecks                   =.true.
    self%expansionFactorMaximum              =0.0d0
    self%timeTurnaround                      =0.0d0
    self%timeMaximum                         =0.0d0
    self%expansionRatePrevious               =-1.0d0
    self%expansionRateExpansionFactorPrevious=-1.0d0
    ! Build root finder for epoch of matter domination.
    if (self%cosmologyParameters_%OmegaDarkEnergy() /= 0.0d0)                                        &
         & self%finderDomination=rootFinder(                                                         &
         &                                  rootFunction     =matterDarkEnergyLogarithmicDomination, &
         &                                  toleranceAbsolute=0.0d+0                               , &
         &                                  toleranceRelative=1.0d-6                                 &
         &                                 )
    ! Build root finder for epoch of matter equality.
    darkEnergyExponentCurrent=self%exponentDarkEnergy(expansionFactor=1.0d0)
    if (darkEnergyExponentCurrent > -3.0d0) then
       ! Dark energy density is decaying less rapidly than matter.
       if (self%cosmologyParameters_%OmegaMatter() < self%cosmologyParameters_%OmegaDarkEnergy()) then
          ! Matter density is less than dark energy density. Search backward for epoch of domination.
          rangeExpandUpward  =1.0d0
          rangeExpandDownward=0.5d0
       else
          ! Matter density is greater than dark energy density. Search forward for epoch of domination.
          rangeExpandUpward  =2.0d0
          rangeExpandDownward=1.0d0
       end if
    else
       ! Dark energy density is decaying more rapidly than matter.
       if (self%cosmologyParameters_%OmegaMatter() < self%cosmologyParameters_%OmegaDarkEnergy()) then
          ! Matter density is less than dark energy density. Search forward for epoch of domination.
          rangeExpandUpward  =2.0d0
          rangeExpandDownward=1.0d0
       else
          ! Matter density is greater than dark energy density. Search backward for epoch of domination.
          rangeExpandUpward  =1.0d0
          rangeExpandDownward=0.50d0
       end if
    end if
    self%finderEquality=rootFinder(                                                           &
         &                         rootFunction       =matterDarkEnergyLogarithmicDomination, &
         &                         toleranceAbsolute  =0.0d0                                , &
         &                         toleranceRelative  =1.0d-6                               , &
         &                         rangeExpandUpward  =rangeExpandUpward                    , &
         &                         rangeExpandDownward=rangeExpandDownward                  , &
         &                         rangeExpandType    =rangeExpandMultiplicative              &
         &                        )
    ! Force a build of the expansion factor table, which will determine if this Universe collapses.
    call self%expansionFactorTabulate()
   return
  end function matterDarkEnergyLogarithmicConstructorInternal

  double precision function matterDarkEnergyLogarithmicCosmicTime(self,expansionFactor,collapsingPhase)
    !!{
    Return the cosmological matter density in units of the critical density at the present day.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (cosmologyFunctionsMatterDarkEnergyLogarithmic), intent(inout)           :: self
    double precision                                               , intent(in   )           :: expansionFactor
    logical                                                        , intent(in   ), optional :: collapsingPhase
    logical                                                                                  :: collapsingPhaseActual

    ! Validate the input.
    call self%epochValidate(                                         &
         &                  expansionFactorIn=expansionFactor      , &
         &                  collapsingIn     =collapsingPhase      , &
         &                  collapsingOut    =collapsingPhaseActual  &
         &                 )
    ! Ensure tabulation is initialized.
    if (.not.self%ageTableInitialized) call self%expansionFactorTabulate(self%ageTableTimeMinimum)
    ! Ensure that the tabulation spans a sufficient range of expansion factors.
    if (collapsingPhaseActual) then
       ! We assume that the universe does not collapse.
       call Error_Report('non-collapsing universe assumed'//{introspection:location})
    else
       ! In expanding phase ensure that sufficiently small and large expansion factors have been reached.
       do while (self%ageTableExpansionFactor(                        1) > expansionFactor)
          self%ageTableTimeMinimum=    self%ageTableTimeMinimum/ageTableIncrementFactor
          call self%expansionFactorTabulate()
       end do
       do while (self%ageTableExpansionFactor(self%ageTableNumberPoints) < expansionFactor)
          self%ageTableTimeMaximum=max(self%ageTableTimeMaximum*ageTableIncrementFactor,self%timeTurnaround)
          call self%expansionFactorTabulate()
       end do
    end if
    ! Interpolate to get cosmic time.
    matterDarkEnergyLogarithmicCosmicTime=self%interpolatorTime%interpolate(expansionFactor)
    return
  end function matterDarkEnergyLogarithmicCosmicTime

  double precision function matterDarkEnergyLogarithmicOmegaDarkEnergyEpochal(self,time,expansionFactor,collapsingPhase)
    !!{
    Return the dark energy density parameter at expansion factor {\normalfont \ttfamily expansionFactor}.
    !!}
    use :: Cosmology_Parameters, only : hubbleUnitsStandard
    use :: Error               , only : Error_Report
    implicit none
    class           (cosmologyFunctionsMatterDarkEnergyLogarithmic), intent(inout)           :: self
    double precision                                               , intent(in   ), optional :: expansionFactor      , time
    logical                                                        , intent(in   ), optional :: collapsingPhase
    double precision                                                                         :: expansionFactorActual
    !$GLC attributes unused :: collapsingPhase

    call self%epochValidate(                                          &
         &                  timeIn            =time                 , &
         &                  expansionFactorIn =expansionFactor      , &
         &                  collapsingIn      =collapsingPhase      , &
         &                  expansionFactorOut=expansionFactorActual  &
         &                 )
    matterDarkEnergyLogarithmicOmegaDarkEnergyEpochal                                                                                      &
         & =                                   self%cosmologyParameters_%OmegaDarkEnergy       (                                         ) &
         &             *expansionFactorActual**self                     %exponentDarkEnergy    (expansionFactor    =expansionFactorActual) &
         &             *(                                                                                                                  &
         &                                     self%cosmologyParameters_%HubbleConstant        (hubbleUnitsStandard                      ) &
         &               /                     self                     %HubbleParameterEpochal(expansionFactor    =expansionFactorActual) &
         &              )**2
    return
  end function matterDarkEnergyLogarithmicOmegaDarkEnergyEpochal

  double precision function matterDarkEnergyLogarithmicDominationEpochMatter(self,dominateFactor)
    use :: Cosmology_Functions_Parameters, only : requestTypeExpansionFactor
    use :: Root_Finder                   , only : rangeExpandMultiplicative
    implicit none
    class           (cosmologyFunctionsMatterDarkEnergyLogarithmic), intent(inout) :: self
    double precision                                               , intent(in   ) :: dominateFactor
    double precision                                                               :: aDominantCurvature , aDominantDarkEnergy      , &
         &                                                                            aMatterEquality    , darkEnergyExponentCurrent, &
         &                                                                            rangeExpandDownward, rangeExpandUpward

    ! Choose present day as default - will be used if no other densities present (i.e. Einstein-de Sitter).
    matterDarkEnergyLogarithmicDominationEpochMatter=1.0d0
    ! Case where dark energy is present.
    if (self%cosmologyParameters_%OmegaDarkEnergy() /= 0.0d0) then
       darkEnergyExponentCurrent=self%exponentDarkEnergy(expansionFactor=1.0d0)
       if (darkEnergyExponentCurrent > -3.0d0) then
          ! Dark energy density is decaying less rapidly than matter.
          if (self%cosmologyParameters_%OmegaMatter() < dominateFactor*self%cosmologyParameters_%OmegaDarkEnergy()) then
             ! Matter density is less than dominated dark energy density. Search backward for epoch of domination.
             rangeExpandUpward  =1.0d0
             rangeExpandDownward=0.5d0
          else
             ! Matter density is greater than dominated dark energy density. Search forward for epoch of domination.
             rangeExpandUpward  =2.0d0
             rangeExpandDownward=1.0d0
          end if
       else
          ! Dark energy density is decaying more rapidly than matter.
          if (self%cosmologyParameters_%OmegaMatter() < dominateFactor*self%cosmologyParameters_%OmegaDarkEnergy()) then
             ! Matter density is less than dominated dark energy density. Search forward for epoch of domination.
             rangeExpandUpward  =2.0d0
             rangeExpandDownward=1.0d0
          else
             ! Matter density is greater than dominated dark energy density. Search backward for epoch of domination.
             rangeExpandUpward  =1.0d0
             rangeExpandDownward=0.50d0
          end if
       end if
       call self%finderDomination%rangeExpand(                                               &
            &                                 rangeExpandUpward  =rangeExpandUpward        , &
            &                                 rangeExpandDownward=rangeExpandDownward      , &
            &                                 rangeExpandType    =rangeExpandMultiplicative  &
            &                                )
       factorDominateCurrent = dominateFactor       
       call self%targetSelf()
       aDominantDarkEnergy  =self%finderDomination%find(rootGuess=1.0d0)
       ! Choose earliest expansion factor.
       matterDarkEnergyLogarithmicDominationEpochMatter=min(matterDarkEnergyLogarithmicDominationEpochMatter,aDominantDarkEnergy)
    end if
    if (self%cosmologyParameters_%OmegaCurvature() /= 0.0d0) then
       ! Find the expansion factor of matter-curvature equality.
       aMatterEquality=self%equalityEpochMatterCurvature(requestTypeExpansionFactor)
       ! Find the earlier expansion factor at which matter dominates by the specified amount (ratio of matter
       ! to curvature density scales as the expansion factor).
       aDominantCurvature=aMatterEquality/dominateFactor
       ! Choose earliest expansion factor.
       matterDarkEnergyLogarithmicDominationEpochMatter=min(matterDarkEnergyLogarithmicDominationEpochMatter,aDominantCurvature)
    end if
    return
  end function matterDarkEnergyLogarithmicDominationEpochMatter

  double precision function matterDarkEnergyLogarithmicHubbleParameterEpochal(self,time,expansionFactor,collapsingPhase)
    !!{
    Returns the Hubble parameter at the request cosmological time, {\normalfont \ttfamily time}, or expansion factor, {\normalfont \ttfamily expansionFactor}.
    !!}
    use :: Cosmology_Parameters, only : hubbleUnitsStandard
    use :: Error               , only : Error_Report
    implicit none
    class           (cosmologyFunctionsMatterDarkEnergyLogarithmic), intent(inout)           :: self
    double precision                                               , intent(in   ), optional :: expansionFactor      , time
    logical                                                        , intent(in   ), optional :: collapsingPhase
    double precision                                                                         :: expansionFactorActual, sqrtArgument
    
    call self%epochValidate(                                          &
         &                  timeIn            =time                 , &
         &                  expansionFactorIn =expansionFactor      , &
         &                  collapsingIn      =collapsingPhase      , &
         &                  expansionFactorOut=expansionFactorActual  &
         &                 )
    ! Compute the Hubble parameter at the specified expansion factor.
    sqrtArgument=                                                                                            &
         &       max(                                                                                        &
         &            self%cosmologyParameters_%OmegaMatter    ()                                            &
         &           /expansionFactorActual**3                                                               &
         &           +self%cosmologyParameters_%OmegaDarkEnergy()                                            &
         &           *expansionFactorActual**self%exponentDarkEnergy(expansionFactor=expansionFactorActual)  &
         &           +self%cosmologyParameters_%OmegaCurvature ()                                            &
         &           /expansionFactorActual**2                                                             , &
         &           0.0d0                                                                                   &
         &          )
    matterDarkEnergyLogarithmicHubbleParameterEpochal=self%cosmologyParameters_%HubbleConstant(hubbleUnitsStandard)*sqrt(sqrtArgument)
    ! Make the Hubble parameter negative if we are in the collapsing phase of the Universe.
    if (self%collapsingUniverse) then
       if    (present(time           )) then
          if    (time>self%timeTurnaround) matterDarkEnergyLogarithmicHubbleParameterEpochal=-matterDarkEnergyLogarithmicHubbleParameterEpochal
       else
          if (present(collapsingPhase)) then
             if (collapsingPhase         ) matterDarkEnergyLogarithmicHubbleParameterEpochal=-matterDarkEnergyLogarithmicHubbleParameterEpochal
          end if
       end if
    end if
    return
  end function matterDarkEnergyLogarithmicHubbleParameterEpochal

  double precision function matterDarkEnergyLogarithmicHubbleParameterRateOfChange(self,time,expansionFactor,collapsingPhase)
    !!{
    Returns the rate of change of the Hubble parameter at the requested cosmological time, {\normalfont \ttfamily time}, or expansion factor, {\normalfont \ttfamily expansionFactor}.
    !!}
    implicit none
    class           (cosmologyFunctionsMatterDarkEnergyLogarithmic), intent(inout)           :: self
    double precision                                               , intent(in   ), optional :: expansionFactor      , time
    logical                                                        , intent(in   ), optional :: collapsingPhase
    double precision                                                                         :: expansionFactorActual

    call self%epochValidate(                                          &
         &                  timeIn            =time                 , &
         &                  expansionFactorIn =expansionFactor      , &
         &                  collapsingIn      =collapsingPhase      , &
         &                  expansionFactorOut=expansionFactorActual  &
         &                 )
    ! Compute the rate of change of the Hubble parameter.
    matterDarkEnergyLogarithmicHubbleParameterRateOfChange                                                     &
         & =0.5d0                                                                                              &
         & *self%hubbleParameterEpochal(expansionFactor=expansionFactorActual,collapsingPhase=collapsingPhase) &
         & *self%expansionRate         (                expansionFactorActual                                ) &
         & /(                                                                                                  &
         &   +self%cosmologyParameters_%OmegaMatter    ()                                                      &
         &   /expansionFactorActual**3                                                                         &
         &   +self%cosmologyParameters_%OmegaDarkEnergy()                                                      &
         &   *expansionFactorActual**self%exponentDarkEnergy(expansionFactor=expansionFactorActual)            &
         &   +self%cosmologyParameters_%OmegaCurvature ()                                                      &
         &   /expansionFactorActual**2                                                                         &
         &  )                                                                                                  &
         & *(                                                                                                  &
         &   -3.0d0*self%cosmologyParameters_%OmegaMatter()                                                    &
         &   /expansionFactorActual**3                                                                         &
         &   +self%cosmologyParameters_%OmegaDarkEnergy  ()                                                    &
         &   *expansionFactorActual**self%exponentDarkEnergy(expansionFactor=expansionFactorActual)            &
         &   *(                                                                                                &
         &     +self%exponentDarkEnergy(expansionFactor=expansionFactorActual)                                 &
         &     +expansionFactorActual                                                                          &
         &     *log(expansionFactorActual)                                                                     &
         &     *self%exponentDarkEnergyDerivative(expansionFactor=expansionFactorActual)                       &
         &    )                                                                                                &
         &   -2.0d0*self%cosmologyParameters_%OmegaCurvature()                                                 &
         &   /expansionFactorActual**2                                                                         &
         & )
    return
  end function matterDarkEnergyLogarithmicHubbleParameterRateOfChange

  double precision function matterDarkEnergyLogarithmicEqualityEpochMatterDarkEnergy(self,requestType)
    !!{
    Return the epoch of matter-dark energy magnitude equality (either expansion factor or cosmic time).
    !!}
    use :: Cosmology_Functions_Parameters, only : requestTypeExpansionFactor, requestTypeTime
    implicit none
    class  (cosmologyFunctionsMatterDarkEnergyLogarithmic), intent(inout)           :: self
    integer                                               , intent(in   ), optional :: requestType
    integer                                                                         :: requestTypeActual

    if (present(requestType)) then
       requestTypeActual=requestType
    else
       requestTypeActual=requestTypeExpansionFactor
    end if
    factorDominateCurrent =  1.0d0
    call self%targetSelf()
    matterDarkEnergyLogarithmicEqualityEpochMatterDarkEnergy=self%finderEquality%find(rootGuess=1.0d0)
    if (present(requestType)) then
       if (requestType == requestTypeTime)                                               &
            &                  matterDarkEnergyLogarithmicEqualityEpochMatterDarkEnergy  &
            & =self%cosmicTime(matterDarkEnergyLogarithmicEqualityEpochMatterDarkEnergy)
    end if
    return
  end function matterDarkEnergyLogarithmicEqualityEpochMatterDarkEnergy

  double precision function matterDarkEnergyLogarithmicDomination(expansionFactor)
    !!{
    Function used in root finding when seeking epoch at which matter dominates over dark energy.
    !!}
    implicit none
    double precision, intent(in   ) :: expansionFactor

    matterDarkEnergyLogarithmicDomination=                                                                            &
         &                                 self_%cosmologyParameters_%OmegaMatter    ()                               &
         &                                /expansionFactor**3                                                         &
         &                                -factorDominateCurrent                                                      &
         &                                *self_%cosmologyParameters_%OmegaDarkEnergy()                               &
         &                                *expansionFactor**self_%exponentDarkEnergy(expansionFactor=expansionFactor)
    return
  end function matterDarkEnergyLogarithmicDomination

  subroutine matterDarkEnergyLogarithmicTargetSelf(self)
    !!{
    Set a module-scope pointer to the current dark energy cosmology functions object.
    !!}
    implicit none
    class(cosmologyFunctionsMatterDarkEnergyLogarithmic), intent(in   ), target :: self

    self_ => self
    return
  end subroutine matterDarkEnergyLogarithmicTargetSelf

  subroutine matterDarkEnergyLogarithmicMakeExpansionFactorTable(self,time)
    !!{
    Builds a table of expansion factor vs. time for dark energy universes.
    !!}
    use :: Cosmology_Parameters, only : hubbleUnitsTime
    use :: Numerical_Ranges    , only : Make_Range     , rangeTypeLogarithmic
    implicit none
    class           (cosmologyFunctionsMatterDarkEnergyLogarithmic)             , intent(inout), target   :: self
    double precision                                                            , intent(in   ), optional :: time
    double precision                                               , parameter                            :: turnaroundTimeTolerance         =1.0d-12
    double precision                                               , parameter                            :: odeToleranceAbsolute            =1.0d-9 , odeToleranceRelative    =1.0d-9
    double precision                                               , allocatable, dimension(:)            :: ageTableExpansionFactorTemporary        , ageTableTimeTemporary
    integer                                                                                               :: iTime                                   , prefixPointCount
    double precision                                                                                      :: OmegaDominant                           , deltaTime                      , &
         &                                                                                                   densityPower                            , timeDominant                   , &
         &                                                                                                   expansionFactorDominant                 , timeActual
    logical                                                                                               :: solutionFound                           , timeExceeded

    ! Find expansion factor early enough that a single component dominates the evolution of the Universe.
    call self%densityScalingEarlyTime(factorDominate,densityPower,expansionFactorDominant,OmegaDominant)
    ! Find the corresponding time.
    timeDominant=-2.0d0/densityPower/self%cosmologyParameters_%HubbleConstant(hubbleUnitsTime)/sqrt(OmegaDominant)/expansionFactorDominant**(0.5d0*densityPower)
    ! Find minimum and maximum times to tabulate.
    if (present(time)) then
       timeActual=time
       do while (self%ageTableTimeMinimum > min(timeActual,timeDominant)/2.0d0)
          self%ageTableTimeMinimum=self%ageTableTimeMinimum/ageTableIncrementFactor
       end do
       do while (self%ageTableTimeMaximum < max(timeActual,timeDominant)*2.0d0)
          self%ageTableTimeMaximum=self%ageTableTimeMaximum*ageTableIncrementFactor
       end do
    else
       do while (self%ageTableTimeMinimum > timeDominant/2.0d0)
          self%ageTableTimeMinimum=self%ageTableTimeMinimum/ageTableIncrementFactor
       end do
       do while (self%ageTableTimeMaximum < timeDominant*2.0d0)
          self%ageTableTimeMaximum=self%ageTableTimeMaximum*ageTableIncrementFactor
       end do
    end if
    ! Determine number of points to tabulate.
    self%ageTableNumberPoints=int(log10(self%ageTableTimeMaximum/self%ageTableTimeMinimum)     *dble(ageTableNPointsPerDecade))+1
    self%ageTableTimeMaximum =self%ageTableTimeMinimum*10.0d0**(dble(self%ageTableNumberPoints)/dble(ageTableNPointsPerDecade))
    ! Assume this Universe does not collapse initially.
    self%collapsingUniverse    =.false.
    self%expansionFactorMaximum=0.0d0
    self%timeTurnaround        =0.0d0
    self%timeMaximum           =0.0d0
    ! Deallocate arrays if currently allocated.
    if (allocated(self%ageTableTime)) then
       ! Determine number of points that are being added at the start of the array.
       prefixPointCount=int(log10(self%ageTableTime(1)/self%ageTableTimeMinimum)*dble(ageTableNPointsPerDecade)+0.5d0)
       call Move_Alloc(self%ageTableTime           ,ageTableTimeTemporary           )
       call Move_Alloc(self%ageTableExpansionFactor,ageTableExpansionFactorTemporary)
       ! Allocate the arrays to current required size.
       allocate(self%ageTableTime           (self%ageTableNumberPoints))
       allocate(self%ageTableExpansionFactor(self%ageTableNumberPoints))
       ! Create set of grid points in time variable.
       self%ageTableTime=Make_Range(self%ageTableTimeMinimum,self%ageTableTimeMaximum,self%ageTableNumberPoints,rangeTypeLogarithmic)
       ! Set the expansion factors to a negative value to indicate they are not yet computed.
       self%ageTableExpansionFactor=-1.0d0
       ! Paste in the previously computed regions.
       self%ageTableTime           (prefixPointCount+1:prefixPointCount+size(ageTableTimeTemporary))=ageTableTimeTemporary
       self%ageTableExpansionFactor(prefixPointCount+1:prefixPointCount+size(ageTableTimeTemporary))=ageTableExpansionFactorTemporary
       ! Deallocate the temporary arrays.
       deallocate(ageTableTimeTemporary           )
       deallocate(ageTableExpansionFactorTemporary)
    else
       ! Allocate the arrays to current required size.
       allocate(self%ageTableTime           (self%ageTableNumberPoints))
       allocate(self%ageTableExpansionFactor(self%ageTableNumberPoints))
       ! Create set of grid points in time variable.
       self%ageTableTime=Make_Range(self%ageTableTimeMinimum,self%ageTableTimeMaximum,self%ageTableNumberPoints,rangeTypeLogarithmic)
       ! Set the expansion factors to a negative value to indicate they are not yet computed.
       self%ageTableExpansionFactor=-1.0d0
    end if
    ! Compute quantities required for table interpolation.
    self%ageTableTimeLogarithmicMinimum=log(self%ageTableTimeMinimum)
    self%ageTableInverseDeltaLogTime   =dble(self%ageTableNumberPoints-1)/log(self%ageTableTimeMaximum/self%ageTableTimeMinimum)
    ! For the initial time, we approximate that we are at sufficiently early times that a single component dominates the
    ! Universe and use the appropriate analytic solution.
    if (self%ageTableExpansionFactor(1) < 0.0d0)                        &
         &    self%ageTableExpansionFactor           (               1) &
         & =(                                                           &
         &   -0.5d0                                                     &
         &   *densityPower                                              &
         &   *self%ageTableTime                       (              1) &
         &   *self%cosmologyParameters_%HubbleConstant(hubbleUnitsTime) &
         &   *sqrt(OmegaDominant)                                       &
         &  )**(-2.0d0/densityPower)
    ! Solve ODE to get corresponding expansion factors.
    self%iTableTurnaround  =  self%ageTableNumberPoints
    call self%targetSelf()
    do iTime=2,self%ageTableNumberPoints
       ! Compute the expansion factor if it is not already computed.
       if (self%ageTableExpansionFactor(iTime) < 0.0d0) then
          self%ageTableExpansionFactor(iTime)=matterDarkEnergyLogarithmicExpansionFactorChange(                                       &
               &                                                                    self%ageTableTime           (iTime-1), &
               &                                                                    self%ageTableTime           (iTime  ), &
               &                                                                    self%ageTableExpansionFactor(iTime-1)  &
               &                                                                   )
          ! Check for a universe which is no longer expanding (i.e. has reached its maximum expansion).
          if (self%ageTableExpansionFactor(iTime) == self%ageTableExpansionFactor(iTime-1)) then
             ! Record that we have a collapsing Universe.
             self%collapsingUniverse=.true.
             ! Record the maximum expansion factor.
             self%expansionFactorMaximum=self%ageTableExpansionFactor(iTime-1)
             ! Find the time of maximum expansion by bisection. Disable checks of epoch ranges while we search for the maximum
             ! expansion factor as we may exceed the maximum while searching.
             self%timeTurnaround   =(self%ageTableTime(iTime-1)+self%ageTableTime(iTime-2))/2.0d0
             deltaTime             =(self%ageTableTime(iTime-1)-self%ageTableTime(iTime-2))/2.0d0
             solutionFound         =.false.
             self%enableRangeChecks=.false.
             do while (.not.solutionFound)
                timeExceeded=(                                                                &
                     &         matterDarkEnergyLogarithmicExpansionFactorChange(                         &
                     &                                 self%ageTableTime           (iTime-2), &
                     &                                 self%timeTurnaround                  , &
                     &                                 self%ageTableExpansionFactor(iTime-2)  &
                     &                                )                                       &
                     &        >=                                                              &
                     &         self%ageTableExpansionFactor(iTime-1)                          &
                     &       )
                solutionFound=timeExceeded .and. deltaTime < turnaroundTimeTolerance*self%timeTurnaround
                if (.not.solutionFound) then
                   deltaTime=0.5d0*deltaTime
                   if (timeExceeded) then
                      self%timeTurnaround=self%timeTurnaround-deltaTime
                   else
                      self%timeTurnaround=self%timeTurnaround+deltaTime
                   end if
                end if
             end do
             self%enableRangeChecks=.true.
             self%timeMaximum      =2.0d0*self%timeTurnaround
             ! Limit the tables to the expanding part of the evolution.
             self%iTableTurnaround=iTime-2
             call Move_Alloc(self%ageTableTime           ,ageTableTimeTemporary           )
             call Move_Alloc(self%ageTableExpansionFactor,ageTableExpansionFactorTemporary)
             self%ageTableNumberPoints=self%iTableTurnaround
             allocate(self%ageTableTime           (self%ageTableNumberPoints))
             allocate(self%ageTableExpansionFactor(self%ageTableNumberPoints))
             self%ageTableTime           =ageTableTimeTemporary           (1:self%ageTableNumberPoints)
             self%ageTableExpansionFactor=ageTableExpansionFactorTemporary(1:self%ageTableNumberPoints)
             exit
          end if
       end if
    end do
    if (allocated(self%interpolatorTime)) deallocate(self%interpolatorTime)
    allocate(self%interpolatorTime)
    self%interpolatorTime=interpolator(self%ageTableExpansionFactor,self%ageTableTime)
    ! Flag that the table is now initialized.
    self%ageTableInitialized=.true.
    return
  end subroutine matterDarkEnergyLogarithmicMakeExpansionFactorTable

  double precision function matterDarkEnergyLogarithmicExpansionFactorChange(timeStart,timeEnd,expansionFactorStart)
    !!{
    Compute the expansion factor at time {\normalfont \ttfamily timeEnd} given an initial value {\normalfont \ttfamily
    expansionFactorStart} at time {\normalfont \ttfamily timeStart}.
    !!}
    use :: Numerical_ODE_Solvers, only : odeSolver
    implicit none
    double precision           , intent(in   ) :: expansionFactorStart       , timeEnd                     , &
         &                                        timeStart
    double precision           , dimension(1)  :: y
    double precision           , parameter     :: odeToleranceAbsolute=1.0d-9, odeToleranceRelative=1.0d-12
    double precision                           :: time
    type            (odeSolver)                :: solver

    time     =timeStart
    y     (1)=expansionFactorStart
    solver   =odeSolver(1_c_size_t,matterDarkEnergyLogarithmicAgeTableODEs,toleranceAbsolute=odeToleranceAbsolute,toleranceRelative=odeToleranceRelative)    
    call solver%solve(time,timeEnd,y)
    matterDarkEnergyLogarithmicExpansionFactorChange=y(1)
    return
  end function matterDarkEnergyLogarithmicExpansionFactorChange

  integer function matterDarkEnergyLogarithmicAgeTableODEs(t,a,dadt)
    !!{
    System of differential equations to solve for expansion factor vs. age.
    !!}
    use :: Interface_GSL, only : GSL_Success
    implicit none
    double precision              , intent(in   ) :: t
    double precision, dimension(:), intent(in   ) :: a
    double precision, dimension(:), intent(  out) :: dadt
    !$GLC attributes unused :: t

    if (a(1) <= 0.0d0) then
       dadt(1)=0.0d0
    else
       dadt(1)=a(1)*self_%expansionRate(a(1))
    end if
    matterDarkEnergyLogarithmicAgeTableODEs=GSL_Success
  end function matterDarkEnergyLogarithmicAgeTableODEs

  double precision function matterDarkEnergyLogarithmicTimeAtDistanceComoving(self,comovingDistance)
    !!{
    Returns the cosmological time corresponding to given {\normalfont \ttfamily comovingDistance}.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (cosmologyFunctionsMatterDarkEnergyLogarithmic), intent(inout) :: self
    double precision                                               , intent(in   ) :: comovingDistance
    !$GLC attributes unused :: self, comovingDistance

    matterDarkEnergyLogarithmicTimeAtDistanceComoving=0.0d0
    call Error_Report('functionality not implemented'//{introspection:location})
    return
  end function matterDarkEnergyLogarithmicTimeAtDistanceComoving

  double precision function matterDarkEnergyLogarithmicDistanceComoving(self,time)
    !!{
    Returns the comoving distance to cosmological time {\normalfont \ttfamily time}.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (cosmologyFunctionsMatterDarkEnergyLogarithmic), intent(inout) :: self
    double precision                                               , intent(in   ) :: time
    !$GLC attributes unused :: self, time

    matterDarkEnergyLogarithmicDistanceComoving=0.0d0
    call Error_Report('functionality not implemented'//{introspection:location})
    return
   end function matterDarkEnergyLogarithmicDistanceComoving

  double precision function matterDarkEnergyLogarithmicDistanceComovingConvert(self,output,distanceLuminosity,distanceModulus,distanceModulusKCorrected,redshift)
    !!{
    Convert between different measures of distance.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (cosmologyFunctionsMatterDarkEnergyLogarithmic), intent(inout)           :: self
    integer                                                        , intent(in   )           :: output
    double precision                                               , intent(in   ), optional :: distanceModulus, distanceModulusKCorrected, &
         &                                                                                      redshift       , distanceLuminosity
    !$GLC attributes unused :: self, output, distanceModulus, distanceModulusKCorrected, redshift, distanceLuminosity

    matterDarkEnergyLogarithmicDistanceComovingConvert=0.0d0
    call Error_Report('functionality not implemented'//{introspection:location})
    return
  end function matterDarkEnergyLogarithmicDistanceComovingConvert

  double precision function matterDarkEnergyLogarithmicEquationOfStateDarkEnergy(self,time,expansionFactor)
    !!{
    Return the dark energy equation of state.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (cosmologyFunctionsMatterDarkEnergyLogarithmic), intent(inout)           :: self
    double precision                                    , intent(in   ), optional :: expansionFactor      , time
    double precision                                                              :: expansionFactorActual

    matterDarkEnergyLogarithmicEquationOfStateDarkEnergy=self%darkEnergyEquationOfStateW0
    if (self%darkEnergyEquationOfStateW1 /= 0.0d0) then
       if (present(expansionFactor)) then
          expansionFactorActual=expansionFactor
       else if (present(time)) then
          expansionFactorActual=self%expansionFactor(time)
       else
          call Error_Report('equation of state is time dependent, but no time given'//{introspection:location})
          expansionFactorActual=1.0d0
       end if
       matterDarkEnergyLogarithmicEquationOfStateDarkEnergy=+matterDarkEnergyLogarithmicEquationOfStateDarkEnergy &
            &                                               +self%darkEnergyEquationOfStateW1                     &
            &                                               *log(expansionFactorActual)
    end if
   return
  end function matterDarkEnergyLogarithmicEquationOfStateDarkEnergy

double precision function matterDarkEnergyLogarithmicExponentDarkEnergy(self,time,expansionFactor)
    !!{
    Return the dark energy exponent for Logarithmic EOS. The xponent takes the form $\mathrm{Exponent}(a)= -3(1+w_0) - \frac{3}{2*w_1\ln{a}$.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (cosmologyFunctionsMatterDarkEnergyLogarithmic), intent(inout)           :: self
    double precision                                               , intent(in   ), optional :: expansionFactor      , time
    double precision                                                                         :: expansionFactorActual

    ! Term 1: -3(1+w0)
    matterDarkEnergyLogarithmicExponentDarkEnergy=-3.0d0*(1.0d0+self%darkEnergyEquationOfStateW0)

    if (self%darkEnergyEquationOfStateW1 /= 0.0d0) then
       if      (present(expansionFactor)) then
          expansionFactorActual=expansionFactor
       else if (present(time           )) then
          expansionFactorActual=self%expansionFactor(time)
       else
          call Error_Report('equation of state is time dependent, but no time given'//{introspection:location})
          expansionFactorActual=1.0d0
       end if
       
       ! Term 2: -1.5 * w1 * ln(a)
       ! Note: The singularity at a=1 vanishes because log(1)=0, so no division check is strictly needed.
       matterDarkEnergyLogarithmicExponentDarkEnergy = matterDarkEnergyLogarithmicExponentDarkEnergy &
            &                                          - 1.5d0 * self%darkEnergyEquationOfStateW1 * log(expansionFactorActual)
    end if
    return
  end function matterDarkEnergyLogarithmicExponentDarkEnergy

  double precision function matterDarkEnergyLogarithmicExponentDarkEnergyDerivative(self,time,expansionFactor)
    !!{
    Return the derivative of the dark energy exponent with respect to expansion factor.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (cosmologyFunctionsMatterDarkEnergyLogarithmic), intent(inout)           :: self
    double precision                                               , intent(in   ), optional :: expansionFactor      , time
    double precision                                                                         :: expansionFactorActual

    if (present(expansionFactor)) then
       expansionFactorActual=expansionFactor
    else if (present(time)) then
       expansionFactorActual=self%expansionFactor(time)
    else
       expansionFactorActual=1.0d0
    end if

    if (expansionFactorActual <= 0.0d0) then
        matterDarkEnergyLogarithmicExponentDarkEnergyDerivative = 0.0d0
    else
        matterDarkEnergyLogarithmicExponentDarkEnergyDerivative = -1.5d0 * self%darkEnergyEquationOfStateW1 / expansionFactorActual
    end if
    return
  end function matterDarkEnergyLogarithmicExponentDarkEnergyDerivative

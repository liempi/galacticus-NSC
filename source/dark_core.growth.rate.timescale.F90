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

  !+    Contributions to this file made by: Matías Liempi

  !!{
  Implementation of the \cite{...} star formation rate law for galactic \glspl{nsc}.
  !!}
  
  !![
  <darkCoreGrowthRates name="darkCoreGrowthRatesTimescale">
   <description>
    A gas inflow rate implementing the model of \citep{...} for galactic \glspl{nsc}.
    \begin{equation}
     \dot{M}_\mathrm{stellar\,BHs}^\mathrm{NSC} = \epsilon_\bullet \dot{M}_\star^\mathrm{NSC},
    \end{equation}    
    where $\epsilon_\bullet=${\normalfont \ttfamily [efficiency]} is a free parameter, and $\dot{M}_\star^\mathrm{NSC}$ is the
    star formation rate of the \glspl{nsc} component.
   </description>
  </darkCoreGrowthRates>
  !!]
  type, extends(darkCoreGrowthRatesClass) :: darkCoreGrowthRatesTimescale
     !!{
     Implementation of the \cite{....} gas inflow rate for galactic \glspl{nsc}.
     !!}
     private
     integer          :: darkCoreRadiusID            , darkCoreGasMassID  , &
       &                 darkCoreVelocityDispersionID, darkCoreTimescaleID
     double precision :: temperature
   contains
     procedure :: rate  => darkCoreTimescaleRate
  end type darkCoreGrowthRatesTimescale

  interface darkCoreGrowthRatesTimescale
     !!{
     Constructors for the {\normalfont \ttfamily ...} gas inflow rate in \glspl{nsc} class.
     !!}
     module procedure darkCoreTimescaleConstructorParameters
     module procedure darkCoreTimescaleConstructorInternal
  end interface darkCoreGrowthRatesTimescale
    
contains

  function darkCoreTimescaleConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily ...} gas inflow rate in \glspl{nsc} class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkCoreGrowthRatesTimescale)                :: self
    type            (inputParameters             ), intent(inout) :: parameters
    double precision                                              :: temperature       

    !![
    <inputParameter>
      <name>temperature</name>
      <defaultValue>100d0</defaultValue>
      <description>The temperature of the gas of the nuclear star cluster.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=darkCoreGrowthRatesTimescale(temperature)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function darkCoreTimescaleConstructorParameters

  function darkCoreTimescaleConstructorInternal(temperature) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily antonini2015} gas inflow rate from NSCs class.
    !!}
    implicit none
    type            (darkCoreGrowthRatesTimescale)                :: self
    double precision                              , intent(in   ) :: temperature       
    !![
      <constructorAssign variables="temperature"/>
      <addMetaProperty component="NSC" name="darkCoreRadius"              id="self%darkCoreRadiusID"                  isEvolvable="no"  isCreator="no"  />
      <addMetaProperty component="NSC" name="darkCoreGasMass"             id="self%darkCoreGasMassID"                 isEvolvable="no"  isCreator="yes" />
      <addMetaProperty component="NSC" name="darkCoreVelocityDispersion"  id="self%darkCoreVelocityDispersionID"      isEvolvable="no"  isCreator="yes" />
      <addMetaProperty component="NSC" name="darkCoreTimescale"           id="self%darkCoreTimescaleID"               isEvolvable="no"  isCreator="yes" />
    !!]
    return
  end function darkCoreTimescaleConstructorInternal

  double precision function darkCoreTimescaleRate(self,node)
    !!{
    Returns the stellar-mass black hole formation rate (in $M_\odot$ Gyr$^{-1}$) onto the galactic \gls{nsc} of {\normalfont \ttfamily
    node}. The rate is assumed to scale with the star formation rate of the \gls{NSC}.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentNSC
    use :: Coordinates                     , only : coordinateSpherical            , assignment(=)
    use :: Numerical_Constants_Math        , only : Pi
    use :: Ideal_Gases_Thermodynamics      , only : Ideal_Gas_Sound_Speed
    use :: Galactic_Structure_Options      , only : componentTypeNuclearStarCluster, massTypeGaseous
    use :: Numerical_Constants_Prefixes    , only : kilo
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal , megaParsec, gigayear
    use :: Mass_Distributions              , only : massDistributionClass
    implicit none
    class           (darkCoreGrowthRatesTimescale), intent(inout), target  :: self
    type            (treeNode                    ), intent(inout)          :: node
    class           (nodeComponentNSC            ),                pointer :: nuclearStarCluster
    class           (massDistributionClass       ), pointer                :: massDistributionNuclearStarCluster_ 
    type            (coordinateSpherical         )                         :: coordinates
    double precision                              , parameter              :: unnormalizedNumberOfStars = 0.079d0              , nuclearStarClusterMeanMass=0.21d0     
    double precision                                                       :: nuclearStarClusterGasFraction                    , soundSpeed                         , &
      &                                                                       frictionForceStar                                , darkCoreRadius                     , &
      &                                                                       nuclearStarClusterVelocity                       , nuclearStarClusterMassBlackHoles   , &
      &                                                                       nuclearStarClusterMassStellar                    , nuclearStarClusterEffectiveVelocity, &
      &                                                                       nuclearStarClusterNumberOfStars                  , nuclearStarClusterRadius           , &  
      &                                                                       nuclearStarClusterGasMassEnclosed                , nuclearStarClusterMassGas          , &
      &                                                                       dynamicalFrictionTimescaleGasNuclearStarCluster  , nuclearStarClusterDensityGas       , &
      &                                                                       dynamicalFrictionTimescaleStarsNuclearStarCluster, nuclearStarClusterDynamicalMass    , &
      &                                                                       dynamicalFrictionTimescale                       , darkCoreMass                       , &
      &                                                                       darkCoreVelocityDispersion                        

    ! Get the nuclear star cluster component.
    nuclearStarCluster              => node              %                      NSC(                     )
    nuclearStarClusterRadius        =  nuclearStarCluster%                   radius(                     ) !Mpc
    nuclearStarClusterMassGas       =  nuclearStarCluster%                  massGas(                     ) 
    nuclearStarClusterMassStellar   =  nuclearStarCluster%              massStellar(                     )
    nuclearStarClusterMassBlackHoles=  nuclearStarCluster%    massStellarBlackHoles(                     )
    darkCoreMass                    =  nuclearStarCluster%             massDarkCore(                     )
    darkCoreRadius                  =  nuclearStarCluster%floatRank0MetaPropertyGet(self%darkCoreRadiusID)
    nuclearStarClusterNumberOfStars =  0.36d0*(nuclearStarClusterMassStellar/unnormalizedNumberOfStars)*nuclearStarClusterMassStellar
    nuclearStarClusterDynamicalMass =  nuclearStarClusterMassStellar+nuclearStarClusterMassGas

    if  (                                         &
      &   nuclearStarClusterMassStellar   <=1.0d2 &
      &   .or.                                    &
      &   nuclearStarClusterRadius        <=0.0d0 &
      &   .or.                                    &
      &   nuclearStarClusterMassGas       <=0.0d0 &
      &   .or.                                    &
      &   nuclearStarClusterMassBlackHoles<=0.0d0 &
      &   .or.                                    &
      &   darkCoreRadius                  <=0.0d0 & 
      & ) then 
       darkCoreTimescaleRate =+0.0d0
       return 
    end if 

    ! Compute the fraction of gas in the nuclear star cluster component. 
    nuclearStarClusterGasFraction       = +nuclearStarClusterMassGas     &
      &                                   /nuclearStarClusterMassStellar
    nuclearStarClusterEffectiveVelocity = sqrt(                                 &
      &                                        +gravitationalConstant_internal  &
      &                                        *nuclearStarClusterDynamicalMass &
      &                                        /nuclearStarClusterRadius        &
      &                                       ) ! km / s
    nuclearStarClusterVelocity          = sqrt(                                 &
      &                                        +gravitationalConstant_internal  &
      &                                        *nuclearStarClusterMassStellar   &
      &                                        /nuclearStarClusterRadius        &
      &                                       )                                 &
      &                                  *    (                                 &
      &                                        +1.0d0                           &
      &                                        +nuclearStarClusterGasFraction   &
      &                                       ) !km / s

    coordinates                         =  [nuclearStarClusterRadius,0.0d0,0.0d0]
    massDistributionNuclearStarCluster_ => node                               %massDistribution    (componentTypeNuclearStarCluster,massTypeGaseous)
    nuclearStarClusterDensityGas        =  massDistributionNuclearStarCluster_%density             (coordinates                                    )
    nuclearStarClusterGasMassEnclosed   =  massDistributionNuclearStarCluster_%massEnclosedBySphere(darkCoreRadius                                 ) 
    ! We do not use this here, but I do not know where to better place this yet.
    if (darkCoreRadius>0.0d0.and.nuclearStarClusterGasMassEnclosed>0.0d0.and.darkCoreMass>0.0d0) then
      PRINT *,darkCoreMass, darkCoreRadius, nuclearStarClusterGasMassEnclosed
      darkCoreVelocityDispersion          = sqrt(                                     &
        &                                         gravitationalConstant_internal      &
        &                                        *darkCoreMass                        &
        &                                        *(                                   &
        &                                          +darkCoreMass                      &
        &                                          +nuclearStarClusterGasMassEnclosed &
        &                                         )                                   &
        &                                        /(                                   &
        &                                          +0.5d0                             &
        &                                          *nuclearStarClusterGasMassEnclosed &
        &                                          *darkCoreRadius                    &
        &                                         )                                   &         
        &                                       )
    else 
      darkCoreVelocityDispersion= 0.0d0
    end if 
    call nuclearStarCluster%floatRank0MetaPropertySet(self%darkCoreVelocityDispersionID,darkCoreVelocityDispersion       )
    call nuclearStarCluster%floatRank0MetaPropertySet(self%darkCoreGasMassID           ,nuclearStarClusterGasMassEnclosed)
    !![
       <objectDestructor name="massDistributionNuclearStarCluster_"/>
    !!]    
    soundSpeed        = Ideal_Gas_Sound_Speed(self%temperature)

    ! Estimates the friction force acting on one star due to the gas as in Schleicher et al. 2020 (https://ui.adsabs.harvard.edu/abs/2022MNRAS.512.6192S/abstract)
    ! The force is in units of M⊙ km s⁻².
    frictionForceStar = +0.11d0*4.0d0*Pi                         & ! Adimensional
      &                 *(gravitationalConstant_internal**2.0d0) & ! Mpc² M⊙⁻² (km s⁻¹)⁴
      &                 *nuclearStarClusterDensityGas            & ! M⊙ Mpc⁻3
      &                 *(nuclearStarClusterMeanMass**2.0d0)     & ! M⊙²
      &                 /(soundSpeed**2.0d0)                     & ! km² s⁻²
      &                 /(megaParsec/kilo)                         ! Convert 1 Mpc to km.
    ! Determinate the gas contribution to the dynamical friction timescale (Eq. 16, Schleicher et al. 2020)
    dynamicalFrictionTimescaleGasNuclearStarCluster   =+0.5d0                                         &
      &                                                *nuclearStarClusterMassGas                     & ! M⊙
      &                                                *(nuclearStarClusterEffectiveVelocity**2.0d0)  & ! (km s⁻¹)²
      &                                                /(                                             &
      &                                                   nuclearStarClusterNumberOfStars             & ! Adimensional 
      &                                                  *frictionForceStar                           & ! M⊙ (km s⁻¹)²
      &                                                  *nuclearStarClusterVelocity                  & ! km s⁻¹
      &                                                  *gigayear                                    & ! Convert from seconds to Gyr.
      &                                                 )
    ! Determinate the stellar contribution to the dynamical friction timescale. (Eq. 17, Schleicher et al. 2020)
    dynamicalFrictionTimescaleStarsNuclearStarCluster =+nuclearStarClusterVelocity        & ! km s⁻¹
      &                                                *nuclearStarClusterMeanMass        & ! M⊙
      &                                                /(                                 &
      &                                                   frictionForceStar               & ! M⊙ km s⁻²
      &                                                  *gigayear                        & ! Convert to Gyr.
      &                                                 )
    ! Compute the dynamical friction timescale in Gyr
    dynamicalFrictionTimescale                        = (                                                             &
      &                                                  +dynamicalFrictionTimescaleGasNuclearStarCluster  **(-1.0d0) &
      &                                                  +dynamicalFrictionTimescaleStarsNuclearStarCluster**(-1.0d0) &
      &                                                 )**(-1.0d0)
   
    call nuclearStarCluster%floatRank0MetaPropertySet(self%darkCoreTimescaleID,dynamicalFrictionTimescale)

    if (dynamicalFrictionTimescale <= 0.0d0) then
      darkCoreTimescaleRate =+0.0d0
    else
      darkCoreTimescaleRate = nuclearStarClusterMassBlackHoles &
        &                    /dynamicalFrictionTimescale
    end if 
    return
  end function darkCoreTimescaleRate


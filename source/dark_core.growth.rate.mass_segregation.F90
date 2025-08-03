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

  use :: Star_Formation_Rates_Nuclear_Star_Clusters, only : starFormationRateNuclearStarClustersClass

  
  !![
  <darkCoreGrowthRates name="darkCoreGrowthRatesMassSegregation">
   <description>
    A dark core mass rate where the mass evolution takes place over a timescale for galactic \glspl{nsc}.
    \begin{equation}
     \dot{M}_\mathrm{stellar\,BHs}^\mathrm{NSC} = \epsilon_\bullet \dot{M}_\star^\mathrm{NSC},
    \end{equation}    
    where $\epsilon_\bullet=${\normalfont \ttfamily [efficiency]} is a free parameter, and $\dot{M}_\star^\mathrm{NSC}$ is the
    star formation rate of the \glspl{nsc} component.
   </description>
  </darkCoreGrowthRates>
  !!]
  type, extends(darkCoreGrowthRatesClass) :: darkCoreGrowthRatesMassSegregation
     !!{
     Implementation of the \cite{....} gas inflow rate for galactic \glspl{nsc}.
     !!}
     private
     class(starFormationRateNuclearStarClustersClass), pointer :: starFormationRateNuclearStarClusters_ => null()

     integer                                                    :: darkCoreRadiusID                               , darkCoreGasMassID          , &
       &                                                           darkCoreVelocityDispersionID                   , darkCoreMassSegregationID  , &
       &                                                           nuclearStarClusterNumberOfStarsID              , nuclearStarClusterDensityID, &
       &                                                           darkCoreTimescaleID
     double precision                                           :: efficiencyBlackHoleFormation                   , boostFactorIMF              , &
       &                                                           fractionBlackHoles
   contains
     final     ::          darkCoreMassSegregationDestructor
     procedure :: rate  => darkCoreMassSegregationRate
  end type darkCoreGrowthRatesMassSegregation

  interface darkCoreGrowthRatesMassSegregation
     !!{
     Constructors for the {\normalfont \ttfamily ...} gas inflow rate in \glspl{nsc} class.
     !!}
     module procedure darkCoreMassSegregationConstructorParameters
     module procedure darkCoreMassSegregationConstructorInternal
  end interface darkCoreGrowthRatesMassSegregation
    
contains

  function darkCoreMassSegregationConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily ...} gas inflow rate in \glspl{nsc} class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkCoreGrowthRatesMassSegregation       )                :: self
    type            (inputParameters                          ), intent(inout) :: parameters
    class           (starFormationRateNuclearStarClustersClass), pointer       :: starFormationRateNuclearStarClusters_
    double precision                                                           :: efficiencyBlackHoleFormation         , boostFactorIMF, &
       &                                                                          fractionBlackHoles

    !![
    <inputParameter>
      <name>efficiencyBlackHoleFormation</name>
      <defaultValue>1.6d-3</defaultValue>
      <description>Efficiency of the stellar mass black hole production in nuclear star clusters.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>boostFactorIMF</name>
      <defaultValue>1.0d0</defaultValue>
      <description>Boost factor to enhance the production of stellar mass black holes in the nuclear star clusters asumming a top-heavy initial mass function </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>fractionBlackHoles</name>
      <defaultValue>0.01d0</defaultValue>
      <description>Boost factor to enhance the production of stellar mass black holes in the nuclear star clusters asumming a top-heavy initial mass function.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="starFormationRateNuclearStarClusters" name="starFormationRateNuclearStarClusters_" source="parameters"/>
    !!]
    self=darkCoreGrowthRatesMassSegregation(efficiencyBlackHoleFormation,boostFactorIMF,fractionBlackHoles,starFormationRateNuclearStarClusters_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="starFormationRateNuclearStarClusters_"/>
    !!]
    return
  end function darkCoreMassSegregationConstructorParameters

  function darkCoreMassSegregationConstructorInternal(efficiencyBlackHoleFormation,boostFactorIMF,fractionBlackHoles,starFormationRateNuclearStarClusters_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily antonini2015} gas inflow rate from NSCs class.
    !!}
    implicit none
    type            (darkCoreGrowthRatesMassSegregation       )                        :: self
    class           (starFormationRateNuclearStarClustersClass), intent(in   ), target :: starFormationRateNuclearStarClusters_
    double precision                                           , intent(in   )         :: efficiencyBlackHoleFormation         , boostFactorIMF, &
       &                                                                                  fractionBlackHoles
    !![
      <constructorAssign variables="efficiencyBlackHoleFormation,boostFactorIMF,fractionBlackHoles,*starFormationRateNuclearStarClusters_"/>
      <addMetaProperty component="NSC" name="darkCoreRadius"                  id="self%darkCoreRadiusID"                  isEvolvable="no"  isCreator="no" />
      <addMetaProperty component="NSC" name="darkCoreGasMass"                 id="self%darkCoreGasMassID"                 isEvolvable="no"  isCreator="yes"/>
      <addMetaProperty component="NSC" name="darkCoreVelocityDispersion"      id="self%darkCoreVelocityDispersionID"      isEvolvable="no"  isCreator="yes"/>
      <addMetaProperty component="NSC" name="darkCoreTimescale"               id="self%darkCoreTimescaleID"               isEvolvable="no"  isCreator="yes"/>
      <addMetaProperty component="NSC" name="nuclearStarClusterNumberOfStars" id="self%nuclearStarClusterNumberOfStarsID" isEvolvable="no"  isCreator="yes"/> 
      <addMetaProperty component="NSC" name="nuclearStarClusterDensity"       id="self%nuclearStarClusterDensityID"       isEvolvable="no"  isCreator="yes"/> 
    !!]
    return
  end function darkCoreMassSegregationConstructorInternal

  subroutine darkCoreMassSegregationDestructor(self)
    !!{
    Destructor for the \refClass{nodePropertyExtractorStarFormationRate} property extractor class.
    !!}
    implicit none
    type (darkCoreGrowthRatesMassSegregation), intent(inout) :: self
  
    !![
    <objectDestructor name="self%starFormationRateNuclearStarClusters_"/>
    !!]
    return
  end subroutine darkCoreMassSegregationDestructor

  double precision function darkCoreMassSegregationRate(self,node)
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
    class           (darkCoreGrowthRatesMassSegregation), intent(inout), target  :: self
    type            (treeNode                          ), intent(inout)          :: node
    class           (nodeComponentNSC                  ),                pointer :: nuclearStarCluster
    class           (massDistributionClass             ), pointer                :: massDistributionNuclearStarCluster_ 
    type            (coordinateSpherical               )                         :: coordinates
    double precision                                    , parameter              :: unnormalizedNumberOfStars = 0.079d0, nuclearStarClusterMeanMass=0.21d0     
    double precision                                                             :: nuclearStarClusterGasFraction      , darkCoreRadius                     , &
      &                                                                             nuclearStarClusterMassBlackHoles   , darkCoreVelocityDispersion         , &
      &                                                                             nuclearStarClusterMassStellar      , darkCoreMass                       , &
      &                                                                             nuclearStarClusterNumberOfStars    , nuclearStarClusterRadius           , &  
      &                                                                             nuclearStarClusterGasMassEnclosed  , nuclearStarClusterMassGas          , &
      &                                                                             nuclearStarClusterDynamicalMass    , nuclearStarClusterStarFormationRate, &
      &                                                                             nuclearStarClusterDensityGas

    ! Get the nuclear star cluster component.
    nuclearStarCluster                 => node              %                      NSC(                     )
    nuclearStarClusterRadius           =  nuclearStarCluster%                   radius(                     ) !Mpc
    nuclearStarClusterMassGas          =  nuclearStarCluster%                  massGas(                     ) 
    nuclearStarClusterMassStellar      =  nuclearStarCluster%              massStellar(                     )
    darkCoreMass                       =  nuclearStarCluster%             massDarkCore(                     )
    darkCoreRadius                     =  nuclearStarCluster%floatRank0MetaPropertyGet(self%darkCoreRadiusID)
    nuclearStarClusterNumberOfStars    =  0.38d0*(nuclearStarClusterMassStellar/unnormalizedNumberOfStars)
    nuclearStarClusterDynamicalMass    =  nuclearStarClusterMassStellar+nuclearStarClusterMassGas
    nuclearStarClusterStarFormationRate=  self%starFormationRateNuclearStarClusters_%rate(node)
    if  (                                         &
      &   nuclearStarClusterMassStellar   <=0.0d0 &
      &   .or.                                    &
      &   nuclearStarClusterMassGas       < 0.0d0 &
      & ) then 
       darkCoreMassSegregationRate =+0.0d0
       return 
    end if 

    ! Compute the fraction of gas in the nuclear star cluster component. 
    nuclearStarClusterGasFraction       = +nuclearStarClusterMassGas     &
      &                                   /nuclearStarClusterMassStellar

    coordinates                         =  [nuclearStarClusterRadius,0.0d0,0.0d0]
    massDistributionNuclearStarCluster_ => node                               %massDistribution    (componentTypeNuclearStarCluster,massTypeGaseous)
    nuclearStarClusterDensityGas        =  massDistributionNuclearStarCluster_%density             (coordinates                                    )
    nuclearStarClusterGasMassEnclosed   =  massDistributionNuclearStarCluster_%massEnclosedBySphere(darkCoreRadius                                 ) 
    ! We do not use this here, but I do not know where to better place this yet.
    if ( darkCoreRadius                  >0.0d0 &
      & .and.                                   &
      & darkCoreMass                     >0.0d0 &
      & .and.                                   &
      & nuclearStarClusterGasMassEnclosed>0.0d0 &
      & ) then
      darkCoreVelocityDispersion          = sqrt(                                      &
        &                                        (                                     &
        &                                          gravitationalConstant_internal      &
        &                                         *darkCoreMass                        &
        &                                         *(                                   &
        &                                           +darkCoreMass                      &
        &                                           +nuclearStarClusterGasMassEnclosed &
        &                                          )                                   &
        &                                         )                                    &
        &                                         /(                                   &
        &                                            darkCoreMass                      &
        &                                           *darkCoreRadius                    &
        &                                          )                                   &         
        &                                        )

    else 
      darkCoreVelocityDispersion= 0.0d0
    end if
    call nuclearStarCluster%floatRank0MetaPropertySet(self%nuclearStarClusterNumberOfStarsID,nuclearStarClusterNumberOfStars  )
    call nuclearStarCluster%floatRank0MetaPropertySet(self%darkCoreVelocityDispersionID     ,darkCoreVelocityDispersion       )
    call nuclearStarCluster%floatRank0MetaPropertySet(self%darkCoreGasMassID                ,nuclearStarClusterGasMassEnclosed)
    !![
       <objectDestructor name="massDistributionNuclearStarCluster_"/>
    !!]    

    call nuclearStarCluster%floatRank0MetaPropertySet(self%darkCoreTimescaleID, 0.0d0)
    if ( nuclearStarClusterStarFormationRate <= 0.0d0) then
      darkCoreMassSegregationRate =+0.0d0
    else
      darkCoreMassSegregationRate = nuclearStarClusterStarFormationRate*self%efficiencyBlackHoleFormation*self%boostFactorIMF*self%fractionBlackHoles
    end if 
    return
  end function darkCoreMassSegregationRate


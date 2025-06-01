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
  <darkCoreGrowthRates name="darkCoreGrowthRatesCoreCollapse">
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
  type, extends(darkCoreGrowthRatesClass) :: darkCoreGrowthRatesCoreCollapse
     !!{
     Implementation of the \cite{....} gas inflow rate for galactic \glspl{nsc}.
     !!}
     private
     integer          :: darkCoreRadiusID                 , darkCoreGasMassID          , &
       &                 darkCoreVelocityDispersionID     , darkCoreTimescaleID        , &
       &                 nuclearStarClusterNumberOfStarsID, nuclearStarClusterDensityID
   contains
     procedure :: rate  => darkCoreCoreCollapseRate
  end type darkCoreGrowthRatesCoreCollapse

  interface darkCoreGrowthRatesCoreCollapse
     !!{
     Constructors for the {\normalfont \ttfamily ...} gas inflow rate in \glspl{nsc} class.
     !!}
     module procedure darkCoreCoreCollapseConstructorParameters
     module procedure darkCoreCoreCollapseConstructorInternal
  end interface darkCoreGrowthRatesCoreCollapse
    
contains

  function darkCoreCoreCollapseConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily ...} gas inflow rate in \glspl{nsc} class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkCoreGrowthRatesCoreCollapse)                :: self
    type            (inputParameters                ), intent(inout) :: parameters

    self=darkCoreGrowthRatesCoreCollapse()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function darkCoreCoreCollapseConstructorParameters

  function darkCoreCoreCollapseConstructorInternal() result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily antonini2015} gas inflow rate from NSCs class.
    !!}
    implicit none
    type (darkCoreGrowthRatesCoreCollapse) :: self
    !![
      <addMetaProperty component="NSC" name="darkCoreRadius"                  id="self%darkCoreRadiusID"                  isEvolvable="no"  isCreator="no" />
      <addMetaProperty component="NSC" name="darkCoreGasMass"                 id="self%darkCoreGasMassID"                 isEvolvable="no"  isCreator="yes"/>
      <addMetaProperty component="NSC" name="darkCoreVelocityDispersion"      id="self%darkCoreVelocityDispersionID"      isEvolvable="no"  isCreator="yes"/>
      <addMetaProperty component="NSC" name="darkCoreTimescale"               id="self%darkCoreTimescaleID"               isEvolvable="no"  isCreator="yes"/>
      <addMetaProperty component="NSC" name="nuclearStarClusterNumberOfStars" id="self%nuclearStarClusterNumberOfStarsID" isEvolvable="no"  isCreator="yes"/> 
      <addMetaProperty component="NSC" name="nuclearStarClusterDensity"       id="self%nuclearStarClusterDensityID"       isEvolvable="no"  isCreator="yes"/> 
    !!]
    return
  end function darkCoreCoreCollapseConstructorInternal

  double precision function darkCoreCoreCollapseRate(self,node)
    !!{
    Returns the stellar-mass black hole formation rate (in $M_\odot$ Gyr$^{-1}$) onto the galactic \gls{nsc} of {\normalfont \ttfamily
    node}. The rate is assumed to scale with the star formation rate of the \gls{NSC}.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentNSC
    use :: Galactic_Structure_Options      , only : componentTypeNuclearStarCluster        , massTypeGaseous            , massTypeStellar
    use :: Coordinates                     , only : coordinateSpherical                    , assignment(=)
    use :: Nuclear_Star_Clusters_Utilities , only : nuclearStarClusterCoreCollapseTimescale
    use :: Mass_Distributions              , only : massDistributionClass                  , kinematicsDistributionClass
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal         , megaParsec, gigayear
    implicit none
    class           (darkCoreGrowthRatesCoreCollapse), intent(inout), target  :: self
    type            (treeNode                       ), intent(inout)          :: node
    class           (nodeComponentNSC               ),                pointer :: nuclearStarCluster
    class           (massDistributionClass          ), pointer                :: massDistributionStellarNuclearStarCluster_, massDistributionGasNuclearStarCluster_ 
    class           (kinematicsDistributionClass    ), pointer                :: kinematicsNuclearStarClusterStellar_
    type            (coordinateSpherical            )                         :: coordinates
    double precision                                 , parameter              :: unnormalizedNumberOfStars = 0.079d0, nuclearStarClusterMeanMass=0.21d0     
    double precision                                                          :: darkCoreRadius                     , darkCoreMass                        , &
      &                                                                          nuclearStarClusterMassBlackHoles   , nuclearStarClusterMassStellar       , &
      &                                                                          nuclearStarClusterRadius           , nuclearStarClusterGasMassEnclosed   , &
      &                                                                          nuclearStarClusterMassGas          , darkCoreVelocityDispersion          , &
      &                                                                          coreCollapseTimescale              , nuclearStarClusterDensityGas        , &
      &                                                                          nuclearStarClusterNumberOfStars    , velocityDispersionNuclearStarCluster, &
      &                                                                          nuclearStarClusterVelDisp          , nuclearStarClusterStellarMassEnclosed
      
      ! Get the nuclear star cluster component.
    nuclearStarCluster              => node              %                      NSC(                     )
    nuclearStarClusterRadius        =  nuclearStarCluster%                   radius(                     ) !Mpc
    nuclearStarClusterMassGas       =  nuclearStarCluster%                  massGas(                     ) 
    nuclearStarClusterMassStellar   =  nuclearStarCluster%              massStellar(                     )
    nuclearStarClusterMassBlackHoles=  nuclearStarCluster%    massStellarBlackHoles(                     )
    darkCoreMass                    =  nuclearStarCluster%             massDarkCore(                     )
    darkCoreRadius                  =  nuclearStarCluster%floatRank0MetaPropertyGet(self%darkCoreRadiusID)
    nuclearStarClusterNumberOfStars =  0.38d0*(nuclearStarClusterMassStellar/unnormalizedNumberOfStars)

    if  (                                         &
      &   nuclearStarClusterMassStellar   <=0.0d0 &
      &   .or.                                    &
      &   nuclearStarClusterMassGas       < 0.0d0 &
      &   .or.                                    &
      &   nuclearStarClusterMassBlackHoles<=0.0d0 &
      & ) then 
       darkCoreCoreCollapseRate =+0.0d0
       return 
    end if 

    if (       nuclearStarClusterRadius>0.0d0 &
      &   .and.          darkCoreRadius>0.0d0 &
      &   .and.            darkCoreMass>0.0d0 &
      & ) then

      coordinates                                =  [nuclearStarClusterRadius,0.0d0,0.0d0]
      massDistributionGasNuclearStarCluster_     => node                                  %massDistribution    (componentTypeNuclearStarCluster,massTypeGaseous)
      nuclearStarClusterDensityGas               =  massDistributionGasNuclearStarCluster_%density             (coordinates                                    )
      nuclearStarClusterGasMassEnclosed          =  massDistributionGasNuclearStarCluster_%massEnclosedBySphere(darkCoreRadius                                 ) 
    
      massDistributionStellarNuclearStarCluster_ => node                                      %massDistribution      (componentTypeNuclearStarCluster,massTypeStellar)
      kinematicsNuclearStarClusterStellar_       => massDistributionStellarNuclearStarCluster_%kinematicsDistribution(                                               )
      nuclearStarClusterStellarMassEnclosed      =  massDistributionStellarNuclearStarCluster_%massEnclosedBySphere  (darkCoreRadius                                 ) 
      nuclearStarClusterVelDisp                  = sqrt(gravitationalConstant_internal*nuclearStarClusterStellarMassEnclosed/darkCoreRadius)*(1+nuclearStarClusterGasMassEnclosed/nuclearStarClusterStellarMassEnclosed)

      coordinates                                =  [darkCoreRadius,0.0d0,0.0d0]
      velocityDispersionNuclearStarCluster       =  kinematicsNuclearStarClusterStellar_      %velocityDispersion1D  (coordinates, massDistributionStellarNuclearStarCluster_, &
        &                                                                                                                          massDistributionStellarNuclearStarCluster_   )
      !![
        <objectDestructor name="massDistributionGasNuclearStarCluster_"/>
        <objectDestructor name="massDistributionStellarNuclearStarCluster_"/>
        <objectDestructor name="kinematicsNuclearStarClusterStellar_"/>
      !!]
    else
       nuclearStarClusterVelDisp = 0.0d0
    end if

    if (nuclearStarClusterVelDisp>0.0d0) then
      darkCoreVelocityDispersion = 26.0d0*nuclearStarClusterVelDisp
    else 
      darkCoreVelocityDispersion = 0.0d0
    end if 

    coreCollapseTimescale = nuclearStarClusterCoreCollapseTimescale(                               &
      &                                                             nuclearStarClusterRadius     , &
      &                                                             nuclearStarClusterMassStellar, &
      &                                                             nuclearStarClusterMassGas    , &
      &                                                             .true.                         &
      &                                                            )
    
    call nuclearStarCluster%floatRank0MetaPropertySet(self%darkCoreTimescaleID              , coreCollapseTimescale)
    call nuclearStarCluster%floatRank0MetaPropertySet(self%nuclearStarClusterNumberOfStarsID, nuclearStarClusterNumberOfStars  )
    call nuclearStarCluster%floatRank0MetaPropertySet(self%nuclearStarClusterDensityID      , nuclearStarClusterDensityGas     )
    call nuclearStarCluster%floatRank0MetaPropertySet(self%darkCoreVelocityDispersionID     , darkCoreVelocityDispersion       )
    call nuclearStarCluster%floatRank0MetaPropertySet(self%darkCoreGasMassID                , nuclearStarClusterGasMassEnclosed)
    

    if (coreCollapseTimescale <= 0.0d0) then
      darkCoreCoreCollapseRate =+0.0d0
    else
      darkCoreCoreCollapseRate = nuclearStarClusterMassBlackHoles &
        &                    /coreCollapseTimescale
    end if 
    return
  end function darkCoreCoreCollapseRate


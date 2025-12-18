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
  Implementation of the \cite{...} dark core formation rate law for galactic \glspl{nsc}.
  !!}

  use :: Star_Formation_Rates_Nuclear_Star_Clusters, only : starFormationRateNuclearStarClustersClass
  
  !![
  <darkCoreGrowthRates name="darkCoreGrowthRatesMassSegregation">
   <description>
    A dark core mass rate where the mass evolution takes place over a timescale for galactic \glspl{nsc}.
    \begin{equation}
     \dot{M}_\mathrm{dark core}^\mathrm{NSC} = \epsilon_\bullet f_\mathrm{IMF}^\mathrm{boost} f_\mathrm{} \dot{M}_\star^\mathrm{NSC},
    \end{equation}    
    where $\epsilon_\bullet=${\normalfont \ttfamily [efficiency]} is a free parameter, and $\dot{M}_\star^\mathrm{NSC}$ is the
    star formation rate of the \glspl{nsc} component.
   </description>
  </darkCoreGrowthRates>
  !!]
  type, extends(darkCoreGrowthRatesClass) :: darkCoreGrowthRatesMassSegregation
     !!{
     Implementation of a dark core formation rate law for galactic \glspl{nsc}.
     !!}
     private
     class(starFormationRateNuclearStarClustersClass), pointer :: starFormationRateNuclearStarClusters_ => null()
     double precision                                          :: efficiencyBlackHoleFormation                   , boostFactorIMF   , &
       &                                                          fractionBlackHoles                             , contractionFactor
   contains
     final     ::          darkCoreMassSegregationDestructor
     procedure :: rate  => darkCoreMassSegregationRate
  end type darkCoreGrowthRatesMassSegregation

  interface darkCoreGrowthRatesMassSegregation
     !!{
     Constructors for the {\normalfont \ttfamily darkCoreGrowthRate} class.
     !!}
     module procedure darkCoreMassSegregationConstructorParameters
     module procedure darkCoreMassSegregationConstructorInternal
  end interface darkCoreGrowthRatesMassSegregation
    
contains

  function darkCoreMassSegregationConstructorParameters(parameters) result(self)
    !!{
    Constructor for the dark core formation rate law class which takes a parameter set as input.
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
      <description>Efficiency of the stellar-mass black holes production in nuclear star clusters.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>boostFactorIMF</name>
      <defaultValue>1.0d0</defaultValue>
      <description>Boost factor to enhance the production of stellar-mass black holes in the nuclear star clusters asumming a top-heavier initial mass function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>fractionBlackHoles</name>
      <defaultValue>0.01d0</defaultValue>
      <description>Factor that takes into account already placed stellar-mass black holes in the center of the nuclear star cluster.</description>
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
    double precision                                           , intent(in   )         :: efficiencyBlackHoleFormation         , boostFactorIMF   , &
       &                                                                                  fractionBlackHoles
    !![
      <constructorAssign variables="efficiencyBlackHoleFormation,boostFactorIMF,fractionBlackHoles,*starFormationRateNuclearStarClusters_"/>
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
    Returns the dark core formation rate (in $M_\odot$ Gyr$^{-1}$) in the galactic \gls{nsc} of {\normalfont \ttfamily
    node}. The rate is assumed to scale with the star formation rate of the \gls{NSC}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentNSC
    implicit none
    class           (darkCoreGrowthRatesMassSegregation), intent(inout), target  :: self
    type            (treeNode                          ), intent(inout)          :: node
    class           (nodeComponentNSC                  ),                pointer :: nuclearStarCluster
    double precision                                                             :: nuclearStarClusterStarFormationRate
    
    darkCoreMassSegregationRate =+0.0d0

    ! Get the nuclear star cluster component.
    nuclearStarCluster                 => node                                      %        NSC(    )
    nuclearStarClusterStarFormationRate=  self%starFormationRateNuclearStarClusters_%       rate(node)
    
    if (nuclearStarClusterStarFormationRate<=0.0d0) return

    darkCoreMassSegregationRate =+self%efficiencyBlackHoleFormation        & 
      &                          *self%boostFactorIMF                      &
      &                          *self%fractionBlackHoles                  &
      &                          *     nuclearStarClusterStarFormationRate
    return
  end function darkCoreMassSegregationRate


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
  <nuclearStarClusterStellarBlackHoleGrowthRates name="nuclearStarClusterStellarBlackHoleGrowthRatesScaled">
   <description>
    A gas inflow rate implementing the model of \citep{...} for galactic \glspl{nsc}.
    \begin{equation}
     \dot{M}_\mathrm{stellar\,BHs}^\mathrm{NSC} = \epsilon_\bullet \dot{M}_\star^\mathrm{NSC},
    \end{equation}    
    where $\epsilon_\bullet=${\normalfont \ttfamily [efficiency]} is a free parameter, and $\dot{M}_\star^\mathrm{NSC}$ is the
    star formation rate of the \glspl{nsc} component.
   </description>
  </nuclearStarClusterStellarBlackHoleGrowthRates>
  !!]
  type, extends(nuclearStarClusterStellarBlackHoleGrowthRatesClass) :: nuclearStarClusterStellarBlackHoleGrowthRatesScaled
     !!{
     Implementation of the \cite{....} gas inflow rate for galactic \glspl{nsc}.
     !!}
     private
     class           (starFormationRateNuclearStarClustersClass), pointer :: starFormationRateNuclearStarClusters_ => null()
     double precision                                                     :: efficiency
     contains
     final     ::          scaledDestructor
     procedure :: rate  => scaledRate
  end type nuclearStarClusterStellarBlackHoleGrowthRatesScaled

  interface nuclearStarClusterStellarBlackHoleGrowthRatesScaled
     !!{
     Constructors for the {\normalfont \ttfamily ...} gas inflow rate in \glspl{nsc} class.
     !!}
     module procedure scaledConstructorParameters
     module procedure scaledConstructorInternal
  end interface nuclearStarClusterStellarBlackHoleGrowthRatesScaled
    
contains

  function scaledConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily ...} gas inflow rate in \glspl{nsc} class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nuclearStarClusterStellarBlackHoleGrowthRatesScaled)                :: self
    type            (inputParameters                                    ), intent(inout) :: parameters
    class           (starFormationRateNuclearStarClustersClass          ), pointer       :: starFormationRateNuclearStarClusters_
    double precision                                                                     :: efficiency       

    !![
    <inputParameter>
      <name>efficiency</name>
      <defaultSource>\citep{antonini_coevolution_2015}</defaultSource>
      <defaultValue>1.0d-2</defaultValue>
      <description>Parameter controlling the rate of the gas inflow onto the nuclear star cluster.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="starFormationRateNuclearStarClusters" name="starFormationRateNuclearStarClusters_" source="parameters"/>
    !!]
    self=nuclearStarClusterStellarBlackHoleGrowthRatesScaled(efficiency,starFormationRateNuclearStarClusters_)
    !![
    <inputParametersValidate source="parameters"        />
    <objectDestructor name="starFormationRateNuclearStarClusters_"/>
    !!]
    return
  end function scaledConstructorParameters

  function scaledConstructorInternal(efficiency,starFormationRateNuclearStarClusters_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily antonini2015} gas inflow rate from NSCs class.
    !!}
    implicit none
    type            (nuclearStarClusterStellarBlackHoleGrowthRatesScaled)                        :: self
    class           (starFormationRateNuclearStarClustersClass          ), intent(in   ), target :: starFormationRateNuclearStarClusters_
    double precision                                                     , intent(in   )         :: efficiency
    !![
    <constructorAssign variables="efficiency, *starFormationRateNuclearStarClusters_"/>
    !!]
    return
  end function scaledConstructorInternal

  subroutine scaledDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily NuclearStarClusterGrowth} class
    !!}
    implicit none
    type(nuclearStarClusterStellarBlackHoleGrowthRatesScaled), intent(inout) :: self
    !![
    <objectDestructor name="self%starFormationRateNuclearStarClusters_"/>
    !!]
    return
  end subroutine scaledDestructor

  double precision function scaledRate(self,node) result(rate)
    !!{
    Returns the stellar-mass black hole formation rate (in $M_\odot$ Gyr$^{-1}$) onto the galactic \gls{nsc} of {\normalfont \ttfamily
    node}. The rate is assumed to scale with the star formation rate of the \gls{NSC}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentNSC
    implicit none
    class           (nuclearStarClusterStellarBlackHoleGrowthRatesScaled), intent(inout), target  :: self
    type            (treeNode                                           ), intent(inout)          :: node
    class           (nodeComponentNSC                                   ),                pointer :: nuclearStarCluster
    double precision                                                                              :: rateStarFormationNuclearStarCluster

    ! Get the nuclear star cluster component.
    nuclearStarCluster                  => node%                                      NSC (   )
    ! Get the star formation rate of the nuclear star cluster component.
    rateStarFormationNuclearStarCluster =  self%starFormationRateNuclearStarClusters_%rate(node)  
    ! Find the rate of mass accretion onto the nuclear star cluster.
    if (rateStarFormationNuclearStarCluster <= 0.0d0) then
      rate    =+0.0d0
    else
       rate   =+self%efficiency                          &
            &  *     rateStarFormationNuclearStarCluster
    end if 
    return
  end function scaledRate

  

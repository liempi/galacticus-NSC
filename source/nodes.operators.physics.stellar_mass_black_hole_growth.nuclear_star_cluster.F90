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
  Implements a node operator class that performs black hole formation in nuclear star cluster.
  !!}

  use :: Nuclear_Star_Cluster_Stellar_Black_Holes_Growth_Rates, only : nuclearStarClusterStellarBlackHoleGrowthRatesClass
  !![
  <nodeOperator name="nodeOperatorStellarBlackHoleGrowthNuclearStarClusters">
   <description>A node operator class that performs stellar-mass black hole formation in \gsl{NSCs}.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorStellarBlackHoleGrowthNuclearStarClusters
     !!{
     A node operator class that performs the black hole evolution in dark cores.
     !!}
     private
     class  (nuclearStarClusterStellarBlackHoleGrowthRatesClass), pointer :: nuclearStarClusterStellarBlackHoleGrowthRates_ => null()
   contains
     final     ::                          stellarMassblackHoleGrowthDestructor
     procedure :: differentialEvolution => stellarMassBlackHoleGrowthDifferentialEvolution
  end type nodeOperatorStellarBlackHoleGrowthNuclearStarClusters
  
  interface nodeOperatorStellarBlackHoleGrowthNuclearStarClusters
     !!{
     Constructors for the {\normalfont \ttfamily blackHoleFormationDarkCore} node operator class.
     !!}
     module procedure stellarMassBlackHoleGrowthConstructorParameters
     module procedure stellarMassBlackHoleGrowthConstructorInternal
  end interface nodeOperatorStellarBlackHoleGrowthNuclearStarClusters
  
contains

  function stellarMassBlackHoleGrowthConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily starFormation} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorStellarBlackHoleGrowthNuclearStarClusters)                :: self
    type (inputParameters                                      ), intent(inout) :: parameters
    class(nuclearStarClusterStellarBlackHoleGrowthRatesClass   ), pointer       :: nuclearStarClusterStellarBlackHoleGrowthRates_
    !![

    <objectBuilder class="nuclearStarClusterStellarBlackHoleGrowthRates" name="nuclearStarClusterStellarBlackHoleGrowthRates_" source="parameters"/>
    !!]
    self=nodeOperatorStellarBlackHoleGrowthNuclearStarClusters(nuclearStarClusterStellarBlackHoleGrowthRates_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="nuclearStarClusterStellarBlackHoleGrowthRates_"/>
    !!]
    return
  end function stellarMassBlackHoleGrowthConstructorParameters

  function stellarMassBlackHoleGrowthConstructorInternal(nuclearStarClusterStellarBlackHoleGrowthRates_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily blackHoleFormationDarkCore} node operator class.
    !!}
    implicit none
    type   (nodeOperatorStellarBlackHoleGrowthNuclearStarClusters)                        :: self
    class  (nuclearStarClusterStellarBlackHoleGrowthRatesClass   ), intent(in   ), target :: nuclearStarClusterStellarBlackHoleGrowthRates_
    !![
    <constructorAssign variables="*nuclearStarClusterStellarBlackHoleGrowthRates_"/>
    !!]
    return
  end function stellarMassBlackHoleGrowthConstructorInternal

  subroutine stellarMassblackHoleGrowthDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily blackHoleFormationDarkCore} node operator class.
    !!}
    implicit none
    type (nodeOperatorStellarBlackHoleGrowthNuclearStarClusters), intent(inout) :: self

    !![
    <objectDestructor name="self%nuclearStarClusterStellarBlackHoleGrowthRates_"/>
    !!]
    return
  end subroutine stellarMassblackHoleGrowthDestructor
  
  subroutine stellarMassBlackHoleGrowthDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Performs the stellar-mass black Hole growth in nuclear star clusters.
    !!}
    use :: Galacticus_Nodes    , only : interruptTask, nodeComponentNSC, propertyInactive, treeNode
    use :: Histories           , only : operator(*)  , history
    use :: Abundances_Structure, only : operator(*)
    implicit none
    class           (nodeOperatorStellarBlackHoleGrowthNuclearStarClusters), intent(inout), target  :: self
    type            (treeNode                                             ), intent(inout), target  :: node
    logical                                                                , intent(inout)          :: interrupt
    procedure       (interruptTask                                        ), intent(inout), pointer :: functionInterrupt
    integer                                                                , intent(in   )          :: propertyType
    class           (nodeComponentNSC                                     )               , pointer :: nuclearStarCluster
    double precision                                                                                :: rateStellarMassBlackHoleGrowth
    type            (history                                              )                         :: historyTransferRate
    
    ! Return immediately if inactive property is requested.
    if (propertyInactive(propertyType)) return
    ! Check for a nuclear star cluster, return immediately if nuclear star cluster is unphysical.
    nuclearStarCluster => node%NSC()
    if (nuclearStarCluster%massStellar() <= 0.0d0) return 
    ! Get the stellar-mass black hole growth rate, return immediately if the rate is negative or zero.
    rateStellarMassBlackHoleGrowth = self%nuclearStarClusterStellarBlackHoleGrowthRates_%rate(node) 
    if (rateStellarMassBlackHoleGrowth <= 0.0d0) return
    ! Adjust quantities.
    call nuclearStarCluster%          massStellarRate(                                        &
        &                                             -rateStellarMassBlackHoleGrowth         &
        &                                            )
    call nuclearStarCluster%massStellarBlackHolesRate(                                        &
        &                                             +rateStellarMassBlackHoleGrowth         &
        &                                            )
    call nuclearStarCluster%    abundancesStellarRate(                                        &
        &                                             -rateStellarMassBlackHoleGrowth         & 
        &                                             *nuclearStarCluster%abundancesStellar() &
        &                                             /nuclearStarCluster%massStellar()       &
        &                                            )
    historyTransferRate=nuclearStarCluster%stellarPropertiesHistory()
    if (historyTransferRate%exists()) &
        & call nuclearStarCluster%stellarPropertiesHistoryRate (-rateStellarMassBlackHoleGrowth &
             &                                                  *historyTransferRate            &
             &                                                 )
    call historyTransferRate%destroy()
    !! Star formation history.
    historyTransferRate=nuclearStarCluster%starFormationHistory()
    if (historyTransferRate%exists()) &
        & call nuclearStarCluster%starFormationHistoryRate     (-rateStellarMassBlackHoleGrowth &
             &                                                  *historyTransferRate            &
             &                                                 )
    call historyTransferRate%destroy()
    return
  end subroutine stellarMassBlackHoleGrowthDifferentialEvolution
  
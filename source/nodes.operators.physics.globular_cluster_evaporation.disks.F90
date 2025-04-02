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
  Implements a node operator class that performs star formation in disks.
  !!}

  use :: Globular_Cluster_Evaporation_Rates_Disks, only : globularClusterEvaporationRateDisksClass

  !![
  <nodeOperator name="nodeOperatorGlobularClusterEvaporationDisks">
   <description>A node operator class that performs star formation.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorGlobularClusterEvaporationDisks
     !!{
     A node operator class that performs star formation.
     !!}
     private
     class  (globularClusterEvaporationRateDisksClass), pointer :: globularClusterEvaporationRateDisks_ => null()
     integer                                                    :: globularClusterStellarMassDiskID 
   contains
     final     ::                                   globularClusterEvaporationDisksDestructor
     procedure :: differentialEvolution          => globularClusterEvaporationDisksDifferentialEvolution
  end type nodeOperatorGlobularClusterEvaporationDisks
  
  interface nodeOperatorGlobularClusterEvaporationDisks
     !!{
     Constructors for the {\normalfont \ttfamily globularClusterFormationDisks} node operator class.
     !!}
     module procedure globularClusterEvaporationDisksConstructorParameters
     module procedure globularClusterEvaporationDisksConstructorInternal
  end interface nodeOperatorGlobularClusterEvaporationDisks
  
contains

  function globularClusterEvaporationDisksConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily globularClusterFormation} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorGlobularClusterEvaporationDisks)                :: self
    type (inputParameters                            ), intent(inout) :: parameters
    class(globularClusterEvaporationRateDisksClass   ), pointer       :: globularClusterEvaporationRateDisks_

    !![
    <objectBuilder class="globularClusterEvaporationRateDisks" name="globularClusterEvaporationRateDisks_" source="parameters"/>
    !!]
    self=nodeOperatorGlobularClusterEvaporationDisks(globularClusterEvaporationRateDisks_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="globularClusterEvaporationRateDisks_"/>
    !!]
    return
  end function globularClusterEvaporationDisksConstructorParameters

  function globularClusterEvaporationDisksConstructorInternal(globularClusterEvaporationRateDisks_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily globularClusterFormationDisks} node operator class.
    !!}
    implicit none
    type (nodeOperatorGlobularClusterEvaporationDisks)                        :: self
    class(globularClusterEvaporationRateDisksClass   ), intent(in   ), target :: globularClusterEvaporationRateDisks_

    !![
    <constructorAssign variables="*globularClusterEvaporationRateDisks_"/>
    !!]
    !![
    <addMetaProperty component="disk" name="globularClusterStellarMassDisk" id="self%globularClusterStellarMassDiskID" isEvolvable="yes" isCreator="no"/>
    !!]
    return
  end function globularClusterEvaporationDisksConstructorInternal

  subroutine globularClusterEvaporationDisksDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily globularClusterFormationDisks} node operator class.
    !!}
    implicit none
    type(nodeOperatorGlobularClusterEvaporationDisks), intent(inout) :: self

    !![
    <objectDestructor name="self%globularClusterEvaporationRateDisks_"/>
    !!]
    return
  end subroutine globularClusterEvaporationDisksDestructor
  
  subroutine globularClusterEvaporationDisksDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Perform globular cluster evaporation  in a disk.
    !!}
    use :: Galacticus_Nodes              , only : propertyInactive, propertyTypeActive, propertyEvaluate, nodeComponentDisk
    implicit none
    class           (nodeOperatorGlobularClusterEvaporationDisks), intent(inout), target  :: self
    type            (treeNode                                   ), intent(inout), target  :: node
    logical                                                      , intent(inout)          :: interrupt
    procedure       (interruptTask                              ), intent(inout), pointer :: functionInterrupt
    integer                                                      , intent(in   )          :: propertyType
    class           (nodeComponentDisk                          )               , pointer :: disk
    double precision                                                                      :: rateGlobularClusterEvaporation
    ! Check for a realistic disk, return immediately if disk is unphysical.
    disk => node%disk()
    if     (     disk%angularMomentum() <= 0.0d0      &
         &  .or. disk%radius         () <= 0.0d0      &
         &  .or. disk%massGas        () <= 0.0d0      &
	       &  .or. disk%massStellar    () <= 0.0d0      &
         & ) return
    if (propertyInactive(propertyType)) return

    rateGlobularClusterEvaporation=self%globularClusterEvaporationRateDisks_%rate(node)

    if (rateGlobularClusterEvaporation <= 0.0d0) return

    call disk%massStellarRate           (                                        & 
          &                              +rateGlobularClusterEvaporation         &
          &                             )

    call disk%floatRank0MetaPropertyRate(                                        &
          &                               self%globularClusterStellarMassDiskID, &
          &                              -rateGlobularClusterEvaporation         &
          &                             )
    return
  end subroutine globularClusterEvaporationDisksDifferentialEvolution
  

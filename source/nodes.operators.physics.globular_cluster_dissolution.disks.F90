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

  use :: Globular_Cluster_Dissolution_Rates_Disks, only : globularClusterDissolutionRateDisksClass

  !![
  <nodeOperator name="nodeOperatorglobularClusterDissolutionDisks">
   <description>A node operator class that performs star formation.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorglobularClusterDissolutionDisks
     !!{
     A node operator class that performs star formation.
     !!}
     private
     class  (globularClusterDissolutionRateDisksClass), pointer :: globularClusterDissolutionRateDisks_ => null()
     integer                                                    :: globularClusterStellarMassDiskID 
   contains
     final     ::                          globularClusterDissolutionDisksDestructor
     procedure :: differentialEvolution => globularClusterDissolutionDisksDifferentialEvolution
  end type nodeOperatorglobularClusterDissolutionDisks
  
  interface nodeOperatorglobularClusterDissolutionDisks
     !!{
     Constructors for the {\normalfont \ttfamily globularClusterFormationDisks} node operator class.
     !!}
     module procedure globularClusterDissolutionDisksConstructorParameters
     module procedure globularClusterDissolutionDisksConstructorInternal
  end interface nodeOperatorglobularClusterDissolutionDisks
  
contains

  function globularClusterDissolutionDisksConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily globularClusterFormation} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorglobularClusterDissolutionDisks)                :: self
    type (inputParameters                            ), intent(inout) :: parameters
    class(globularClusterDissolutionRateDisksClass   ), pointer       :: globularClusterDissolutionRateDisks_

    !![
    <objectBuilder class="globularClusterDissolutionRateDisks" name="globularClusterDissolutionRateDisks_" source="parameters"/>
    !!]
    self=nodeOperatorglobularClusterDissolutionDisks(globularClusterDissolutionRateDisks_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="globularClusterDissolutionRateDisks_"/>
    !!]
    return
  end function globularClusterDissolutionDisksConstructorParameters

  function globularClusterDissolutionDisksConstructorInternal(globularClusterDissolutionRateDisks_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily globularClusterFormationDisks} node operator class.
    !!}
    implicit none
    type (nodeOperatorglobularClusterDissolutionDisks)                        :: self
    class(globularClusterDissolutionRateDisksClass   ), intent(in   ), target :: globularClusterDissolutionRateDisks_

    !![
    <constructorAssign variables="*globularClusterDissolutionRateDisks_"/>
    !!]
    !![
    <addMetaProperty component="disk" name="globularClusterStellarMassDisk" id="self%globularClusterStellarMassDiskID" isEvolvable="yes" isCreator="no"/>
    !!]
    return
  end function globularClusterDissolutionDisksConstructorInternal

  subroutine globularClusterDissolutionDisksDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily globularClusterFormationDisks} node operator class.
    !!}
    implicit none
    type(nodeOperatorglobularClusterDissolutionDisks), intent(inout) :: self

    !![
    <objectDestructor name="self%globularClusterDissolutionRateDisks_"/>
    !!]
    return
  end subroutine globularClusterDissolutionDisksDestructor
  
  subroutine globularClusterDissolutionDisksDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Perform globular cluster dissolution in a disk.
    !!}
    use :: Galacticus_Nodes, only : propertyInactive, propertyTypeActive, propertyEvaluate, nodeComponentDisk
    implicit none
    class           (nodeOperatorglobularClusterDissolutionDisks), intent(inout), target  :: self
    type            (treeNode                                   ), intent(inout), target  :: node
    logical                                                      , intent(inout)          :: interrupt
    procedure       (interruptTask                              ), intent(inout), pointer :: functionInterrupt
    integer                                                      , intent(in   )          :: propertyType
    class           (nodeComponentDisk                          )               , pointer :: disk
    double precision                                                                      :: globularClusterMass, rateGlobularClusterDissolution
    ! Check for a realistic disk, return immediately if disk is unphysical.
    disk => node%disk()
    if     (         disk%angularMomentum() <= 0.0d0 &
         &       .or.disk%radius         () <= 0.0d0 &
         &       .or.disk%massGas        () <= 0.0d0 &
         &       .or.disk%massStellar    () <= 0.0d0 &
         & ) return
    if (propertyInactive(propertyType)) return

    globularClusterMass           = disk%floatRank0MetaPropertyGet(self%globularClusterStellarMassDiskID)
    rateGlobularClusterDissolution=self%globularClusterDissolutionRateDisks_%rate(node)

    if (rateGlobularClusterDissolution<=0.0d0.or.globularClusterMass<=0.0d0) return

    call disk%massStellarRate           (+rateGlobularClusterDissolution)

    call disk%floatRank0MetaPropertyRate(                                        &
          &                               self%globularClusterStellarMassDiskID, &
          &                                   -rateGlobularClusterDissolution    &
          &                              )
    return
  end subroutine globularClusterDissolutionDisksDifferentialEvolution

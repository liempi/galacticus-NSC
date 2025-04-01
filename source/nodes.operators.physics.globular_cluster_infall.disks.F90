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

  use :: Globular_Cluster_Infall_Rates_Disks, only : globularClusterInfallRateDisksClass

  !![
  <nodeOperator name="nodeOperatorglobularClusterInfallDisks">
   <description>A node operator class that performs star formation.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorglobularClusterInfallDisks
     !!{
     A node operator class that performs star formation.
     !!}
     private
     class  (globularClusterInfallRateDisksClass), pointer :: globularClusterInfallRateDisks_ => null()
     integer                                               :: globularClusterStellarMassDiskID
     double precision                                      :: fraction 
   contains
     final     ::                          globularClusterInfallDisksDestructor
     procedure :: differentialEvolution => globularClusterInfallDisksDifferentialEvolution
  end type nodeOperatorglobularClusterInfallDisks
  
  interface nodeOperatorglobularClusterInfallDisks
     !!{
     Constructors for the {\normalfont \ttfamily globularClusterFormationDisks} node operator class.
     !!}
     module procedure globularClusterInfallDisksConstructorParameters
     module procedure globularClusterInfallDisksConstructorInternal
  end interface nodeOperatorglobularClusterInfallDisks
  
contains

  function globularClusterInfallDisksConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily globularClusterFormation} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorglobularClusterInfallDisks)                :: self
    type (inputParameters                       ), intent(inout) :: parameters
    class(globularClusterInfallRateDisksClass   ), pointer       :: globularClusterInfallRateDisks_
    double precision                                             :: fraction

    !![
    <inputParameter>
      <name>fraction</name>
      <defaultValue>0.5d0</defaultValue>
      <source>parameters</source>
      <description>Fraction of the the globular clusters that survives to merge with the nuclear star cluster.</description>
    </inputParameter>
    <objectBuilder class="globularClusterInfallRateDisks" name="globularClusterInfallRateDisks_" source="parameters"/>
    !!]
    self=nodeOperatorglobularClusterInfallDisks(fraction,globularClusterInfallRateDisks_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="globularClusterInfallRateDisks_"/>
    !!]
    return
  end function globularClusterInfallDisksConstructorParameters

  function globularClusterInfallDisksConstructorInternal(fraction,globularClusterInfallRateDisks_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily globularClusterFormationDisks} node operator class.
    !!}
    implicit none
    type (nodeOperatorglobularClusterInfallDisks)                        :: self
    class(globularClusterInfallRateDisksClass   ), intent(in   ), target :: globularClusterInfallRateDisks_
    double precision                             , intent(in  )          :: fraction
    !![
    <constructorAssign variables="fraction,*globularClusterInfallRateDisks_"/>
    !!]
    !![
    <addMetaProperty component="disk" name="globularClusterStellarMassDisk" id="self%globularClusterStellarMassDiskID" isEvolvable="yes" isCreator="no"/>
    !!]
    return
  end function globularClusterInfallDisksConstructorInternal

  subroutine globularClusterInfallDisksDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily globularClusterFormationDisks} node operator class.
    !!}
    implicit none
    type(nodeOperatorglobularClusterInfallDisks), intent(inout) :: self

    !![
    <objectDestructor name="self%globularClusterInfallRateDisks_"/>
    !!]
    return
  end subroutine globularClusterInfallDisksDestructor
  
  subroutine globularClusterInfallDisksDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Perform globular cluster dissolution in a disk.
    !!}
    use :: Galacticus_Nodes, only : propertyInactive, propertyTypeActive      , propertyEvaluate, nodeComponentDisk, &
            &                       nodeComponentNSC, nodeComponentNSCStandard
    implicit none
    class           (nodeOperatorglobularClusterInfallDisks), intent(inout), target  :: self
    type            (treeNode                              ), intent(inout), target  :: node
    logical                                                 , intent(inout)          :: interrupt
    procedure       (interruptTask                         ), intent(inout), pointer :: functionInterrupt
    integer                                                 , intent(in   )          :: propertyType
    class           (nodeComponentDisk                     )               , pointer :: disk
    class           (nodeComponentNSC                      )               , pointer :: nuclearStarCluster
    double precision                                                                 :: rateGlobularClusterInfall
    ! Check for a realistic disk, return immediately if disk is unphysical.
    disk => node%disk()
    if     (         disk%angularMomentum() <= 0.0d0 &
         &       .or.disk%radius         () <= 0.0d0 &
         &       .or.disk%massGas        () <= 0.0d0 &
         &       .or.disk%massStellar    () <= 0.0d0 &
         & ) return
    if (propertyInactive(propertyType)) return

    rateGlobularClusterInfall=self%fraction*self%globularClusterInfallRateDisks_%rate(node)

    if (rateGlobularClusterInfall <= 0.0d0) return

    ! Get nuclear star cluster component
    nuclearStarCluster => node%NSC()

    ! Detect nuclear star cluster component type.
    select type (nuclearStarCluster)
    type is (nodeComponentNSC)
       ! Generic type - interrupt and create a nuclear star cluster.
       interrupt         =  .true.
       functionInterrupt => nuclearStarClusterCreate
       return
    class default
       ! A nuclear star cluster exists - continue processing.
       ! Remove gas from the spheroid component and add to the nuclear star cluster component.
       call disk             %floatRank0MetaPropertyRate(                                       &
          &                                              self%globularClusterStellarMassDiskID, &
          &                                                  -rateGlobularClusterInfall         &
          &                                             )
       ! WARNING, here we need to adjust stellar histories in both, disk and nuclear star cluster components.
       call nuclearStarCluster%massStellarRate(+rateGlobularClusterInfall)
    end select
    return
  end subroutine globularClusterInfallDisksDifferentialEvolution

  subroutine nuclearStarClusterCreate(node,timeEnd)
    !!{
    Creates the nuclear star cluster via interrupt.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentNSC, treeNode
    implicit none
    type             (treeNode        ), intent(inout), target  :: node
    double precision                   , intent(in   ), optional:: timeEnd
    class            (nodeComponentNSC),                pointer :: nuclearStarCluster
    !$GLC attributes unused :: timeEnd

    nuclearStarCluster => node%NSC(autoCreate=.true.)
    return 
  end subroutine nuclearStarClusterCreate

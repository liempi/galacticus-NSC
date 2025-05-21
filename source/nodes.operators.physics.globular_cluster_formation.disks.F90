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

  use :: Globular_Cluster_Formation_Rates_Disks, only : globularClusterFormationRateDisksClass
  use :: Satellite_Merging_Mass_Movements      , only : mergerMassMovementsClass

  !![
  <nodeOperator name="nodeOperatorGlobularClusterFormationDisks">
   <description>A node operator class that performs star formation.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorGlobularClusterFormationDisks
     !!{
     A node operator class that performs star formation.
     !!}
     private
     class  (globularClusterFormationRateDisksClass), pointer :: globularClusterFormationRateDisks_ => null()
     class  (mergerMassMovementsClass              ), pointer :: mergerMassMovements_               => null()
     integer                                                  :: globularClusterStellarMassDiskID            , globularClusterStellarMassSpheroidID
   contains
     final     ::                                   globularClusterFormationDisksDestructor
     procedure :: differentialEvolutionScales    => globularClusterFormationDisksDifferentialEvolutionScales
     procedure :: differentialEvolutionInactives => globularClusterFormationDisksDifferentialEvolutionInactives
     procedure :: differentialEvolution          => globularClusterFormationDisksDifferentialEvolution
     procedure :: galaxiesMerge                  => globularClusterFormationDisksGalaxiesMerge
  end type nodeOperatorGlobularClusterFormationDisks
  
  interface nodeOperatorGlobularClusterFormationDisks
     !!{
     Constructors for the {\normalfont \ttfamily globularClusterFormationDisks} node operator class.
     !!}
     module procedure globularClusterFormationDisksConstructorParameters
     module procedure globularClusterFormationDisksConstructorInternal
  end interface nodeOperatorGlobularClusterFormationDisks
  
contains

  function globularClusterFormationDisksConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily globularClusterFormation} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorGlobularClusterFormationDisks)                :: self
    type (inputParameters                          ), intent(inout) :: parameters
    class(globularClusterFormationRateDisksClass   ), pointer       :: globularClusterFormationRateDisks_
    class(mergerMassMovementsClass                 ), pointer       :: mergerMassMovements_

    !![
    <objectBuilder class="globularClusterFormationRateDisks" name="globularClusterFormationRateDisks_" source="parameters"/>
    <objectBuilder class="mergerMassMovements"               name="mergerMassMovements_"               source="parameters"/>
    !!]
    self=nodeOperatorGlobularClusterFormationDisks(globularClusterFormationRateDisks_,mergerMassMovements_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="globularClusterFormationRateDisks_"/>
    <objectDestructor name="mergerMassMovements_"              />
    !!]
    return
  end function globularClusterFormationDisksConstructorParameters

  function globularClusterFormationDisksConstructorInternal(globularClusterFormationRateDisks_,mergerMassMovements_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily globularClusterFormationDisks} node operator class.
    !!}
    implicit none
    type (nodeOperatorGlobularClusterFormationDisks)                        :: self
    class(globularClusterFormationRateDisksClass   ), intent(in   ), target :: globularClusterFormationRateDisks_
    class(mergerMassMovementsClass                 ), intent(in   ), target :: mergerMassMovements_

    !![
    <constructorAssign variables="*globularClusterFormationRateDisks_"/>
    <constructorAssign variables="*mergerMassMovements_"/>
    !!]
    !![
    <addMetaProperty component="disk"     name="globularClusterStellarMassDisk"     id="self%globularClusterStellarMassDiskID"     isEvolvable="yes" isCreator="yes"/>
    <addMetaProperty component="spheroid" name="globularClusterStellarMassSpheroid" id="self%globularClusterStellarMassSpheroidID" isEvolvable="yes" isCreator="no"/>
    !!]
    return
  end function globularClusterFormationDisksConstructorInternal

  subroutine globularClusterFormationDisksDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily globularClusterFormationDisks} node operator class.
    !!}
    implicit none
    type(nodeOperatorGlobularClusterFormationDisks), intent(inout) :: self

    !![
    <objectDestructor name="self%globularClusterFormationRateDisks_"/>
    <objectDestructor name="self%mergerMassMovements_"              />
    !!]
    return
  end subroutine globularClusterFormationDisksDestructor
  
  subroutine globularClusterFormationDisksDifferentialEvolutionInactives(self,node)
    !!{
    Mark disk as inactive for ODE solving.    
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk
    implicit none
    class(nodeOperatorGlobularClusterFormationDisks), intent(inout) :: self
    type (treeNode                                 ), intent(inout) :: node
    class(nodeComponentDisk                        ), pointer       :: disk
    
    ! Get disk component.
    disk => node%disk()
    ! Mark as inactive.
    select type (disk)
    type is (nodeComponentDisk)
       ! Disk does not yet exist - nothing to do here.
    class default
       call disk%floatRank0MetaPropertyInactive(self%globularClusterStellarMassDiskID)
    end select
    return
  end subroutine globularClusterFormationDisksDifferentialEvolutionInactives

  subroutine globularClusterFormationDisksDifferentialEvolutionScales(self,node)
    !!{
    Set absolute ODE solver scale for the unweighted and time-weighted stellar masses.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk
    implicit none
    class           (nodeOperatorGlobularClusterFormationDisks), intent(inout) :: self
    type            (treeNode                                 ), intent(inout) :: node
    class           (nodeComponentDisk                        ), pointer       :: disk
    double precision                                           , parameter     :: massMinimum=1.0d+0, scaleRelative=1.0d-3
    double precision                                                           :: mass
    ! Get disk component
    disk =>  node%disk()
    ! Set scale for masses.
    mass = disk%massGas()+disk%massStellar()
    ! Set scales.
    select type (disk)
    type is (nodeComponentDisk)
       ! Disk does not yet exist - nothing to do here.
    class default
       call disk%floatRank0MetaPropertyScale(self%globularClusterStellarMassDiskID, massMinimum)
    end select
    return
  end subroutine globularClusterFormationDisksDifferentialEvolutionScales
  
  subroutine globularClusterFormationDisksDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Perform globular cluster formation in a disk.
    !!}
    use :: Galacticus_Nodes, only : propertyInactive, propertyTypeActive, propertyEvaluate, nodeComponentDisk
    implicit none
    class           (nodeOperatorGlobularClusterFormationDisks), intent(inout), target  :: self
    type            (treeNode                                 ), intent(inout), target  :: node
    logical                                                    , intent(inout)          :: interrupt
    procedure       (interruptTask                            ), intent(inout), pointer :: functionInterrupt
    integer                                                    , intent(in   )          :: propertyType
    class           (nodeComponentDisk                        )               , pointer :: disk
    double precision                                                                    :: rateGlobularClusterFormation
    ! Check for a realistic disk, return immediately if disk is unphysical.
    disk => node%disk()
    if     (     disk%angularMomentum() < 0.0d0      &
         &  .or. disk%radius         () < 0.0d0      &
         &  .or. disk%massGas        () < 0.0d0      &
         & ) return
    if (propertyInactive(propertyType)) return

    rateGlobularClusterFormation=self%globularClusterFormationRateDisks_%rate(node)

    if (rateGlobularClusterFormation <= 0.0d0) return

    call disk%massStellarRate           (                                        & 
          &                              -rateGlobularClusterFormation           &
          &                             )

    call disk%floatRank0MetaPropertyRate(                                        &
          &                               self%globularClusterStellarMassDiskID, &
          &                              +rateGlobularClusterFormation           &
          &                             )
    return
  end subroutine globularClusterFormationDisksDifferentialEvolution
  
  subroutine globularClusterFormationDisksGalaxiesMerge(self,node)
    !!{
    Combine integrals of star formation rate when galaxies merge.
    !!}
    use :: Error                           , only : Error_Report
    use :: Galacticus_Nodes                , only : nodeComponentDisk    , nodeComponentSpheroid
    use :: Satellite_Merging_Mass_Movements, only : destinationMergerDisk, destinationMergerSpheroid, destinationMergerUnmoved, enumerationDestinationMergerType
    implicit none
    class  (nodeOperatorGlobularClusterFormationDisks), intent(inout) :: self
    type   (treeNode                                 ), intent(inout) :: node
    type   (treeNode                                 ), pointer       :: nodeHost
    class  (nodeComponentDisk                        ), pointer       :: disk                   , diskHost
    class  (nodeComponentSpheroid                    ), pointer       :: spheroid               , spheroidHost
    type   (enumerationDestinationMergerType         )                :: destinationGasSatellite, destinationStarsSatellite, &
         &                                                               destinationGasHost     , destinationStarsHost
    logical                                                           :: mergerIsMajor
    ! Find the node to merge with.
    nodeHost               => node    %mergesWith(                 )
    disk                   => node    %disk      (autoCreate=.true.)
    spheroid               => node    %spheroid  (autoCreate=.true.)
    diskHost               => nodeHost%disk      (autoCreate=.true.)
    spheroidHost           => nodeHost%spheroid  (autoCreate=.true.)

    ! Get mass movement descriptors.
    call self%mergerMassMovements_%get(node,destinationGasSatellite,destinationStarsSatellite,destinationGasHost,destinationStarsHost,mergerIsMajor)
    ! Move star formation rates within the host if necessary.
    select case (destinationStarsHost%ID)
    case (destinationMergerDisk    %ID)
       call diskHost                     %floatRank0MetaPropertySet(                                        self%globularClusterStellarMassDiskID     , &
            &                                                       +diskHost    %floatRank0MetaPropertyGet(self%globularClusterStellarMassDiskID    )  &
            &                                                       +spheroidHost%floatRank0MetaPropertyGet(self%globularClusterStellarMassSpheroidID)  &
            &                                                      )
    case (destinationMergerSpheroid%ID)
       call spheroidHost                 %floatRank0MetaPropertySet(                                        self%globularClusterStellarMassSpheroidID , &
            &                                                       +diskHost    %floatRank0MetaPropertyGet(self%globularClusterStellarMassDiskID    )  &
            &                                                       +spheroidHost%floatRank0MetaPropertyGet(self%globularClusterStellarMassSpheroidID)  & 
            &                                                      )
    case (destinationMergerUnmoved %ID)
       ! Do nothing.
    case default
       call Error_Report('unrecognized movesTo descriptor'//{introspection:location})
    end select

    select case (destinationStarsSatellite%ID)
    case (destinationMergerDisk    %ID)
       call diskHost    %floatRank0MetaPropertySet(                                              self%globularClusterStellarMassDiskID     , &
            &                                      +diskHost          %floatRank0MetaPropertyGet(self%globularClusterStellarMassDiskID    )  &
            &                                      +disk              %floatRank0MetaPropertyGet(self%globularClusterStellarMassDiskID    )  &
            &                                      +spheroid          %floatRank0MetaPropertyGet(self%globularClusterStellarMassSpheroidID)  &
            &                                     )
    case (destinationMergerSpheroid%ID)
      call spheroidHost%floatRank0MetaPropertySet(                                               self%globularClusterStellarMassSpheroidID , &
            &                                      +spheroidHost      %floatRank0MetaPropertyGet(self%globularClusterStellarMassSpheroidID)  &
            &                                      +disk              %floatRank0MetaPropertyGet(self%globularClusterStellarMassDiskID    )  &
            &                                      +spheroid          %floatRank0MetaPropertyGet(self%globularClusterStellarMassSpheroidID)  &
            &                                     )
    case default
       call Error_Report('unrecognized movesTo descriptor'//{introspection:location})
    end select
    ! Zero rates in the secondary,
    call    disk                     %floatRank0MetaPropertySet(                                 self%globularClusterStellarMassDiskID     , &
         &                                                      +0.0d0                                                                       &
         &                                                     )
    call    spheroid                 %floatRank0MetaPropertySet(                                 self%globularClusterStellarMassSpheroidID , &
         &                                                      +0.0d0                                                                       &
         &                                                     )
    return
  end subroutine globularClusterFormationDisksGalaxiesMerge
  

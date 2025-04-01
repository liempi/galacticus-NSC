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

  use :: Globular_Cluster_Formation_Rates_Spheroids, only : globularClusterFormationRateSpheroidsClass
  use :: Satellite_Merging_Mass_Movements          , only : mergerMassMovementsClass

  !![
  <nodeOperator name="nodeOperatorGlobularClusterFormationSpheroids">
   <description>A node operator class that performs star formation.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorGlobularClusterFormationSpheroids
     !!{
     A node operator class that performs star formation.
     !!}
     private
     class  (globularClusterFormationRateSpheroidsClass), pointer :: globularClusterFormationRateSpheroids_ => null()
     class  (mergerMassMovementsClass                  ), pointer :: mergerMassMovements_                   => null()
     integer                                                      :: globularClusterStellarMassDiskID                , globularClusterStellarMassSpheroidID
   contains
     final     ::                                   globularClusterFormationSpheroidsDestructor
     procedure :: differentialEvolutionScales    => globularClusterFormationSpheroidsDifferentialEvolutionScales
     procedure :: differentialEvolutionInactives => globularClusterFormationSpheroidsDifferentialEvolutionInactives
     procedure :: differentialEvolution          => globularClusterFormationSpheroidsDifferentialEvolution
     procedure :: galaxiesMerge                  => globularClusterFormationSpheroidsGalaxiesMerge
  end type nodeOperatorGlobularClusterFormationSpheroids
  
  interface nodeOperatorGlobularClusterFormationSpheroids
     !!{
     Constructors for the {\normalfont \ttfamily globularClusterFormationSpheroids} node operator class.
     !!}
     module procedure globularClusterFormationSpheroidsConstructorParameters
     module procedure globularClusterFormationSpheroidsConstructorInternal
  end interface nodeOperatorGlobularClusterFormationSpheroids
  
contains

  function globularClusterFormationSpheroidsConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily globularClusterFormation} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorGlobularClusterFormationSpheroids)                :: self
    type (inputParameters                              ), intent(inout) :: parameters
    class(globularClusterFormationRateSpheroidsClass   ), pointer       :: globularClusterFormationRateSpheroids_
    class(mergerMassMovementsClass                     ), pointer       :: mergerMassMovements_

    !![
    <objectBuilder class="globularClusterFormationRateSpheroids" name="globularClusterFormationRateSpheroids_" source="parameters"/>
    <objectBuilder class="mergerMassMovements"                   name="mergerMassMovements_"                   source="parameters"/>
    !!]
    self=nodeOperatorGlobularClusterFormationSpheroids(globularClusterFormationRateSpheroids_,mergerMassMovements_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="globularClusterFormationRateSpheroids_"/>
    <objectDestructor name="mergerMassMovements_"              />
    !!]
    return
  end function globularClusterFormationSpheroidsConstructorParameters

  function globularClusterFormationSpheroidsConstructorInternal(globularClusterFormationRateSpheroids_,mergerMassMovements_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily globularClusterFormationSpheroids} node operator class.
    !!}
    implicit none
    type (nodeOperatorGlobularClusterFormationSpheroids)                        :: self
    class(globularClusterFormationRateSpheroidsClass   ), intent(in   ), target :: globularClusterFormationRateSpheroids_
    class(mergerMassMovementsClass                     ), intent(in   ), target :: mergerMassMovements_

    !![
    <constructorAssign variables="*globularClusterFormationRateSpheroids_"/>
    <constructorAssign variables="*mergerMassMovements_"/>
    !!]
    !![
    <addMetaProperty component="disk"     name="globularClusterStellarMassDisk"     id="self%globularClusterStellarMassDiskID"     isEvolvable="yes" isCreator="no" />
    <addMetaProperty component="spheroid" name="globularClusterStellarMassSpheroid" id="self%globularClusterStellarMassSpheroidID" isEvolvable="yes" isCreator="yes"/>
    !!]
    return
  end function globularClusterFormationSpheroidsConstructorInternal

  subroutine globularClusterFormationSpheroidsDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily globularClusterFormationSpheroids} node operator class.
    !!}
    implicit none
    type(nodeOperatorGlobularClusterFormationSpheroids), intent(inout) :: self

    !![
    <objectDestructor name="self%globularClusterFormationRateSpheroids_"/>
    <objectDestructor name="self%mergerMassMovements_"                  />
    !!]
    return
  end subroutine globularClusterFormationSpheroidsDestructor
  
  subroutine globularClusterFormationSpheroidsDifferentialEvolutionInactives(self,node)
    !!{
    Mark disk as inactive for ODE solving.    
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSpheroid
    implicit none
    class(nodeOperatorGlobularClusterFormationSpheroids), intent(inout) :: self
    type (treeNode                                     ), intent(inout) :: node
    class(nodeComponentSpheroid                        ), pointer       :: spheroid
    
    ! Get spheroid component.
    spheroid => node%spheroid()
    ! Mark as inactive.
    select type (spheroid)
    type is (nodeComponentSpheroid)
       ! Disk does not yet exist - nothing to do here.
    class default
       call spheroid%floatRank0MetaPropertyInactive(self%globularClusterStellarMassSpheroidID)
    end select
    return
  end subroutine globularClusterFormationSpheroidsDifferentialEvolutionInactives

  subroutine globularClusterFormationSpheroidsDifferentialEvolutionScales(self,node)
    !!{
    Set absolute ODE solver scale for the unweighted and time-weighted stellar masses.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSpheroid
    implicit none
    class           (nodeOperatorGlobularClusterFormationSpheroids), intent(inout) :: self
    type            (treeNode                                     ), intent(inout) :: node
    class           (nodeComponentSpheroid                        ), pointer       :: spheroid
    double precision                                               , parameter     :: massMinimum=1.0d+0, scaleRelative=1.0d-3
    double precision                                                               :: mass
    ! Get spheroid component
    spheroid =>  node%spheroid()
    ! Set scale for masses.
    mass = spheroid%massGas()+spheroid%massStellar()
    ! Set scales.
    select type (spheroid)
    type is (nodeComponentSpheroid)
       ! Spheroid does not yet exist - nothing to do here.
    class default
       call spheroid%floatRank0MetaPropertyScale(self%globularClusterStellarMassSpheroidID, max(scaleRelative*mass,massMinimum))
    end select
    return
  end subroutine globularClusterFormationSpheroidsDifferentialEvolutionScales

  subroutine globularClusterFormationSpheroidsDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Perform globular cluster formation in a spheroid.
    !!}
    use :: Galacticus_Nodes, only : propertyInactive, propertyTypeActive, propertyEvaluate, nodeComponentSpheroid
    implicit none
    class           (nodeOperatorGlobularClusterFormationSpheroids), intent(inout), target  :: self
    type            (treeNode                                     ), intent(inout), target  :: node
    logical                                                        , intent(inout)          :: interrupt
    procedure       (interruptTask                                ), intent(inout), pointer :: functionInterrupt
    integer                                                        , intent(in   )          :: propertyType
    class           (nodeComponentSpheroid                        )               , pointer :: spheroid
    double precision                                                                        :: rateGlobularClusterFormation
    ! Check for a realistic disk, return immediately if disk is unphysical.
    spheroid =>  node    %spheroid()
    if  (        spheroid%angularMomentum() < 0.0d0      &
         &  .or. spheroid%radius         () < 0.0d0      &
         &  .or. spheroid%massGas        () < 0.0d0      &
         & ) return
    if (propertyInactive(propertyType)) return

    rateGlobularClusterFormation=self%globularClusterFormationRateSpheroids_%rate(node)

    if (rateGlobularClusterFormation <= 0.0d0) return

    call spheroid%massStellarRate           (                                            & 
          &                                  -rateGlobularClusterFormation               &
          &                                 )

    call spheroid%floatRank0MetaPropertyRate(                                            &
          &                                   self%globularClusterStellarMassSpheroidID, &
          &                                  +rateGlobularClusterFormation               &
          &                                 )
    return
  end subroutine globularClusterFormationSpheroidsDifferentialEvolution
  
  subroutine globularClusterFormationSpheroidsGalaxiesMerge(self,node)
    !!{
    Combine integrals of star formation rate when galaxies merge.
    !!}
    use :: Error                           , only : Error_Report
    use :: Galacticus_Nodes                , only : nodeComponentDisk    , nodeComponentSpheroid
    use :: Satellite_Merging_Mass_Movements, only : destinationMergerDisk, destinationMergerSpheroid, destinationMergerUnmoved, enumerationDestinationMergerType
    implicit none
    class  (nodeOperatorGlobularClusterFormationSpheroids), intent(inout) :: self
    type   (treeNode                                     ), intent(inout) :: node
    type   (treeNode                                     ), pointer       :: nodeHost
    class  (nodeComponentDisk                            ), pointer       :: disk                   , diskHost
    class  (nodeComponentSpheroid                        ), pointer       :: spheroid               , spheroidHost
    type   (enumerationDestinationMergerType             )                :: destinationGasSatellite, destinationStarsSatellite, &
         &                                                                   destinationGasHost     , destinationStarsHost
    logical                                                               :: mergerIsMajor
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
  end subroutine globularClusterFormationSpheroidsGalaxiesMerge
  

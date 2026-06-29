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
  Implements a node operator class that computes the stellar mass-weighted ages of disk, spheroid and nuclear star cluster
  components.
  !!}
  
  use :: Star_Formation_Rates_Nuclear_Star_Clusters, only : starFormationRateNuclearStarClustersClass
  use :: Satellite_Merging_Mass_Movements          , only : mergerMassMovementsClass
  use :: Satellite_Merging_Nuclear_Star_Clusters   , only : nuclearStarClusterMovementsClass
  !![
  <nodeOperator name="nodeOperatornuclearStarClusterFormationTime">
    <description>
      A node operator class that computes the stellar mass-weighted ages of disk, spheroid and nuclear star cluster components. Intended to be paired
      with the \refClass{nodePropertyExtractornuclearStarClusterFormationTime} class to extract those ages for output.
    </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatornuclearStarClusterFormationTime
     !!{
     A node operator class that computes the stellar mass-weighted ages of disk and spheroid components.
     !!}
     private
     class           (mergerMassMovementsClass                 ), pointer :: mergerMassMovements_                  => null()
     class           (nuclearStarClusterMovementsClass         ), pointer :: nuclearStarClusterMovements_          => null()
     class           (starFormationRateNuclearStarClustersClass), pointer :: starFormationRateNuclearStarClusters_ => null()
     integer                                                              :: nuclearStarClusterFormationTimeID              , nuclearStarClusterFormationStarFormationRateID, &
         &                                                                   nuclearStarClusterFormationStellarMassID
     double precision                                                     :: massMinimum
   contains
     final     ::                          nuclearStarClusterFormationTimeDestructor
     procedure :: differentialEvolution => nuclearStarClusterFormationTimeDifferentialEvolution
     procedure :: galaxiesMerge         => nuclearStarClusterFormationTimeGalaxiesMerge
  end type nodeOperatornuclearStarClusterFormationTime
  
  interface nodeOperatornuclearStarClusterFormationTime
     !!{
     Constructors for the \refClass{nodeOperatornuclearStarClusterFormationTime} node operator class.
     !!}
     module procedure nuclearStarClusterFormationTimeConstructorParameters
     module procedure nuclearStarClusterFormationTimeConstructorInternal
  end interface nodeOperatornuclearStarClusterFormationTime
  
contains

  function nuclearStarClusterFormationTimeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatornuclearStarClusterFormationTime} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatornuclearStarClusterFormationTime)                :: self
    type (inputParameters                             ), intent(inout) :: parameters
    class(mergerMassMovementsClass                    ), pointer       :: mergerMassMovements_
    class(nuclearStarClusterMovementsClass            ), pointer       :: nuclearStarClusterMovements_  
    class(starFormationRateNuclearStarClustersClass   ), pointer       :: starFormationRateNuclearStarClusters_ => null()     
    double precision                                                   :: massMinimum
    !![
    <inputParameter>
    <name>massMinimum</name>
    <defaultValue>1.0d5</defaultValue>
    <description>Minimum mass of the nuclear star cluster to be considered as a nuclear star cluster</description>
    <source>parameters</source>
    </inputParameter>

    <objectBuilder class="mergerMassMovements"                  name="mergerMassMovements_"                  source="parameters"/>
    <objectBuilder class="nuclearStarClusterMovements"          name="nuclearStarClusterMovements_"          source="parameters"/>
    <objectBuilder class="starFormationRateNuclearStarClusters" name="starFormationRateNuclearStarClusters_" source="parameters"/>
    !!]
    self=nodeOperatornuclearStarClusterFormationTime(massMinimum,mergerMassMovements_,nuclearStarClusterMovements_,starFormationRateNuclearStarClusters_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="mergerMassMovements_"                 />
    <objectDestructor name="nuclearStarClusterMovements_"         />
    <objectDestructor name="starFormationRateNuclearStarClusters_"/>
    !!]
    return
  end function nuclearStarClusterFormationTimeConstructorParameters

  function nuclearStarClusterFormationTimeConstructorInternal(massMinimum,mergerMassMovements_,nuclearStarClusterMovements_,starFormationRateNuclearStarClusters_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatornuclearStarClusterFormationTime} node operator class.
    !!}
    implicit none
    type             (nodeOperatornuclearStarClusterFormationTime)                         :: self
    class            (mergerMassMovementsClass                    ), intent(in   ), target :: mergerMassMovements_
    class            (nuclearStarClusterMovementsClass            ), intent(in   ), target :: nuclearStarClusterMovements_
    class            (starFormationRateNuclearStarClustersClass   ), intent(in   ), target :: starFormationRateNuclearStarClusters_     
    double precision                                               , intent(in   )         :: massMinimum
    !![
    <constructorAssign variables="massMinimum, *mergerMassMovements_, *nuclearStarClusterMovements_, *starFormationRateNuclearStarClusters_"/>
    !!]
    
    !![
    <addMetaProperty component="NSC" name="nuclearStarClusterFormationTime"              id="self%nuclearStarClusterFormationTimeID"              isEvolvable="no" isCreator="yes"/>
    <addMetaProperty component="NSC" name="nuclearStarClusterFormationStarFormationRate" id="self%nuclearStarClusterFormationStarFormationRateID" isEvolvable="no" isCreator="yes"/>
    <addMetaProperty component="NSC" name="nuclearStarClusterFormationStellarMass"       id="self%nuclearStarClusterFormationStellarMassID"       isEvolvable="no" isCreator="yes"/>
    !!]
    return
  end function nuclearStarClusterFormationTimeConstructorInternal

  subroutine nuclearStarClusterFormationTimeDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatornuclearStarClusterFormationTime} node operator class.
    !!}
    implicit none
    type(nodeOperatornuclearStarClusterFormationTime), intent(inout) :: self

    !![
    <objectDestructor name="self%mergerMassMovements_"                 />
    <objectDestructor name="self%nuclearStarClusterMovements_"         />
    <objectDestructor name="self%starFormationRateNuclearStarClusters_"/>
    !!]
    return
  end subroutine nuclearStarClusterFormationTimeDestructor

  subroutine nuclearStarClusterFormationTimeDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Integrates unweighted and time-weighted star formation rates in nuclear star cluster components.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentNSC, nodeComponentBasic, propertyActive
    implicit none
    class           (nodeOperatornuclearStarClusterFormationTime), intent(inout), target  :: self
    type            (treeNode                                   ), intent(inout), target  :: node
    logical                                                      , intent(inout)          :: interrupt
    procedure       (interruptTask                              ), intent(inout), pointer :: functionInterrupt
    integer                                                      , intent(in   )          :: propertyType
    class           (nodeComponentBasic                         )               , pointer :: basic
    class           (nodeComponentNSC                           )               , pointer :: nuclearStarCluster
    double precision                                                                      :: stellarMassNuclearStarCluster  , time                               , &
        &                                                                                    formationTimeNuclearStarCluster, nuclearStarClusterStarFormationRate 
    logical                                                                               :: ageIsNonZero
    !$GLC attributes unused :: interrupt, functionInterrupt, propertyType

    ! Return immediately if active variables are requested.
    if (propertyActive(propertyType)) return
    ! Get the star formation rates.
    nuclearStarCluster => node%NSC()

    ageIsNonZero                       =.false.
    stellarMassNuclearStarCluster      =nuclearStarCluster%massStellar()
    nuclearStarClusterStarFormationRate=self%starFormationRateNuclearStarClusters_%rate(node)

    ! Find the current cosmic time.
    basic => node %basic()
    time  =  basic%time ()

    select type (nuclearStarCluster)
    type is (nodeComponentNSC)
       ! NSC does not yet exist - nothing to do here. class default
    class default
    if (nuclearStarClusterStarFormationRate> 0.0d0) then
        formationTimeNuclearStarCluster=nuclearStarCluster%floatRank0MetaPropertyGet(self%nuclearStarClusterFormationTimeID)-stellarMassNuclearStarCluster/nuclearStarClusterStarFormationRate
    else
        formationTimeNuclearStarCluster=nuclearStarCluster%floatRank0MetaPropertyGet(self%nuclearStarClusterFormationTimeID)
    end if
    if (formationTimeNuclearStarCluster > 0.0d0) ageIsNonZero=.true.

    if (.not.ageIsNonZero) then
        if (stellarMassNuclearStarCluster>= self%massMinimum) then 
            call nuclearStarCluster%floatRank0MetaPropertySet(self%nuclearStarClusterFormationTimeID             , time                               )
            call nuclearStarCluster%floatRank0MetaPropertySet(self%nuclearStarClusterFormationStarFormationRateID, nuclearStarClusterStarFormationRate)
            call nuclearStarCluster%floatRank0MetaPropertySet(self%nuclearStarClusterFormationStellarMassID      , stellarMassNuclearStarCluster      )
        else
            call nuclearStarCluster%floatRank0MetaPropertySet(self%nuclearStarClusterFormationTimeID             , 0.0d0)
            call nuclearStarCluster%floatRank0MetaPropertySet(self%nuclearStarClusterFormationStarFormationRateID, 0.0d0)
            call nuclearStarCluster%floatRank0MetaPropertySet(self%nuclearStarClusterFormationStellarMassID      , 0.0d0)
        end if
    end if 
    end select
    return
  end subroutine nuclearStarClusterFormationTimeDifferentialEvolution

  subroutine nuclearStarClusterFormationTimeGalaxiesMerge(self,node)
    !!{
    Combine integrals of star formation rate when galaxies merge.
    !!}
    use :: Error                           , only : Error_Report
    use :: Galacticus_Nodes                , only : nodeComponentDisk    , nodeComponentSpheroid    , nodeComponentNSC
    use :: Satellite_Merging_Mass_Movements, only : destinationMergerDisk, destinationMergerSpheroid, destinationMergerUnmoved, enumerationDestinationMergerType
    implicit none
    class  (nodeOperatornuclearStarClusterFormationTime), intent(inout) :: self
    type   (treeNode                                   ), intent(inout) :: node
    type   (treeNode                                   ), pointer       :: nodeHost
    class  (nodeComponentDisk                          ), pointer       :: disk                                      , diskHost
    class  (nodeComponentSpheroid                      ), pointer       :: spheroid                                  , spheroidHost
    class  (nodeComponentNSC                           ), pointer       :: nuclearStarClusterSatellite               , nuclearStarClusterHost
    type   (enumerationDestinationMergerType           )                :: destinationGasSatellite                   , destinationStarsSatellite                      , &
         &                                                                 destinationGasHost                        , destinationStarsHost
    logical                                                             :: mergerIsMajor                             , haveNuclearStarClusterSatellite                , &
         &                                                                 haveNuclearStarClusterHost                , nuclearStarClusterSatelliteIsDestroyed
    double precision                                                    :: timeFormationNuclearStarClusterHost       , timeFormationNuclearStarClusterSatellite       , &
         &                                                                 starFormationRateNuclearStarClusterHost   , starFormationRateNuclearStarClusterSatellite                 , &
         &                                                                 stellarMassFormationNuclearStarClusterHost, stellarMassFormationNuclearStarClusterSatellite

    ! Find the node to merge with.
    nodeHost                    => node    %mergesWith(                 )
    disk                        => node    %disk      (autoCreate=.true.)
    spheroid                    => node    %spheroid  (autoCreate=.true.)
    nuclearStarClusterSatellite => node    %NSC       (                 )
    diskHost                    => nodeHost%disk      (autoCreate=.true.)
    spheroidHost                => nodeHost%spheroid  (autoCreate=.true.)
    nuclearStarClusterHost      => nodeHost%NSC       (                 )

    select type (nuclearStarClusterSatellite)
    type is (nodeComponentNSC)
       haveNuclearStarClusterSatellite=.false.
    class default
       haveNuclearStarClusterSatellite=.true.
    end select 

    select type (nuclearStarClusterHost)
    type is (nodeComponentNSC)
       haveNuclearStarClusterHost=.false.
    class default
       haveNuclearStarClusterHost=.true.
    end select

    ! Get mass movement descriptors.
    call self%mergerMassMovements_        %get        (node,destinationGasSatellite,destinationStarsSatellite,destinationGasHost,destinationStarsHost,mergerIsMajor)
    ! This descriptor return if the nuclear star cluster in the satellite is destroyed.
    call self%nuclearStarClusterMovements_%isDestroyed(node,nuclearStarClusterSatelliteIsDestroyed                                                                 )

    if (haveNuclearStarClusterHost     ) then
       timeFormationNuclearStarClusterHost       =nuclearStarClusterHost%floatRank0MetaPropertyGet(self%nuclearStarClusterFormationTimeID             )
       starFormationRateNuclearStarClusterHost   =nuclearStarClusterHost%floatRank0MetaPropertyGet(self%nuclearStarClusterFormationStarFormationRateID)
       stellarMassFormationNuclearStarClusterHost=nuclearStarClusterHost%floatRank0MetaPropertyGet(self%nuclearStarClusterFormationStellarMassID      )
    else
       timeFormationNuclearStarClusterHost       =0.0d0
       starFormationRateNuclearStarClusterHost   =0.0d0
       stellarMassFormationNuclearStarClusterHost=0.0d0
    end if

    if (haveNuclearStarClusterSatellite) then
       timeFormationNuclearStarClusterSatellite       =nuclearStarClusterSatellite%floatRank0MetaPropertyGet(self%nuclearStarClusterFormationTimeID             )
       starFormationRateNuclearStarClusterSatellite   =nuclearStarClusterSatellite%floatRank0MetaPropertyGet(self%nuclearStarClusterFormationStarFormationRateID)
       stellarMassFormationNuclearStarClusterSatellite=nuclearStarClusterSatellite%floatRank0MetaPropertyGet(self%nuclearStarClusterFormationStellarMassID      )
    else
       timeFormationNuclearStarClusterSatellite       =0.0d0
       starFormationRateNuclearStarClusterSatellite   =0.0d0
       stellarMassFormationNuclearStarClusterSatellite=0.0d0
    end if

    select case (destinationStarsHost%ID)
    ! If the host contains a nuclear star cluster, it should remains as it is. If it does not, do nothing.
    case (destinationMergerDisk    %ID) 
        if (haveNuclearStarClusterHost) then  
            call nuclearStarClusterHost%floatRank0MetaPropertySet(self%nuclearStarClusterFormationTimeID             , timeFormationNuclearStarClusterHost       )
            call nuclearStarClusterHost%floatRank0MetaPropertySet(self%nuclearStarClusterFormationStarFormationRateID, starFormationRateNuclearStarClusterHost   )
            call nuclearStarClusterHost%floatRank0MetaPropertySet(self%nuclearStarClusterFormationStellarMassID      , stellarMassFormationNuclearStarClusterHost)
        end if
    case (destinationMergerSpheroid%ID)
        if (haveNuclearStarClusterHost) then  
            call nuclearStarClusterHost%floatRank0MetaPropertySet(self%nuclearStarClusterFormationTimeID             , timeFormationNuclearStarClusterHost       )
            call nuclearStarClusterHost%floatRank0MetaPropertySet(self%nuclearStarClusterFormationStarFormationRateID, starFormationRateNuclearStarClusterHost   )
            call nuclearStarClusterHost%floatRank0MetaPropertySet(self%nuclearStarClusterFormationStellarMassID      , stellarMassFormationNuclearStarClusterHost)
        end if
    case (destinationMergerUnmoved%ID)
     ! Do nothing.
    case default
        call Error_Report('unrecognized movesTo descriptor'//{introspection:location})
    end select

    ! Now, we should check what happens. If the nuclear star cluster in the satellite is destroyed, then we must set the time formation as zero.
    select case (destinationStarsSatellite%ID)
    ! The nuclear star cluster is not destroyed. 
    case (destinationMergerDisk    %ID)
    if (.not.nuclearStarClusterSatelliteIsDestroyed) then
        ! Auto create a nuclear star cluster in the host if needed.
        if (.not.haveNuclearStarClusterHost) nuclearStarClusterHost => nodeHost%NSC(autoCreate=.true.) 
        call nuclearStarClusterHost%floatRank0MetaPropertySet(self%nuclearStarClusterFormationTimeID             , timeFormationNuclearStarClusterHost             &
                &                                                                                                 +timeFormationNuclearStarClusterSatellite        &
                &                                            )
        call nuclearStarClusterHost%floatRank0MetaPropertySet(self%nuclearStarClusterFormationStarFormationRateID, starFormationRateNuclearStarClusterHost         &
                &                                                                                                 +starFormationRateNuclearStarClusterSatellite    &
                &                                            )
        call nuclearStarClusterHost%floatRank0MetaPropertySet(self%nuclearStarClusterFormationStellarMassID      , stellarMassFormationNuclearStarClusterHost      &
                &                                                                                                 +stellarMassFormationNuclearStarClusterSatellite &
                &                                            )
    end if
    case (destinationMergerSpheroid%ID)
    if (.not.nuclearStarClusterSatelliteIsDestroyed) then
        ! Auto create a nuclear star cluster in the host if needed.
        if (.not.haveNuclearStarClusterHost) nuclearStarClusterHost => nodeHost%NSC(autoCreate=.true.) 
        call nuclearStarClusterHost%floatRank0MetaPropertySet(self%nuclearStarClusterFormationTimeID             , timeFormationNuclearStarClusterHost             &
                &                                                                                                 +timeFormationNuclearStarClusterSatellite        &
                &                                            )
        call nuclearStarClusterHost%floatRank0MetaPropertySet(self%nuclearStarClusterFormationStarFormationRateID, starFormationRateNuclearStarClusterHost         &
                &                                                                                                 +starFormationRateNuclearStarClusterSatellite    &
                &                                            )
        call nuclearStarClusterHost%floatRank0MetaPropertySet(self%nuclearStarClusterFormationStellarMassID      , stellarMassFormationNuclearStarClusterHost      &
                &                                                                                                 +stellarMassFormationNuclearStarClusterSatellite &
                &                                            )

    end if
    case default
        call Error_Report('unrecognized movesTo descriptor'//{introspection:location})
    end select
    
    if (haveNuclearStarClusterSatellite) then 
        call nuclearStarClusterSatellite%floatRank0MetaPropertySet(self%nuclearStarClusterFormationTimeID             ,+0.0d0)
        call nuclearStarClusterSatellite%floatRank0MetaPropertySet(self%nuclearStarClusterFormationStarFormationRateID,+0.0d0)
        call nuclearStarClusterSatellite%floatRank0MetaPropertySet(self%nuclearStarClusterFormationStellarMassID      ,+0.0d0)
    end if 
    return
  end subroutine nuclearStarClusterFormationTimeGalaxiesMerge
  
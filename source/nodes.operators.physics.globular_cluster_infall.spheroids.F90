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
  Implements a node operator class that performs star formation in spheroids.
  !!}

  use :: Globular_Cluster_Infall_Rates_Spheroids, only : globularClusterInfallRateSpheroidsClass

  !![
  <nodeOperator name="nodeOperatorglobularClusterInfallSpheroids">
   <description>A node operator class that performs star formation.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorglobularClusterInfallSpheroids
     !!{
     A node operator class that performs star formation.
     !!}
     private
     class  (globularClusterInfallRateSpheroidsClass), pointer :: globularClusterInfallRateSpheroids_ => null()
     integer                                                   :: globularClusterStellarMassSpheroidID
     double precision                                          :: fraction 
   contains
     final     ::                          globularClusterInfallSpheroidsDestructor
     procedure :: differentialEvolution => globularClusterInfallSpheroidDifferentialEvolution
  end type nodeOperatorglobularClusterInfallSpheroids
  
  interface nodeOperatorglobularClusterInfallSpheroids
     !!{
     Constructors for the {\normalfont \ttfamily globularClusterFormationSpheroids} node operator class.
     !!}
     module procedure globularClusterInfallSpheroidsConstructorParameters
     module procedure globularClusterInfallSpheroidsConstructorInternal
  end interface nodeOperatorglobularClusterInfallSpheroids
  
contains

  function globularClusterInfallSpheroidsConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily globularClusterFormation} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorglobularClusterInfallSpheroids)                :: self
    type (inputParameters                           ), intent(inout) :: parameters
    class(globularClusterInfallRateSpheroidsClass   ), pointer       :: globularClusterInfallRateSpheroids_
    double precision                                                 :: fraction

    !![
    <inputParameter>
      <name>fraction</name>
      <defaultValue>0.5d0</defaultValue>
      <source>parameters</source>
      <description>Fraction of the the globular clusters that survives to merge with the nuclear star cluster.</description>
    </inputParameter>
    <objectBuilder class="globularClusterInfallRateSpheroids" name="globularClusterInfallRateSpheroids_" source="parameters"/>
    !!]
    self=nodeOperatorglobularClusterInfallSpheroids(fraction,globularClusterInfallRateSpheroids_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="globularClusterInfallRateSpheroids_"/>
    !!]
    return
  end function globularClusterInfallSpheroidsConstructorParameters

  function globularClusterInfallSpheroidsConstructorInternal(fraction,globularClusterInfallRateSpheroids_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily globularClusterFormationDisks} node operator class.
    !!}
    implicit none
    type (nodeOperatorglobularClusterInfallSpheroids)                        :: self
    class(globularClusterInfallRateSpheroidsClass   ), intent(in   ), target :: globularClusterInfallRateSpheroids_
    double precision                                 , intent(in   )         :: fraction
    !![
    <constructorAssign variables="fraction,*globularClusterInfallRateSpheroids_"/>
    !!]
    !![
    <addMetaProperty component="spheroid" name="globularClusterStellarMassSpheroid" id="self%globularClusterStellarMassSpheroidID" isEvolvable="yes" isCreator="no"/>
    !!]
    return
  end function globularClusterInfallSpheroidsConstructorInternal

  subroutine globularClusterInfallSpheroidsDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily globularClusterFormationDisks} node operator class.
    !!}
    implicit none
    type(nodeOperatorglobularClusterInfallSpheroids), intent(inout) :: self

    !![
    <objectDestructor name="self%globularClusterInfallRateSpheroids_"/>
    !!]
    return
  end subroutine globularClusterInfallSpheroidsDestructor
  
  subroutine globularClusterInfallSpheroidDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Perform globular cluster dissolution in an spheroid.
    !!}
    use :: Galacticus_Nodes, only : propertyInactive, propertyTypeActive      , propertyEvaluate, nodeComponentSpheroid, &
            &                       nodeComponentNSC, nodeComponentNSCStandard
    implicit none
    class           (nodeOperatorglobularClusterInfallSpheroids), intent(inout), target  :: self
    type            (treeNode                                  ), intent(inout), target  :: node
    logical                                                     , intent(inout)          :: interrupt
    procedure       (interruptTask                             ), intent(inout), pointer :: functionInterrupt
    integer                                                     , intent(in   )          :: propertyType
    class           (nodeComponentSpheroid                     )               , pointer :: spheroid
    class           (nodeComponentNSC                          )               , pointer :: nuclearStarCluster
    double precision                                                                     :: rateGlobularClusterInfall
    ! Check for a realistic disk, return immediately if spheroid is unphysical.
    spheroid => node%spheroid()
    if     (         spheroid%angularMomentum() <= 0.0d0 &
         &       .or.spheroid%radius         () <= 0.0d0 &
         &       .or.spheroid%massGas        () <= 0.0d0 &
         &       .or.spheroid%massStellar    () <= 0.0d0 &
         & ) return
    if (propertyInactive(propertyType)) return

    rateGlobularClusterInfall=self%fraction*self%globularClusterInfallRateSpheroids_%rate(node)

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
       call spheroid          %floatRank0MetaPropertyRate(                                          &
          &                                              self%globularClusterStellarMassSpheroidID, &
          &                                                  -rateGlobularClusterInfall             &
          &                                             )
       ! WARNING, here we need to adjust stellar histories in both, spheroid and nuclear star cluster components.
       call nuclearStarCluster%massStellarRate(+rateGlobularClusterInfall)
    end select
    return
  end subroutine globularClusterInfallSpheroidDifferentialEvolution

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

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
  Implements a node operator class that performs globular cluster evaporation in spheroids.
  !!}

  use :: Globular_Cluster_Evaporation_Rates_Spheroids, only : globularClusterEvaporationRateSpheroidsClass

  !![
  <nodeOperator name="nodeOperatorGlobularClusterEvaporationSpheroids">
   <description>A node operator class that performs star formation.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorGlobularClusterEvaporationSpheroids
     !!{
     A node operator class that performs star formation.
     !!}
     private
     class  (globularClusterEvaporationRateSpheroidsClass), pointer :: globularClusterEvaporationRateSpheroids_ => null()
     integer                                                        :: globularClusterStellarMassSpheroidID 
   contains
     final     ::                                   globularClusterEvaporationSpheroidsDestructor
     procedure :: differentialEvolution          => globularClusterEvaporationSpheroidsDifferentialEvolution
  end type nodeOperatorGlobularClusterEvaporationSpheroids
  
  interface nodeOperatorGlobularClusterEvaporationSpheroids
     !!{
     Constructors for the {\normalfont \ttfamily globularClusterFormationSpheroids} node operator class.
     !!}
     module procedure globularClusterEvaporationSpheroidsConstructorParameters
     module procedure globularClusterEvaporationSpheroidsConstructorInternal
  end interface nodeOperatorGlobularClusterEvaporationSpheroids
  
contains

  function globularClusterEvaporationSpheroidsConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily globularClusterFormation} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorGlobularClusterEvaporationSpheroids)                :: self
    type (inputParameters                                ), intent(inout) :: parameters
    class(globularClusterEvaporationRateSpheroidsClass   ), pointer       :: globularClusterEvaporationRateSpheroids_

    !![
    <objectBuilder class="globularClusterEvaporationRateSpheroids" name="globularClusterEvaporationRateSpheroids_" source="parameters"/>
    !!]
    self=nodeOperatorGlobularClusterEvaporationSpheroids(globularClusterEvaporationRateSpheroids_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="globularClusterEvaporationRateSpheroids_"/>
    !!]
    return
  end function globularClusterEvaporationSpheroidsConstructorParameters

  function globularClusterEvaporationSpheroidsConstructorInternal(globularClusterEvaporationRateSpheroids_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily globularClusterFormationSpheroids} node operator class.
    !!}
    implicit none
    type (nodeOperatorGlobularClusterEvaporationSpheroids)                        :: self
    class(globularClusterEvaporationRateSpheroidsClass   ), intent(in   ), target :: globularClusterEvaporationRateSpheroids_

    !![
    <constructorAssign variables="*globularClusterEvaporationRateSpheroids_"/>
    !!]
    !![
    <addMetaProperty component="spheroid" name="globularClusterStellarMassSpheroid" id="self%globularClusterStellarMassSpheroidID" isEvolvable="yes" isCreator="no"/>
    !!]
    return
  end function globularClusterEvaporationSpheroidsConstructorInternal

  subroutine globularClusterEvaporationSpheroidsDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily globularClusterFormationSpheroids} node operator class.
    !!}
    implicit none
    type(nodeOperatorGlobularClusterEvaporationSpheroids), intent(inout) :: self

    !![
    <objectDestructor name="self%globularClusterEvaporationRateSpheroids_"/>
    !!]
    return
  end subroutine globularClusterEvaporationSpheroidsDestructor
  
  subroutine globularClusterEvaporationSpheroidsDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Perform globular cluster evaporation  in an spheroid.
    !!}
    use :: Galacticus_Nodes, only : propertyInactive, propertyTypeActive, propertyEvaluate, nodeComponentSpheroid
    implicit none
    class           (nodeOperatorGlobularClusterEvaporationSpheroids), intent(inout), target  :: self
    type            (treeNode                                       ), intent(inout), target  :: node
    logical                                                          , intent(inout)          :: interrupt
    procedure       (interruptTask                                  ), intent(inout), pointer :: functionInterrupt
    integer                                                          , intent(in   )          :: propertyType
    class           (nodeComponentSpheroid                          )               , pointer :: spheroid
    double precision                                                                          :: rateGlobularClusterEvaporation
    ! Check for a realistic spheroid, return immediately if spheroid is unphysical.
    spheroid => node%spheroid()
    if     (         spheroid%angularMomentum() <= 0.0d0 &
         &       .or.spheroid%radius         () <= 0.0d0 &
         &       .or.spheroid%massGas        () <= 0.0d0 &
         &       .or.spheroid%massStellar    () <= 0.0d0 &
         & ) return
    if (propertyInactive(propertyType)) return

    rateGlobularClusterEvaporation=self%globularClusterEvaporationRateSpheroids_%rate(node)

    if (rateGlobularClusterEvaporation <= 0.0d0) return

    call spheroid%massStellarRate           (+rateGlobularClusterEvaporation)

    call spheroid%floatRank0MetaPropertyRate(                                            &
          &                                   self%globularClusterStellarMassSpheroidID, &
          &                                  -rateGlobularClusterEvaporation             &
          &                                 )
    return
  end subroutine globularClusterEvaporationSpheroidsDifferentialEvolution

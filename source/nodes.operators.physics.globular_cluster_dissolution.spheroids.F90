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

  use :: Globular_Cluster_Dissolution_Rates_Spheroids, only : globularClusterDissolutionRateSpheroidsClass

  !![
  <nodeOperator name="nodeOperatorGlobularClusterDissolutionSpheroids">
   <description>A node operator class that performs star formation.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorGlobularClusterDissolutionSpheroids
     !!{
     A node operator class that performs star formation.
     !!}
     private
     class  (globularClusterDissolutionRateSpheroidsClass), pointer :: globularClusterDissolutionRateSpheroids_ => null()
     integer                                                        :: globularClusterStellarMassSpheroidID 
   contains
     final     ::                                   globularClusterDissolutionSpheroidsDestructor
     procedure :: differentialEvolution          => globularClusterDissolutionSpheroidsDifferentialEvolution
  end type nodeOperatorGlobularClusterDissolutionSpheroids
  
  interface nodeOperatorGlobularClusterDissolutionSpheroids
     !!{
     Constructors for the {\normalfont \ttfamily globularClusterFormationSpheroids} node operator class.
     !!}
     module procedure globularClusterDissolutionSpheroidsConstructorParameters
     module procedure globularClusterDissolutionSpheroidsConstructorInternal
  end interface nodeOperatorGlobularClusterDissolutionSpheroids
  
contains

  function globularClusterDissolutionSpheroidsConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily globularClusterFormation} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorGlobularClusterDissolutionSpheroids)                :: self
    type (inputParameters                                ), intent(inout) :: parameters
    class(globularClusterDissolutionRateSpheroidsClass   ), pointer       :: globularClusterDissolutionRateSpheroids_

    !![
    <objectBuilder class="globularClusterDissolutionRateSpheroids" name="globularClusterDissolutionRateSpheroids_" source="parameters"/>
    !!]
    self=nodeOperatorGlobularClusterDissolutionSpheroids(globularClusterDissolutionRateSpheroids_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="globularClusterDissolutionRateSpheroids_"/>
    !!]
    return
  end function globularClusterDissolutionSpheroidsConstructorParameters

  function globularClusterDissolutionSpheroidsConstructorInternal(globularClusterDissolutionRateSpheroids_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily globularClusterFormationSpheroids} node operator class.
    !!}
    implicit none
    type (nodeOperatorGlobularClusterDissolutionSpheroids)                        :: self
    class(globularClusterDissolutionRateSpheroidsClass   ), intent(in   ), target :: globularClusterDissolutionRateSpheroids_

    !![
    <constructorAssign variables="*globularClusterDissolutionRateSpheroids_"/>
    !!]
    !![
    <addMetaProperty component="spheroid" name="globularClusterStellarMassSpheroid" id="self%globularClusterStellarMassSpheroidID" isEvolvable="yes" isCreator="no"/>
    !!]
    return
  end function globularClusterDissolutionSpheroidsConstructorInternal

  subroutine globularClusterDissolutionSpheroidsDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily globularClusterFormationSpheroids} node operator class.
    !!}
    implicit none
    type(nodeOperatorGlobularClusterDissolutionSpheroids), intent(inout) :: self

    !![
    <objectDestructor name="self%globularClusterDissolutionRateSpheroids_"/>
    !!]
    return
  end subroutine globularClusterDissolutionSpheroidsDestructor
  
  subroutine globularClusterDissolutionSpheroidsDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Perform globular cluster dissolution in an spheroid.
    !!}
    use :: Galacticus_Nodes, only : propertyInactive, propertyTypeActive, propertyEvaluate, nodeComponentSpheroid
    implicit none
    class           (nodeOperatorGlobularClusterDissolutionSpheroids), intent(inout), target  :: self
    type            (treeNode                                       ), intent(inout), target  :: node
    logical                                                          , intent(inout)          :: interrupt
    procedure       (interruptTask                                  ), intent(inout), pointer :: functionInterrupt
    integer                                                          , intent(in   )          :: propertyType
    class           (nodeComponentSpheroid                          )               , pointer :: spheroid
    double precision                                                                          :: globularClusterMass, rateGlobularClusterDissolution
    ! Check for a realistic spheroid, return immediately if spheroid is unphysical.
    spheroid => node%spheroid()
    if     (         spheroid%angularMomentum() <= 0.0d0 &
         &       .or.spheroid%radius         () <= 0.0d0 &
         &       .or.spheroid%massGas        () <= 0.0d0 &
         &       .or.spheroid%massStellar    () <= 0.0d0 &
         & ) return
    if (propertyInactive(propertyType)) return

    globularClusterMass           = spheroid%floatRank0MetaPropertyGet(self%globularClusterStellarMassSpheroidID)
    rateGlobularClusterDissolution= self%globularClusterDissolutionRateSpheroids_%rate(node)

    if (rateGlobularClusterDissolution<=0.0d0.or.globularClusterMass<=0.0d0) return

    call spheroid%massStellarRate           (+rateGlobularClusterDissolution)

    call spheroid%floatRank0MetaPropertyRate(                                            &
          &                                   self%globularClusterStellarMassSpheroidID, &
          &                                  -rateGlobularClusterDissolution             &
          &                                 )
    return
  end subroutine globularClusterDissolutionSpheroidsDifferentialEvolution

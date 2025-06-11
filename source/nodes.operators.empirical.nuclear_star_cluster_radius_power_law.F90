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
  Implements a node operator class that implements an an empirical power law relationship between disk stellar mass and
  stellar radius.
  !!}

  !![
  <nodeOperator name="nodeOperatorNuclearStarClusterRadiusPowerLaw">
   <description>
    A node operator that implements an an empirical power law relationship between disk stellar radius and stellar mass.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorNuclearStarClusterRadiusPowerLaw
     !!{
     Implements a power law prescription for the stellar mass--stellar radius relation of disks. Specifically:
     \begin{equation}
       r_\mathrm{s} = \radiusPivot \left( \frac{M_\star+M_\mathrm{gas}}{\massPivot} \right)^\radiusPivot, 
     \end{equation}
     where $r_\mathrm{s}$ is the nuclear star cluster scale radius, $M_\star$ is the stellar mass of the nuclear star cluster, 
     $M_\mathrm{gas}$ is the gaseous mass of the nuclear star cluster, and $\radiusPivot$, and $\massPivot$,
     are free parameters.
     !!}
     private
     double precision :: radiusPivot, massPivot, &
      &                  exponent
   contains
     !![
     <methods>
       <method method="update" description="Update the nuclear star cluster radius to be consistent with its stellar mass."/>
     </methods>
     !!]
     procedure :: update                              => nuclearStarClusterRadiusPowerLawUpdate
     procedure :: differentialEvolutionSolveAnalytics => nuclearStarClusterRadiusPowerLawSolveAnalytics
     procedure :: nodesMerge                          => nuclearStarClusterRadiusPowerLawNodesMerge
  end type nodeOperatorNuclearStarClusterRadiusPowerLaw
  
  interface nodeOperatorNuclearStarClusterRadiusPowerLaw
     !!{
     Constructors for the \refClass{nodeOperatorNuclearStarClusterRadiusPowerLaw} node operator class.
     !!}
     module procedure nuclearStarClusterRadiusPowerLawConstructorParameters
     module procedure nuclearStarClusterRadiusPowerLawConstructorInternal
  end interface nodeOperatorNuclearStarClusterRadiusPowerLaw
  
contains

  function nuclearStarClusterRadiusPowerLawConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorNuclearStarClusterRadiusPowerLaw} class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorNuclearStarClusterRadiusPowerLaw)                :: self
    type            (inputParameters                             ), intent(inout) :: parameters
    double precision                                                              :: radiusPivot, massPivot, &
      &                                                                              exponent

    !![
    <inputParameter>
      <name>radiusPivot</name>
      <source>parameters</source>
      <description>Exponent $\radiusPivot$ in the power law fit.</description> 
      <defaultValue>3.3d-6</defaultValue>
      <defaultSource>\cite{Neumayer2020}</defaultSource>
    </inputParameter>
    <inputParameter>
      <name>massPivot</name>
      <source>parameters</source>
      <description>The $\massPivot$ in the power law fit.</description>
      <defaultValue>1.0d6</defaultValue>
      <defaultSource>\cite{Neumayer2020}</defaultSource>
    </inputParameter>
    <inputParameter>
      <name>exponent</name>
      <source>parameters</source>
      <description>Exponent in the power law fit.</description>
      <defaultValue>0.5d0</defaultValue>
      <defaultSource>\cite{Neumayer2020}</defaultSource>
    </inputParameter>
    !!]
    self=nuclearStarClusterRadiusPowerLawConstructorInternal(radiusPivot,massPivot,exponent)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
  end function nuclearStarClusterRadiusPowerLawConstructorParameters

  function nuclearStarClusterRadiusPowerLawConstructorInternal(radiusPivot,massPivot,exponent) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorNuclearStarClusterRadiusPowerLaw} node operator class.
    !!}
    implicit none
    type            (nodeOperatorNuclearStarClusterRadiusPowerLaw)             :: self
    double precision                                              , intent(in) :: radiusPivot, massPivot, &
      &                                                                           exponent 
    !![
    <constructorAssign variables="radiusPivot, massPivot, exponent"/>
    !!]
    return
  end function nuclearStarClusterRadiusPowerLawConstructorInternal

  subroutine nuclearStarClusterRadiusPowerLawUpdate(self,node)
    !!{
    Update radius of the nuclear star cluster.
    !!} 
    use :: Galacticus_Nodes, only : nodeComponentNSC, nodeComponentNSCStandard
    implicit none
    class           (nodeOperatorNuclearStarClusterRadiusPowerLaw), intent(inout) :: self
    type            (treeNode                                    ), intent(inout) :: node
    class           (nodeComponentNSC                            ), pointer       :: nuclearStarCluster
    double precision                                                              :: radius

    nuclearStarCluster=> node%NSC()

    select type(nuclearStarCluster)

    type is (nodeComponentNSC)
      return

    class default
    radius            = +self%radiusPivot                    &
         &              *(                                   &
                          ( nuclearStarCluster%massStellar() &
         &                 +nuclearStarCluster%massGas    () &
         &                )/self%massPivot                   &
         &               )**self%exponent
    call nuclearStarCluster%radiusSet(radius)
    end select
    return  
  end subroutine nuclearStarClusterRadiusPowerLawUpdate
   
  subroutine nuclearStarClusterRadiusPowerLawSolveAnalytics(self,node,time)
    !!{
    Set the radius of the nuclear star cluster.
    !!}
    implicit none
    class           (nodeOperatorNuclearStarClusterRadiusPowerLaw), intent(inout) :: self
    type            (treeNode                                    ), intent(inout) :: node
    double precision                                              , intent(in   ) :: time
    !$GLC attributes unused :: time

    call self%update(node)
    return
  end subroutine nuclearStarClusterRadiusPowerLawSolveAnalytics

  subroutine nuclearStarClusterRadiusPowerLawNodesMerge(self,node)
    !!{
    Update the radius of the nuclear star cluster after a merger.
    !!}
    implicit none
    class(nodeOperatorNuclearStarClusterRadiusPowerLaw), intent(inout) :: self
    type (treeNode                                    ), intent(inout) :: node

    call self%update(node)
    return
  end subroutine nuclearStarClusterRadiusPowerLawNodesMerge

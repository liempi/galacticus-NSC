!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024
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
  Implements a node operator class that inserts an empirical model of the radius which scales with the radius of the NSC.
  !!}

  !![
  <nodeOperator name="nodeOperatorDarkCoreRadius">
   <description>A node operator class that inserts an empirical model of the formation history of a massive elliptical galaxy.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorDarkCoreRadius
     !!{     
     A node operator class that inserts an empirical model for the evolution of the dark core radius. Specifically, this node operator
     assumes that the radius scales with the radius of the host nuclear star cluster
     \begin{equation}
     r_{\rm DC} = \epsilon r_{\rm NSC}
     \end{equation}
     where $\epsilon=${\normalfont \ttfamily [efficiency]} is an efficiency free parameter.
     !!}
     private
     integer          :: darkCoreRadiusID      
     double precision :: efficiency
   contains
     !![
     <methods>
       <method method="update" description="Update the dark core radius to be consistent with its host nuclear star cluster."/>
     </methods>
     !!]
     procedure :: update                              => darkCoreRadiusUpdate
     procedure :: differentialEvolutionInactives      => darkCoreRadiusDifferentialEvolutionInactives
     procedure :: differentialEvolutionSolveAnalytics => darkCoreRadiusDifferentialEvolutionSolveAnalytics
     procedure :: nodesMerge                          => darkCoreRadiusNodesMerge
  end type nodeOperatorDarkCoreRadius
  
  interface nodeOperatorDarkCoreRadius
     !!{
     Constructors for the {\normalfont \ttfamily darkCoreRadius} node operator class.
     !!}
     module procedure darkCoreRadiusConstructorParameters
     module procedure darkCoreRadiusConstructorInternal
  end interface nodeOperatorDarkCoreRadius
  
contains

  function darkCoreRadiusConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily starFormation} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorDarkCoreRadius)                :: self
    type            (inputParameters           ), intent(inout) :: parameters
    double precision                                            :: efficiency

    !![
    <inputParameter>
      <name>efficiency</name>
      <source>parameters</source>
      <defaultValue>0.1d0</defaultValue>
      <description>The assumed value with dark Core radius scales with the nuclear star cluster radius.</description>
    </inputParameter>
    !!]
    self=nodeOperatorDarkCoreRadius(efficiency)

    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function darkCoreRadiusConstructorParameters

  function darkCoreRadiusConstructorInternal(efficiency) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily darkCoreRadius} node operator class.
    !!}
    implicit none
    type            (nodeOperatorDarkCoreRadius)                :: self
    double precision                            , intent(in   ) :: efficiency
    !![
    <constructorAssign variables="efficiency"/>
    <addMetaProperty component="NSC" name="darkCoreRadius" id="self%darkCoreRadiusID" isEvolvable="no" isCreator="yes"/>
    !!]
    return
  end function darkCoreRadiusConstructorInternal
  
  subroutine darkCoreRadiusUpdate(self,node)
    !!{
    Update radius of the disk.
    !!} 
    use :: Galacticus_Nodes, only : nodeComponentNSC
    implicit none
    class           (nodeOperatorDarkCoreRadius), intent(inout) :: self
    type            (treeNode                  ), intent(inout) :: node
    class           (nodeComponentNSC          ), pointer       :: nuclearStarCluster
    double precision                                         :: nuclearStarClusterRadius      , darkCoreRadius, &
      &                                                         nuclearStarClusterDarkCoreMass
    nuclearStarCluster             => node              %   NSC()
    nuclearStarClusterRadius       =  nuclearStarCluster%radius()
    nuclearStarClusterDarkCoreMass =  nuclearStarCluster%massDarkCore()

    if (nuclearStarClusterRadius> 0.0d0.and.nuclearStarClusterDarkCoreMass>0.0d0) then
      darkCoreRadius = self%efficiency*nuclearStarClusterRadius
      call nuclearStarCluster%floatRank0MetaPropertySet(self%darkCoreRadiusID, darkCoreRadius)
    end if
    return
  end subroutine darkCoreRadiusUpdate

  subroutine darkCoreRadiusDifferentialEvolutionInactives(self,node)
    !!{
    Mark radius as inactive for ODE solving.    
    !!}
    use :: Galacticus_Nodes, only : nodeComponentNSC
    implicit none
    class(nodeOperatorDarkCoreRadius), intent(inout) :: self
    type (treeNode                  ), intent(inout) :: node
    class(nodeComponentNSC          ), pointer       :: nuclearStarCluster
    
    ! Get nuclear star cluster component.
    nuclearStarCluster => node%NSC     ()
    ! Mark as inactive.
    select type (nuclearStarCluster)
    type is (nodeComponentNSC     )
        ! Nuclear star cluster does not yet exist - nothing to do here.
    class default
       call nuclearStarCluster%floatRank0MetaPropertyInactive(self%darkCoreRadiusID)
    end select
    return
  end subroutine darkCoreRadiusDifferentialEvolutionInactives
  
  subroutine darkCoreRadiusDifferentialEvolutionSolveAnalytics(self,node,time)
      !!{
        Set the radius of the dark core.
      !!}
    implicit none
    class(nodeOperatorDarkCoreRadius), intent(inout) :: self
    type (treeNode                  ), intent(inout) :: node
    double precision                 , intent(in   ) :: time
    !$GLC attributes unused :: time
    
    call self%update(node)
    return 
  end subroutine darkCoreRadiusDifferentialEvolutionSolveAnalytics

    subroutine darkCoreRadiusNodesMerge(self,node)
    !!{
    Update the radius of the dark core after a merger.
    !!}
    implicit none
    class(nodeOperatorDarkCoreRadius), intent(inout) :: self
    type (treeNode                  ), intent(inout) :: node

    call self%update(node)
    return
  end subroutine darkCoreRadiusNodesMerge

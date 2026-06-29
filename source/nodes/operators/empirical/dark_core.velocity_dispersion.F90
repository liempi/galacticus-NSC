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
  Implements a node operator class that computes the dark core velocity dispersion assuming virial
  equilibrium and a contraction as result of gas inflows.
  !!}

  !![
  <nodeOperator name="nodeOperatorDarkCoreVelocityDispersion">
   <description>
    A node operator that implements an analytical relationship between dark core mass and its velocity dispersion.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorDarkCoreVelocityDispersion
     !!{
      Computes the velocity dispersion of the dark core assuming virial equilibrum and a contraction due to gas inflows.
       \begin{equation}
       \sigma_\mathrm{dark\, core} =  \sqrt{ \frac{2GM_\mathrm{dark\, core}(M_\mathrm{dark\, core}+M_\mathrm{gas}^{\rm NSC} (<r_\mathrm{dark\,core}))}{k_\mathrm{contraction} r_\mathrm{dark\, core}M_\mathrm{dark\, core}}  }
     \end{equation}
     where $\sigma_\mathrm{dark\, core}$ is the velocity dispersion of the dark core, $M_\mathrm{dark\, core}$ is the mass of the dark core,
     $M_\mathrm{gas}^{\rm NSC} (<r_\mathrm{dark\,core})$ is the enclosed nuclear star cluster gas, and $k_\mathrm{contraction}${\normalfont \ttfamily [contractionFactor]} a free parameter that mimics the contraction as 
     response to the gas inflow.
     !!}
     private
     double precision :: contractionFactor
     integer          :: darkCoreRadiusID            , darkCoreGasMassID, &
       &                 darkCoreVelocityDispersionID
   contains
     !![
     <methods>
       <method method="update" description="Update the velocity dispersion of the dark core."/>
     </methods>
     !!]
     procedure :: update                              => darkCoreVelocityDispersionUpdate
     procedure :: nodeInitialize                      => darkCoreVelocityDispersionNodeInitialize
     procedure :: differentialEvolutionSolveAnalytics => darkCoreVelocityDispersionSolveAnalytics
     procedure :: nodesMerge                          => darkCoreVelocityDispersionNodesMerge
  end type nodeOperatorDarkCoreVelocityDispersion
  
  interface nodeOperatorDarkCoreVelocityDispersion
     !!{
     Constructors for the \refClass{nodeOperatorDarkCoreVelocityDispersion} node operator class.
     !!}
     module procedure darkCoreVelocityDispersionConstructorParameters
     module procedure darkCoreVelocityDispersionConstructorInternal
  end interface nodeOperatorDarkCoreVelocityDispersion
  
contains

  function darkCoreVelocityDispersionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorDarkCoreVelocityDispersion} class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nodeOperatorDarkCoreVelocityDispersion)                :: self
    type            (inputParameters                       ), intent(inout) :: parameters
    double precision                                                        :: contractionFactor
    !![
    <inputParameter>
      <name>contractionFactor</name>
      <defaultValue>0.1d0</defaultValue>
      <description>Parameter which mimics the contraction of the dark core as result of gas inflows.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=nodeOperatorDarkCoreVelocityDispersion(contractionFactor)
    !![
      <inputParametersValidate source="parameters"/>
    !!]
  end function darkCoreVelocityDispersionConstructorParameters

  function darkCoreVelocityDispersionConstructorInternal(contractionFactor) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily darkCoreVelocityDispersion} class.
    !!}
    implicit none
    type            (nodeOperatorDarkCoreVelocityDispersion)                :: self
    double precision                                        , intent(in   ) :: contractionFactor
    !![
      <constructorAssign variables="contractionFactor"/>
      <addMetaProperty component="NSC" name="darkCoreRadius"             id="self%darkCoreRadiusID"             isEvolvable="no" isCreator="no" />
      <addMetaProperty component="NSC" name="darkCoreGasMass"            id="self%darkCoreGasMassID"            isEvolvable="no" isCreator="yes"/>
      <addMetaProperty component="NSC" name="darkCoreVelocityDispersion" id="self%darkCoreVelocityDispersionID" isEvolvable="no" isCreator="yes"/>
    !!]
    return
  end function darkCoreVelocityDispersionConstructorInternal

  subroutine darkCoreVelocityDispersionUpdate(self,node)
    !!{
    Update velocity dispersion of the dark core asumming virial equilibrum.
    !!} 
    use :: Galacticus_Nodes                , only : nodeComponentNSC
    use :: Mass_Distributions              , only : massDistributionClass
    use :: Galactic_Structure_Options      , only : componentTypeNuclearStarCluster, massTypeGaseous
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (nodeOperatorDarkCoreVelocityDispersion), intent(inout) :: self
    type            (treeNode                              ), intent(inout) :: node
    class           (nodeComponentNSC                      ), pointer       :: nuclearStarCluster
    class           (massDistributionClass                 ), pointer       :: massDistributionNuclearStarCluster_ 
    double precision                                                        :: darkCoreMass                       , darkCoreRadius                   , &
      &                                                                        darkCoreVelocityDispersion         , nuclearStarClusterGasMassEnclosed
    ! Get the nuclear star cluster component.
    nuclearStarCluster=> node%NSC()

    select type(nuclearStarCluster)
    type is (nodeComponentNSC)
      ! Nothing to do here.
      return
    class default
      nuclearStarCluster=> node              %                      NSC(                     )
      darkCoreMass      =  nuclearStarCluster%             massDarkCore(                     )
      darkCoreRadius    =  nuclearStarCluster%floatRank0MetaPropertyGet(self%darkCoreRadiusID)
    
      massDistributionNuclearStarCluster_ => node                               %massDistribution    (componentTypeNuclearStarCluster,massTypeGaseous)
      nuclearStarClusterGasMassEnclosed   =  massDistributionNuclearStarCluster_%massEnclosedBySphere(darkCoreRadius                                 ) 
    
      if ( darkCoreRadius                  > 0.0d0 &
        & .and.                                    &
        & darkCoreMass                     > 0.0d0 &
        & .and.                                    &
        & nuclearStarClusterGasMassEnclosed>=0.0d0 &
        & ) then
        darkCoreVelocityDispersion= sqrt(                                      &
          &                              (                                     &
          &                               +2.0d0                               &
          &                               *gravitationalConstant_internal      &
          &                               *darkCoreMass                        &
          &                               *(                                   &
          &                                 +darkCoreMass                      &
          &                                 +nuclearStarClusterGasMassEnclosed &
          &                                )                                   &
          &                               )                                    &
          &                               /(                                   &
          &                                  darkCoreMass                      &
          &                                 *darkCoreRadius                    &
          &                                 *self%contractionFactor            &
          &                                )                                   &         
          &                              )
       !![
          <objectDestructor name="massDistributionNuclearStarCluster_"/>
       !!]    
      else 
        darkCoreVelocityDispersion= 0.0d0
      end if

      call nuclearStarCluster%floatRank0MetaPropertySet(self%darkCoreVelocityDispersionID     ,darkCoreVelocityDispersion       )
      call nuclearStarCluster%floatRank0MetaPropertySet(self%darkCoreGasMassID                ,nuclearStarClusterGasMassEnclosed)
    end select
    return  
  end subroutine darkCoreVelocityDispersionUpdate

  subroutine darkCoreVelocityDispersionNodeInitialize(self,node)
    !!{
    Initialize the nuclear star cluster.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentNSC
    implicit none
    class(nodeOperatorDarkCoreVelocityDispersion), intent(inout), target  :: self
    type (treeNode                              ), intent(inout), target  :: node
    class(nodeComponentNSC                      ),                pointer :: nuclearStarCluster

    nuclearStarCluster => node%NSC()
    call self%update(node)
    return
  end subroutine darkCoreVelocityDispersionNodeInitialize
   
  subroutine darkCoreVelocityDispersionSolveAnalytics(self,node,time)
    !!{
    Set the  velocity dispersion of the dark core.
    !!}
    implicit none
    class           (nodeOperatorDarkCoreVelocityDispersion), intent(inout) :: self
    type            (treeNode                              ), intent(inout) :: node
    double precision                                        , intent(in   ) :: time
    !$GLC attributes unused :: time

    call self%update(node)
    return
  end subroutine darkCoreVelocityDispersionSolveAnalytics

  subroutine darkCoreVelocityDispersionNodesMerge(self,node)
    !!{
    Update the velocity dispersion of the dark core after a merger.
    !!}
    implicit none
    class(nodeOperatorDarkCoreVelocityDispersion), intent(inout) :: self
    type (treeNode                              ), intent(inout) :: node

    call self%update(node)
    return
  end subroutine darkCoreVelocityDispersionNodesMerge

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
     double precision :: efficiency       
   contains
     procedure :: differentialEvolutionSolveAnalytics => darkCoreRadiusDifferentialEvolutionSolveAnalytics
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
      <description>The assumed value with Dark Core radius scales with the NSC radius.</description>
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
    type            (nodeOperatorDarkCoreRadius)                 :: self
    double precision                            , intent(in   )  :: efficiency
    !![
    <constructorAssign variables="efficiency"/>
    !!]
    return
  end function darkCoreRadiusConstructorInternal
  
  subroutine darkCoreRadiusDifferentialEvolutionSolveAnalytics(self,node,time)
      !!{
        Compute the nuclear star cluster gas mass rate change.
      !!}
    use :: Galacticus_Nodes , only : nodeComponentNSC             , nodeComponentDarkCore, nodeComponentNSCStandard, &
                &                    nodeComponentDarkCoreStandard, treeNode
    implicit none
    class(nodeOperatorDarkCoreRadius), intent(inout)         :: self
    type (treeNode                  ), intent(inout)         :: node
    double precision                 , intent(in   )         :: time
    class (nodeComponentNSC         ),                pointer:: NSC
    class (nodeComponentDarkCore    ),                pointer:: darkCore
    double precision                                         :: radiusNSC, radiusDarkCore
    !$GLC attributes unused :: time

    NSC      =>  node%     NSC()
    darkCore =>  node%darkCore()

    radiusNSC = NSC%radius()
    if (radiusNSC > 0.0d0) then
      !Check if the Dark Core component has already initialized.
      select type (darkCore)
        type is (nodeComponentDarkCore)
          !Null Dark Core class. Do nothing.
          return
        class default 
          !default class, compute the dark core radius and set.
          radiusDarkCore = self%efficiency * radiusNSC
          call darkCore%radiusSet(radiusDarkCore)
          return
      end select
    end if 
    return
  end subroutine darkCoreRadiusDifferentialEvolutionSolveAnalytics

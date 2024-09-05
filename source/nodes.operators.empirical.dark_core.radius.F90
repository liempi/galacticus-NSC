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
     A node operator class that inserts an empirical model of the formation history of a massive elliptical galaxy. The galaxy is
     assumed to grow in the main branch of the tree with a constant specific star formation rate, such that it mass is given by:
     \begin{equation}
       M_\star(t) = M_{\star,0} \exp(-\phi_\star [t-t_0]),
     \end{equation}
     where $M_{\star,0}=${\normalfont \ttfamily [massStellarFinal]} is the stellar mass in the root node of the tree,
     $\phi_\star=${\normalfont \ttfamily [rateStarFormationSpecific]}, and $t_0$ is the cosmic time at the root node of the tree.
     !!}
     private
     double precision :: efficiency       
   contains
     procedure :: differentialEvolution => darkCoreRadiusDifferentialEvolution
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
  
    subroutine darkCoreRadiusDifferentialEvolution(self,node,interrupt,functioninterrupt,propertyType)
      !!{
        Compute the nuclear star cluster gas mass rate change.
      !!}
    use :: Galacticus_Nodes , only : interruptTask           , nodeComponentNSC             , nodeComponentDarkCore, &
                                 &   nodeComponentNSCStandard, nodeComponentDarkCoreStandard, propertyInactive     , & 
                                 &   treeNode
    implicit none
    class(nodeOperatorDarkCoreRadius), intent(inout), target :: self
    type (treeNode                  ), intent(inout), target :: node
    logical                          , intent(inout)         :: interrupt
    procedure (interruptTask        ), intent(inout), pointer:: functioninterrupt
    integer                          , intent(in   )         :: propertyType
    class (nodeComponentNSC         ),                pointer:: NSC
    class (nodeComponentDarkCore    ),                pointer:: darkCore
    double precision                                         :: radiusNSC, radiusDarkCore

    ! Return immediately if inactive variables are requested.
    if (propertyInactive(propertyType)) return

    NSC      =>  node%     NSC()
    darkCore =>  node%darkCore()

    !Check if there is a NSC in the node.
    select type (NSC)
      type is  (nodeComponentNSC        )
        ! Generic NSC component. Nothing to do here.
        return
      class is (nodeComponentNSCStandard)
        ! Standard NSC class. Get the radius and check if it is positive.
        radiusNSC = NSC%radius()
        if (radiusNSC > 0.0d0) then
          !Check if the Dark Core component has already initialized.
          select type (darkCore)
            type is (nodeComponentDarkCore         )
            ! Null Dark Core class in this node, but there is a NSC class.
            ! Let's interrupt the ODE solver and create a Dark Core standard  
            interrupt = .true.
            ! point to the interrupt method to create the Dark Core
            functionInterrupt => DarkCoreCreate
            return
            class is (nodeComponentDarkCoreStandard)
            ! Standard class, compute the dark core radius and set.
              radiusDarkCore = self%efficiency * radiusNSC
              call darkCore%radiusSet(radiusDarkCore)
              return
          end select
        end if 
    end select 
    return
  end subroutine darkCoreRadiusDifferentialEvolution

  subroutine DarkCoreCreate(node,timeEnd)
  !!{
    Creates the Dark Core via interrupt.
  !!}
    use :: Galacticus_Nodes, only : interruptTask   , nodeComponentDarkCore, nodeComponentDarkCoreStandard, &
          &                         propertyInactive, treeNode
    implicit none
    type (treeNode             ), intent(inout), target  :: node
    double precision            , intent(in   ), optional:: timeEnd
    class(nodeComponentDarkCore),                pointer :: darkCore
    !$GLC attributes unused :: timeEnd
    darkCore => node%darkCore(autoCreate=.true.)
    return 
  end subroutine DarkCoreCreate

!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

  !% Implements a node operator class that triggers destruction of satellites based on their bound mass.

  !# <nodeOperator name="nodeOperatorSatelliteDestructionMassThreshold">
  !#  <description>A node operator class that triggers destruction of satellites based on their bound mass.</description>
  !# </nodeOperator>
  type, extends(nodeOperatorClass) :: nodeOperatorSatelliteDestructionMassThreshold
     !% A node operator class that triggers destruction of satellites based on their mass.
     private
     double precision :: massDestruction, massDestructionFractional
   contains
     !@ <objectMethods>
     !@   <object>nodeOperatorSatelliteDestructionMassThreshold</object>
     !@   <objectMethod>
     !@     <method>massDestroy</method>
     !@     <arguments>\textcolor{red}{\textless type(treeNode)\textgreater} node\arginout</arguments>
     !@     <type>\doublezero</type>
     !@     <description>Compute the mass at which the satellite will be destroyed.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: differentialEvolution => satelliteDestructionMassThresholdDifferentialEvolution
     procedure :: massDestroy           => satelliteDestructionMassThresholdMassDestroy
  end type nodeOperatorSatelliteDestructionMassThreshold
  
  interface nodeOperatorSatelliteDestructionMassThreshold
     !% Constructors for the {\normalfont \ttfamily satelliteDestructionMassThreshold} node operator class.
     module procedure satelliteDestructionMassThresholdConstructorParameters
     module procedure satelliteDestructionMassThresholdConstructorInternal
  end interface nodeOperatorSatelliteDestructionMassThreshold

  ! Submodule-scope pointer to self, used in callback functions.
  class(nodeOperatorSatelliteDestructionMassThreshold), pointer :: self_
  !$omp threadprivate(self_)
  
contains

  function satelliteDestructionMassThresholdConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily satelliteDestructionMassThreshold} node operator class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorSatelliteDestructionMassThreshold)                :: self
    type            (inputParameters                              ), intent(inout) :: parameters
    double precision                                                               :: massDestruction, massDestructionFractional

    !# <inputParameter>
    !#   <name>massDestruction</name>
    !#   <defaultValue>0.00d0</defaultValue>
    !#   <description>The absolute mass below which satellites are destroyed.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>massDestructionFractional</name>
    !#   <defaultValue>0.01d0</defaultValue>
    !#   <description>The fraction of the infall mass below which satellites are destroyed.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    self=nodeOperatorSatelliteDestructionMassThreshold(massDestruction, massDestructionFractional)
    !# <inputParametersValidate source="parameters"/>
    return
  end function satelliteDestructionMassThresholdConstructorParameters

  function satelliteDestructionMassThresholdConstructorInternal(massDestruction, massDestructionFractional) result(self)
    !% Internal constructor for the {\normalfont \ttfamily satelliteDestructionMassThreshold} node operator class.
    implicit none
    type            (nodeOperatorSatelliteDestructionMassThreshold)                        :: self
    double precision                                               , intent(in   )         :: massDestruction, massDestructionFractional
    !# <constructorAssign variables="massDestruction, massDestructionFractional"/>
    
    return
  end function satelliteDestructionMassThresholdConstructorInternal

  subroutine satelliteDestructionMassThresholdDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !% Trigger destruction of a satellite halo based on its bound mass.
    use :: Galacticus_Nodes, only : nodeComponentSatellite
    implicit none
    class           (nodeOperatorSatelliteDestructionMassThreshold), intent(inout), target  :: self
    type            (treeNode                                     ), intent(inout)          :: node
    logical                                                        , intent(inout)          :: interrupt
    procedure       (interruptTask                                ), intent(inout), pointer :: functionInterrupt
    integer                                                        , intent(in   )          :: propertyType
    class           (nodeComponentSatellite                       )               , pointer :: satellite
    double precision                                                                        :: massSatellite
    !$GLC attributes unused :: propertyType
    
    if (.not.node%isSatellite()) return
    satellite     => node     %satellite()
    massSatellite =  satellite%boundMass()
    if     (                                        &
         &   massSatellite > 0.0d0                  &
         &  .and.                                   &
         &   massSatellite < self%massDestroy(node) &
         & ) then
       ! Destruction criterion met - trigger an interrupt.
       interrupt         =  .true.
       functionInterrupt => destructionTrigger
       self_             => self
       return
    end if
    return
  end subroutine satelliteDestructionMassThresholdDifferentialEvolution
  
  subroutine destructionTrigger(node)
    !% Trigger destruction of the satellite by setting the time until destruction to zero.
    use :: Galacticus_Nodes, only : nodeComponentSatellite, treeNode
    implicit none
    type (treeNode              ), intent(inout), target  :: node
    class(nodeComponentSatellite)               , pointer :: satellite

    satellite => node%satellite()
    if (satellite%boundMass() < self_%massDestroy(node)) &
         & call satellite%destructionTimeSet(0.0d0)
    return
  end subroutine destructionTrigger

  double precision function satelliteDestructionMassThresholdMassDestroy(self,node)
    !% Compute the detruction mass for a node.
    use :: Galacticus_Nodes, only : treeNode, nodeComponentBasic
    implicit none
    class(nodeOperatorSatelliteDestructionMassThreshold), intent(inout) :: self
    type (treeNode                                     ), intent(inout) :: node
    class(nodeComponentBasic                           ), pointer       :: basic
    
    satelliteDestructionMassThresholdMassDestroy=self%massDestruction
    if (self%massDestructionFractional > 0.0d0) then
       basic                                        => node%basic()
       satelliteDestructionMassThresholdMassDestroy =  max(                                                       &
            &                                              +      satelliteDestructionMassThresholdMassDestroy  , &
            &                                              +basic%mass                                        ()  &
            &                                              *self %massDestructionFractional                       &
            &                                             )
    end if
    return
  end function satelliteDestructionMassThresholdMassDestroy
!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

  !+    Contributions to this file made by:  Alex Merson.
  
  !% Implements calculations of satellite merging times that are chosen to occur randomly between snapshots.

  !# <satelliteMergingTimescales name="satelliteMergingTimescalesRandom">
  !#  <description>Returns a random timescale for merging.</description>
  !# </satelliteMergingTimescales>   

  type, extends(satelliteMergingTimescalesClass) :: satelliteMergingTimescalesRandom
     !% A class implementing satellite merging timescales that are chosen to occur randomly between snapshots.
     private
   contains
     procedure :: timeUntilMerging => randomTimeUntilMerging
  end type satelliteMergingTimescalesRandom

  interface satelliteMergingTimescalesRandom
     !% Constructors for the {\normalfont \ttfamily random} satellite merging timescale class.                                                                                
     module procedure randomConstructorParameters
  end interface satelliteMergingTimescalesRandom

contains

  function randomConstructorParameters(parameters) result(self)
    !% A constructor for the {\normalfont \ttfamily random} satellite merging timescale class which builds the object from a                                                  
    !% parameter set.
    use Input_Parameters
    implicit none
    type(satelliteMergingTimescalesRandom)                  :: self
    type(inputParameters                 ), intent(inout) :: parameters
    !GCC$ attributes unused :: parameters
    
    self=satelliteMergingTimescalesRandom()
    return
  end function randomConstructorParameters
  
  double precision function randomTimeUntilMerging(self,node,orbit)
    !% Return a randomly chosen timescale for merging satellites.
    use Galacticus_Nodes
    use Kepler_Orbits
    use Satellite_Orbits
    use Pseudo_Random
    use FGSL
    implicit none
    class           (satelliteMergingTimescalesRandom), intent(inout) :: self
    type            (treeNode                        ), intent(inout) :: node
    type            (keplerOrbit                     ), intent(inout) :: orbit
    type            (treeNode                        ), pointer       :: nodeParent
    class           (nodeComponentBasic              ), pointer       :: basicParent, basic
    double precision                                                  :: time       , timeParent
    !GCC$ attributes unused :: self, orbit

    ! Get the parent (a.k.a descendant) node.
    nodeParent             =>  node       %parent
    ! Get basic components of node and its parent.
    basic                  =>  node       %basic ()
    basicParent            =>  nodeParent %basic ()
    ! Extract the times at which the node and its parent exist.
    time                   =   basic      %time  ()
    timeParent             =   basicParent%time  ()
    ! Set time until merging.
    randomTimeUntilMerging =  +nodeParent%hostTree%randomNumberGenerator%uniformSample() &
         &                    *(                                                         &
         &                      +time                                                    &
         &                      -timeParent                                              &
         &                     )                                                         &
         &                    +  timeParent    
    return
  end function randomTimeUntilMerging
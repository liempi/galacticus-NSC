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
Contains a module that implements calculations of the infall radius for cooling calculations.
!!}

module Cooling_Infall_Radii
  !!{
  Provides a class that implements calculations of the infall radius for cooling calculations.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  implicit none
  private

  !![
  <functionClass>
   <name>coolingInfallRadius</name>
   <descriptiveName>Cooling Infall Radius</descriptiveName>
   <description>
    Class providing models of the infall radii for gas cooling in the hot atmosphere surrounding a galaxy.
   </description>
   <default>coolingRadius</default>
   <method name="radius" >
    <description>Return the infall radius for {\normalfont \ttfamily node} (in units of Mpc).</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
   <method name="radiusIncreaseRate" >
    <description>Return the rate at which the infall radius grows for {\normalfont \ttfamily node} (in units of Mpc/Gyr).</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
  </functionClass>
  !!]

end module Cooling_Infall_Radii

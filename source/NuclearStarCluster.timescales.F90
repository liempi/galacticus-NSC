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
Provides a class that implements timescales for star formation.
!!}

module NSC_Timescales
  !!{
  Provides a class that implements calculations of timescales for dark cores.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  private

  !![
  <functionClass>
   <name>NSCTimescale</name>
   <descriptiveName>Timescales for NSC</descriptiveName>
   <description>
    Class providing models of timescales for nuclear star clusters.
   </description>
   <default>dynamicalFrictionTime</default>
   <method name="timescale" >
    <description>Returns the timescale (in Gyr) for Nuclear Star Clusters {\normalfont \ttfamily component}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node </argument>
   </method>
  </functionClass>
  !!]

end module NSC_Timescales

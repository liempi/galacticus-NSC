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

  !+    Contributions to this file made by: Matías Liempi

!!{
Contains a module which provides a class that implements rates for the dark core mass within nuclear star clusters.
!!}

module Dark_Core_Growth_Rates
  !!{
  Provides a class that implements calculations of rates of gas inflows onto nuclear star clusters.
  !!}
  use :: Galacticus_Nodes, only : treeNode
  private

  !![
  <functionClass>
   <name>darkCoreGrowthRates</name>
   <descriptiveName>Rates for the dark core mass within nuclear star clusters</descriptiveName>
   <description>Class providing models of rates for the dark core mass within nuclear star clusters.</description>
   <default>timescale</default>
   <method name="rate" >
    <description>Returns the rate (in units of $\mathrm{M}_\odot$ Gyr$^{-1}$) of stellar mass black holes for falling into the dark core within nuclear star clusters {\normalfont \ttfamily node}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
  </functionClass>
  !!]

end module Dark_Core_Growth_Rates

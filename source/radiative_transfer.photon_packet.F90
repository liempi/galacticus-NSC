!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

!% Contains a module which provides photon packets for radiative transfer calculations.

module Radiative_Transfer_Photon_Packet
  !% Provides a class that implements photon packets for radiative transfer calculations.
  private

  !# <functionClass>
  !#  <name>radiativeTransferPhotonPacket</name>
  !#  <descriptiveName>Radiative Transfer Photon Packets</descriptiveName>
  !#  <description>Class providing photon packets for radiative transfer calculations.</description>
  !#  <default>simple</default>
  !#  <method name="wavelengthSet" >
  !#   <description>Set the wavelength (in \AA) of the photon packet.</description>
  !#   <type>void</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ) :: wavelength</argument>
  !#  </method>
  !#  <method name="wavelength" >
  !#   <description>Get the wavelength (in \AA) of the photon packet.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#  </method>
  !#  <method name="wavelengthMinimumSet" >
  !#   <description>Set the minimum wavelength (in \AA) of the photon packet.</description>
  !#   <type>void</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ) :: wavelength</argument>
  !#  </method>
  !#  <method name="wavelengthMinimum" >
  !#   <description>Get the minimum wavelength (in \AA) of the photon packet.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#  </method>
  !#  <method name="wavelengthMaximumSet" >
  !#   <description>Set the maximum wavelength (in \AA) of the photon packet.</description>
  !#   <type>void</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ) :: wavelength</argument>
  !#  </method>
  !#  <method name="wavelengthMaximum" >
  !#   <description>Get the maximum wavelength (in \AA) of the photon packet.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#  </method>
  !#  <method name="luminositySet" >
  !#   <description>Set the luminosity (in $L_\odot$) of the photon packet.</description>
  !#   <type>void</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ) :: luminosity</argument>
  !#  </method>
  !#  <method name="luminosity" >
  !#   <description>Get the luminosity (in $L_\odot$) of the photon packet.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#  </method>
  !#  <method name="positionSet" >
  !#   <description>Set the position of the photon packet.</description>
  !#   <type>void</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ), dimension(3) :: position</argument>
  !#  </method>
  !#  <method name="position" >
  !#   <description>Get the position of the photon packet.</description>
  !#   <type>double precision, dimension(3)</type>
  !#   <pass>yes</pass>
  !#  </method>
  !#  <method name="directionSet" >
  !#   <description>Set the direction of the photon packet.</description>
  !#   <type>void</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ), dimension(3) :: direction</argument>
  !#  </method>
  !#  <method name="direction" >
  !#   <description>Get the direction of the photon packet.</description>
  !#   <type>double precision, dimension(3)</type>
  !#   <pass>yes</pass>
  !#  </method>
  !# </functionClass>

end module Radiative_Transfer_Photon_Packet
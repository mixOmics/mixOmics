# Author : F.Rohart,  Australian Institute for Bioengineering and Nanotechnology, The University of Queensland, Brisbane, QLD
# created:  23-11-2015
# last modified 18-02-2016
# Copyright (C) 2014
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


.onAttach <- function(libname, pkgname){ packageStartupMessage("\nLoaded mixOmics ",as.character(packageDescription("mixOmics")[["Version"]]),
    
    "\n\nVisit our website http://www.mixOmics.org for more details about our methods.",
    "\nAny bug reports or comments? Send us an email at mixomics at math.univ-toulouse.fr or https://bitbucket.org/klecao/package-mixomics/issues",
    "\n\nThank you for using mixOmics!"
    
    )}

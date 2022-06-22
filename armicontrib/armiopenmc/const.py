# Copyright 2022 TerraPower, LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Some constants associated with OpenMC"""
from enum import Enum

from armi.physics import neutronics

GENERALIZED_BC = [
    neutronics.EXTRAPOLATED,
    neutronics.ZERO_INWARD_CURRENT,
    neutronics.GENERAL_BC,
]


class SolutionType(Enum):
    """Type of flux computed."""

    REAL = 0
    ADJOINT = 1
    REAL_AND_ADJOINT = 2

    @staticmethod
    def fromGlobalFluxOpts(opts):
        """Infer solution type from GlobalFluxOptions object."""
        if opts.real and opts.adjoint:
            return SolutionType.REAL_AND_ADJOINT
        elif opts.real:
            return SolutionType.REAL
        else:
            return SolutionType.ADJOINT


class ProblemType(Enum):
    """Problem type."""

    KEFF = 0
    FIXED_SOURCE = 1

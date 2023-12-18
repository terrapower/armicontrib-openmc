# Copyright 2021 TerraPower, LLC
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

"""
Contains components that prepare, execute, and read output from OpenMC.

Anything wanting OpenMC results should instantiate these objects and use them.
The objects here are intended to be short lived e.g. created, used, and destroyed.

The ``Executer`` orchestrates the tasks of interpreting relevant case settings,
writing an input file, executing OpenMC upon that file, and reading the resultant
output files back in. The results read from the output files are then applied to
the reactor model. In addition to dealing specifically with OpenMC itself, this
sub-classes the :py:class:`GlobalFluxExecuter
<armi:armi.physics.neutronics.globalFlux.globalFluxInterface.GlobalFluxExecuter>`, which
assists in performing geometry transformations, flux-based parameter computation, and
other tasks.

See Also
--------
armi.physics.neutronics.globalFlux.globalFluxInterface.GlobalFluxExecuter : ARMI-provided base class to ``GlobalFluxExecuter``.

~armicontrib.armiopenmc.executionOptions.OpenMCOptions : Options that control the behavior of the Executer and its components.

~armicontrib.armiopenmc.schedulers.OpenMCCoreEvaluation : One particular client of this code
that runs OpenMC on a core.
"""

import os
import subprocess

from armi.utils import directoryChangers
from armi.utils import outputCache
from armi.reactor import composites
from armi.physics import executers
from armi.physics.neutronics.globalFlux import globalFluxInterface

from . import inputWriters
from . import outputReaders
from . import executionOptions
from armi.utils import codeTiming

import openmc.lib


class OpenMCExecuter(globalFluxInterface.GlobalFluxExecuter):
    """
    A short-lived object that coordinates the prep, execution, and processing of a flux solve.

    This is an implementation of ARMI's ``DefaultExecuter`` and thus follows that
    pattern closely.

    Parameters
    ----------
    options: executionOptions.OpenMCOptions
        run settings
    armiObj : composite.ArmiObject
        The object representing the scope of the OpenMC run (e.g. the Core or SFP, etc.)
    """

    def _performGeometryTransformations(self):
        """Disable geometry transformations for MC"""
        pass

    @codeTiming.timed
    def writeInput(self):
        """
        Write the input for OpenMC

        Compose an InputWriter from options and send it off to do its work.

        Actually just builds the OpenMC model via the OpenMC API for now.
        """
        writer = inputWriters.OpenMCWriter(self.r, self.options)
        writer.write()

    @codeTiming.timed
    def _execute(self):
        """
        Execute OpenMC on an existing input model. Produce output files.
        """
        globalFluxInterface.GlobalFluxExecuter._execute(self)

        if self.options.executablePath is None:
            raise ValueError(
                f"Cannot find executable at {self.options.executablePath}. " f"Update run settings."
            )

        openmc.run(threads=self.options.nOMPThreads,
                   mpi_args=['mpiexec','-n',str(self.options.nMPIProcesses), '--bind-to', 'none'],
                   openmc_exec=self.options.executablePath)

        return True

    @codeTiming.timed
    def _readOutput(self):
        """Read output."""
        reader = outputReaders.OpenMCReader(self.options)
        reader.apply(self.r)
        return reader

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
Code that schedules OpenMC executions during ARMI cases.

This module contains the :py:class:`OpenMCInterface` class, which integrates OpenMC runs
into ARMI runs using the ARMI :py:class:`~armi.interfaces.Interface` API.  The
:py:class:`OpenMCInterface` is implemented on top of the :py:class:`OpenMCCoreEvaluation`
class, which handles the actual task of running OpenMC. This separation between
scheduling the OpenMC analysis and performing the analysis aids in code re-use and
flexibility; rather than binding the process of running OpenMC analysis to the Interface
interaction hooks, the ``OpenMCCoreEvaluation`` can be used as a standard MPI action in
custom scripts, or be driven directly by a reactivity coefficients interface, for
instance.

.. pyreverse:: armicontrib.armiopenmc.schedulers -A
    :align: center
    :width: 90%
"""

from armi import runLog
from armi import interfaces
from armi.physics import neutronics
from armi.nuclearDataIO import xsLibraries
from armi import mpiActions

from . import executers
from . import executionOptions


class OpenMCInterface(interfaces.Interface):
    """ "Schedules activities related to OpenMC during ARMI run."""

    name = "openmc"  # name is required for all interfaces

    def interactEveryNode(self, cycle=None, node=None):
        """
        Run OpenMC on the Core.

        Builds a OpenMCCoreEvaluation (MpiAction) and invokes it.
        """
        runLog.info("Running OpenMC to update keff, flux, and power.")

        runner = OpenMCCoreEvaluation()
        # Run on any potential worker mpi nodes
        runner.broadcast()
        # Run on this process as well
        runner.invoke(self.o, self.r, self.cs)


class OpenMCCoreEvaluation(mpiActions.MpiAction):
    """
    Sets up a OpenMC run on the Core object and executes it.

    This is one specific client of the OpenMC Executor that happens to be used
    in a normal ARMI run to run OpenMC on the Core.

    MpiActions have access to the operator and the reactor and can therefore
    reach in as appropriate to select which blocks to execute on in a
    analysis-specific way.

    This default implementation just runs one case, on the ``core``, though it
    can be extended to run other cases (e.g. on branch-searched copies of the
    core, the core + sfp, etc.)
    """

    def invokeHook(self):
        """Perform OpenMC calculation for the blocks assigned to this process."""
        for obj in self.mpiIter(selectObjsToRun(self.o)):
            executer = self._buildExecuterForObj(obj)
            executer.run()

        if self.parallel:
            self.gather()
            self.r.syncMpiState()

    def _buildExecuterForObj(self, obj):
        """Build options and executers for a block."""
        opts = executionOptions.OpenMCOptions(
            label=f"openmc-{obj.getName()}-{self.r.p.cycle}-{self.r.p.timeNode}"
        )

        opts.fromReactor(self.r)
        opts.fromUserSettings(self.cs)

        return executers.OpenMCExecuter(opts, obj)


def selectObjsToRun(o):
    """
    Choose objects that will be passed for OpenMC analysis.

    We choose the reactor rather than the core because the globalflux code makes
    some assumptions that it's using the reactor for geometry conversions (should
    be changed in the framework!)
    """
    return [o.r]

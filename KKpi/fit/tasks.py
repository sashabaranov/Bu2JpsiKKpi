import pyutilib
from B2JpsiKKpi12 import *


class SelectionTestTask(pyutilib.workflow.Task):
    def __init__(self, *args, **kwds):
        """Constructor."""
        pyutilib.workflow.Task.__init__(self, *args, **kwds)
        self.inputs.declare('inputdata')
        self.inputs.declare('year')

    def execute(self):
    """Compute the sum of the inputs."""
        inputdata = []
        configure(self.inputdata, self.year)
        run(-1)

if __name__ == "__main__":
    S = SelectionTestTask()
    w = pyutilib.workflow.Workflow()
    w.add(S)

    w(inputdata=['/tmp/albarano/00024671_00000002_1.psix.mdst'], year='2012')
from pyutilib import workflow

class FitHistoTask(pyutilib.workflow.Task):
    def __init__(self, *args, **kwds):
        """Constructor."""
        pyutilib.workflow.Task.__init__(self, *args, **kwds)
        self.inputs.declare('hist')
        self.inputs.declare(’y’)
        self.outputs.declare(’z’)

    def execute(self):
        """Compute the sum of the inputs."""
        self.z = self.x + self.y

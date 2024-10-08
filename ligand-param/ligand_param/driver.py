from stagegausrotation import StageGaussianRotation
from stageinitialize import StageInitialize

class Driver:
    def __init__(self):
        self.stages = []
        return
    
    def add_stage(self, stage):
        self.stages.append(stage.append_stage(stage))
        return
    
    def execute(self, dry_run=False):
        for stage in self.stages:
            stage.execute(dry_run=dry_run)
        return
    
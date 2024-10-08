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
    

if __name__ == "__main__":

    driver = Driver()
    driver.add_stage(StageInitialize("Initialize", base_name="F3KRP"))
    driver.add_stage(StageGaussianRotation("Rotation1", alpha=0.0, beta=0.0))
    driver.add_stage(StageGaussianRotation("Rotation2", alpha=0.1, beta=0.1))
    driver.execute()
    
    print("Done!")
    
    exit(0)
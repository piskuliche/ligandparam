class Driver:
    def __init__(self):
        self.stages = []
        return
    
    def add_stage(self, stage):
        self.stages.append(stage.append_stage(stage))
        return
    
    def execute(self, dry_run=False):
        for stage in self.stages:
            try:
                stage._execute(dry_run=dry_run)
            except Exception as e:
                print(f"Error in stage {stage.name}: {e}")
                print("Exiting")
                raise e
        return
    
    
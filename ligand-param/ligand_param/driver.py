class Driver:
    def __init__(self):
        """ Initialize the Driver class object.

        This class is the base class for all parametrizations. It is designed to be subclassed, and the subclass should
        implement the stages that are needed to complete the parametrization.


        """
        self.stages = []
        return
    
    def add_stage(self, stage):
        """ Add a stage to the list of stages to run.

        Adding a stage is done by appending a stage object to the list of stages to run. Adding a stage
        does not run the stage, it only adds it to the list of stages to run. To run the stages, the execute
        method must be called.

        Stages are designed using the AbstractStage class, and should be subclassed to implement the desired
        behavior. The stages should be added in the order that they should be run.
    

        Parameters
        ----------
        stage : Stage
            The stage object to add to the list of stages to run.
        
        Returns
        -------
        None

        """
        self.stages.append(stage.append_stage(stage))
        self.list_stages()
        return
    
    def execute(self, dry_run=False):
        """ Execute the stages in the list of stages to run. 
        
        This function executes the stages in the list of stages to run. The stages are executed in the order that they
        were added to the list. If a stage fails, the function will print an error message and exit. The stages are
         
        Parameters
        ----------
        dry_run : bool, optional
            If True, the stages will not be executed, but the function will print the commands that would be executed.
        
        Returns
        -------
        None
        """
        for stage in self.stages:
            try:
                stage.execute(dry_run=dry_run)
            except Exception as e:
                print(f"Error in stage {stage.name}: {e}")
                print("Exiting")
                raise e
        return
    
    def clean(self):
        """ Clean up the files created by the stages. 
        
        This function cleans up the files created by the stages. The stages are executed in the
        reverse order that they were added to the list. If a stage fails, the function will print an error message and
        exit and be skipped.

        Returns
        -------
        None
        """

        for stage in reversed(self.stages):
            try:
                stage.clean()
            except NotImplementedError:
                print(f"Clean method not implemented for stage {stage.name}")
                print("Skipping...")
                continue
            except Exception as e:
                print(f"Error in stage {stage.name}: {e}")
                print("Exiting")
                raise e
        return
    
    def list_stages(self):
        """ Print out the list of stages to run. 
        
        This function prints out the list of stages that are in the list of stages to run. The stages are printed in the
        order that they were added to the list.

        Returns
        -------
        None
        """
        print("********************************")
        print("List of Stages to Run")
        for stage in self.stages:
            print(f"-->{stage.name}")
        return
    
    def remove_stage(self, stage_name):
        """ Remove a stage from the list of stages to run. 
        
        This function removes a stage from the list of stages to run. If the stage is not in the list, the function will
        print an error message and exit.

        Parameters
        ----------
        stage_name : str
            The name of the stage to remove from the list of stages to run.
        
        Returns
        -------
        None
        """
        for stage in self.stages:
            if stage.name == stage_name:
                self.stages.remove(stage)
                print(f"Stage {stage_name} removed.")
                self.list_stages()
                return
        print(f"Stage {stage_name} not found in list of stages.")
        return
    
    def insert_stage(self, newstage, stage_name):
        """ Insert a stage into the list of stages to run before the specified stage.
        
        This function inserts a stage into the list of stages to run before the specified stage. If the specified stage
        is not in the list, the function will print an error message and exit.
        
        Parameters
        ----------
        stage_name : str
            The name of the stage to insert into the list of stages to run.
        """
        idx = -1
        for stage in self.stages:
            if stage.name == stage_name:
                idx = self.stages.index(stage)
                self.stages.insert(idx, newstage)
                print(f"Stage {newstage.name} inserted before {stage_name}")
                self.list_stages()
                return
        print(f"Stage {stage_name} not found in list of stages.")
        return

    
    
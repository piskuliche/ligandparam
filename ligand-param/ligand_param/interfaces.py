class SimpleInterface:

    def __init__(self) -> None:
        return
    
    def set_method(self, method):
        self.method = method
        return
    
    def call(self, **kwargs):
        import subprocess
        run = False
        if "run" in kwargs:
            run = kwargs['run']
            del kwargs['run']
    
        command = [self.method]
        for key, value in kwargs.items():
            if value is not None:
                command.extend([f'-{key}', str(value)])
        if run:
            subprocess.run(command)
        else:
            print(command)
        return
   
class Antechamber(SimpleInterface):
    def __init__(self) -> None:
        self.set_method('antechamber')
        return
    
class ParmChk(SimpleInterface):
    def __init__(self) -> None:
        self.set_method('parmchk2')
        return
    
class Leap(SimpleInterface):
    def __init__(self) -> None:
        self.set_method('tleap')
        return
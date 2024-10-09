class SimpleInterface:
    """ This class is a simple interface to call external programs.
    
    This class is a simple interface to call external programs. It is designed to be subclassed
    and the method attribute set to the desired program. The call method will then call the program
    with the specified arguments.
    
    """

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
    """ This class is a simple interface to call the Antechamber program. """
    def __init__(self) -> None:
        self.set_method('antechamber')
        return
    
class ParmChk(SimpleInterface):
    """ This class is a simple interface to call the ParmChk program. """
    def __init__(self) -> None:
        self.set_method('parmchk2')
        return
    
class Leap(SimpleInterface):
    """ This class is a simple interface to call the Leap program. """
    def __init__(self) -> None:
        self.set_method('tleap')
        return
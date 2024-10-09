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
        dry_run = False
        if "dry_run" in kwargs:
            dry_run = kwargs['dry_run']
            del kwargs['dry_run']
    
        command = [self.method]
        for key, value in kwargs.items():
            if key is "inp_pipe":
                command.append(f'< {value}')
            elif key is "out_pipe":
                command.append(f'> {value}')
            else:
                if value is not None:
                    command.extend([f'-{key}', str(value)])
        if dry_run:
            print(command)
        else:
            subprocess.run(command)
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

class Gaussian(SimpleInterface):
    """ This class is a simple interface to call the Gaussian program. """
    def __init__(self) -> None:
        self.set_method('g16')
        return
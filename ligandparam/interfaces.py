class SimpleInterface:

    def __init__(self) -> None:
        """ This class is a simple interface to call external programs.

        This class is a simple interface to call external programs. It is designed to be subclassed
        and the method attribute set to the desired program. The call method will then call the program
        with the specified arguments.
        
        """
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
        shell=False
        for key, value in kwargs.items():
            if key is "inp_pipe":
                command.extend([f'<', str(value)])
                shell=True
            elif key is "out_pipe":
                command.extend([f'>', str(value)])
                shell=True
            else:
                if value is not None:
                    command.extend([f'-{key}', str(value)])

        if dry_run:
            print(command)
        else:
            print("Executing command")
            proc = subprocess.run(command, shell=shell, encoding='utf-8', stdout=subprocess.PIPE)
            for line in proc.stdout.split('\n'):
                print(f"[{command[0]}] -> {line}")
            print(f"Command {command} executed")
        return
   
class Antechamber(SimpleInterface):
    def __init__(self) -> None:
        """ This class is a simple interface to call the Antechamber program. """
        self.set_method('antechamber')
        return
    
class ParmChk(SimpleInterface):
    def __init__(self) -> None:
        """ This class is a simple interface to call the ParmChk program. """
        self.set_method('parmchk2')
        return
    
class Leap(SimpleInterface):
    def __init__(self) -> None:
        """ This class is a simple interface to call the Leap program. """
        self.set_method('tleap')
        return

class Gaussian(SimpleInterface):
    def __init__(self) -> None:
        """ This class is a simple interface to call the Gaussian program. """
        self.set_method('g16')
        return
    
    def call(self, **kwargs):
        """ This function calls the Gaussian program with the specified arguments,
        however, it works slightly differently than the other interfaces. The Gaussian
        interface for some reason isn't compatible with the subprocess.run() function
        so we instead write a bash script to call the program and then execute the script."""
        import subprocess
        dry_run = False
        if "dry_run" in kwargs:
            dry_run = kwargs['dry_run']
            del kwargs['dry_run']
    
        command = [self.method]
        shell=False
        for key, value in kwargs.items():
            if key is "inp_pipe":
                command.extend([f'<', str(value)])
                shell=True
            elif key is "out_pipe":
                command.extend([f'>', str(value)])
                shell=True
            else:
                if value is not None:
                    command.extend([f'-{key}', str(value)])
                    
        self.write_bash(' '.join(command))
        bashcommand = 'bash temp_gaussian_sub.sh'

        if dry_run:
            print(bashcommand)
        else:
            print("Executing command")
            subprocess.run(bashcommand, shell=shell)
            print(f"Command Executed:")
            print(f"{' '.join(bashcommand)}")
        return
    
    def write_bash(self, command):
        """ This function writes a bash script to call the Gaussian program
        with the specified arguments. """
        with open('temp_gaussian_sub.sh', 'w') as f:
            f.write('#!/bin/bash\n\n')
            f.write(command)
        return
    

            
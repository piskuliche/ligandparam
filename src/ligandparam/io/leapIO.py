class LeapWriter:
    def __init__(self, name):
        """ Class to write tleap input files

        This class writes a simple tleap input file, with the ability to 
        add leaprc files and lines to the file. This is a fairly simple interface,
        and does not currently check for errors in the input file.


        Parameters
        ----------
        name : str
            The name of the LeapWriter object

        """
        
        self.name = name
        self.leaprc = []
        self.lines = []
        return
    
    def add_leaprc(self, leaprc):
        """ Add a leaprc file to the tleap input file

        Parameters
        ----------
        leaprc : str
            The name of the leaprc file to add
        """
        self.leaprc.append(leaprc)
        return
    
    def gen_leap(self):
        """ Generate the tleap input file """
        raise DeprecationWarning("gen_leap is deprecated. Use write instead.")
        for leap in self.leaprc:
            self.lines.append(f"source {leap}")
        self.lines.append("")
        return
    
    def add_line(self, line):
        """ Add a line to the tleap input file. 
        
        These lines are printed AFTER the leaprc files.
        
        Parameters
        ----------
        line : str
            The line to add to the tleap input file
        """
        self.lines.append(line)
        return

    def remove_line(self, string):
        """ Remove a line from the tleap input file. 
        
        Parameters
        ----------
        line : str
            The line to remove from the tleap input file
        """
        for line in self.lines:
            if string in line:
                self.lines.remove(line)
        return
    
    def write(self, out_filepath):
        """ Write the tleap input file to disk.

        This method writes the tleap input file to disk. The filename
        is generated from the name of the LeapWriter object.

        """
        with open(out_filepath, 'w') as f:
            for line in self.leaprc:
                f.write(f"source {line}\n")
            for line in self.lines:
                f.write(f"{line}\n")
    

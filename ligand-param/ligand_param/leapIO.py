class LeapWriter:
    def __init__(self, name):
        self.name = name
        self.leaprc = []
        self.lines = []
        return
    
    def add_leaprc(self, leaprc):
        self.leaprc.append(leaprc)
        return
    
    def gen_leap(self):
        for leap in self.leaprc:
            self.lines.append(f"source {leap}")
        self.lines.append("")
        return
    
    def add_line(self, line):
        self.lines.append(line)
        return
    
    def write(self):
        with open(f"tleap.{self.name}.in", 'w') as f:
            for line in self.lines:
                f.write(f"{line}\n")
    

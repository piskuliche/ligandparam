"""HPC / local job-script helpers for multi-RESP fragment MD."""

from __future__ import annotations

class Computer(object):
    """ Base class for computer objects """
    def __init__(self,numnodes):
        """ Initialize the computer object 
        
        Parameters
        ----------
        numnodes : int
            The number of nodes to use
        """
        self.mpirun=""
        self.num_nodes = numnodes
        self.exclude = ""
        self.cores_per_node = 1
        self.amberhome="${AMBERHOME}"
        self.gpu=False
        self.array = []
        self.array_max_running = 1

        
    def get_array(self):
        alist = list(sorted(set(self.array)))
        rs=[]
        if len(alist) > 0:
            rs = [ (alist[0],alist[0]) ]
            for a in alist[1:]:
                if a == rs[-1][1]+1:
                    rs[-1] = ( rs[-1][0], a )
                else:
                    rs.append( (a,a) )
        sarr = []
        for r in rs:
            if r[0] != r[1]:
                sarr.append( "%i-%i"%(r[0],r[1]) )
            else:
                sarr.append( "%i"%(r[0]) )
        return ",".join(sarr)


    def write_array(self,fh):
        """ Write the array to the file handle 
        
        Parameters
        ----------
        fh : file handle
            The file handle to write to
        """
        if len(self.array) > 0:
            fh.write("#SBATCH --array=%s"%(self.get_array()))
            if self.array_max_running > 0:
                fh.write("%%%i"%(self.array_max_running))
            fh.write("\n")


    def unset_amberhome(self):
        """ Unset the amberhome variable """
        self.amberhome = None

    def set_exclude(self,x):
        """ Set the exclude variable """
        self.exclude = x

    def use_gpu(self,x=True):
        """ Use the gpu 
        
        Parameters
        ----------
        x : bool, optional
            If True, use the gpu
        """
        self.gpu=x

    def get_num_cores(self):
        """ Get the number of cores """
        return self.num_nodes * self.cores_per_node



class BASH(Computer):
    """ Generates a bash script """
    def __init__(self,numnodes):
        """ Initialize the bash object

        Parameters
        ----------
        numnodes : int
            The number of nodes to use
        """

        Computer.__init__(self,numnodes)
        self.mpirun="mpirun -n %i"%(self.get_num_cores())
                    
    def open(self,fname):
        """ Open the bash file

        Parameters
        ----------
        fname : str
            The name of the file to open
        """
        fh = open(fname,"w")
        fh.write("#!/bin/bash\n\n")
        return fh

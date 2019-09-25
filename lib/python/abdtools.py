"""
Python data structures and methods for dealing with Almost Block Diagonal
(ABD) style matrices.
"""
import colors as c

import sys
import numpy as np

import scipy.sparse as sp

NOTHING = object()

class ABDMatrix:
    """
    Structure for holding data in an ABDMatrix.

    - This provides ABDMatrix.to_coo_matrix() and ABDMatrix.to_dense() so that either Scipy's sparse matrix tools or numpy arrays can be used from this.
    """

    def __init__(self,fname=None,ftype="abd"):
        if fname:
            if ftype == 'abd':
                self.read_abd_file(fname)
            else:
                print('ABDMatrix only supports file type "abd"')
        else:
            self.has_data = False

    def read_abd_file(self, fname=None):
        self.filename = fname
        with open(fname) as f:
            firstline = f.readline().split()

            self.nrow_top = int(firstline[0])
            self.ncol_top = int(firstline[1])
            self.siztop = self.nrow_top*self.ncol_top
            self.topblk = read_column_major_from_file(f,self.nrow_top,self.ncol_top)

            line = f.readline().split()
            self.nblocks = int(line[0])
            self.nrow_blk = int(line[1])
            self.ncol_blk = int(line[2])
            self.sizblk = self.nrow_blk*self.ncol_blk
            self.blk = []
            for i in range(self.nblocks):
                self.blk.append(read_column_major_from_file(f,self.nrow_blk,self.ncol_blk))
            print("hello")

            line = f.readline().split()
            self.nrow_bot = int(line[0])
            self.ncol_bot = int(line[1])
            self.sizbot = self.nrow_bot*self.ncol_bot
            self.botblk = read_column_major_from_file(f,self.nrow_bot,self.ncol_bot)

        self.nnz = self.siztop + self.nblocks*self.sizblk + self.sizbot
        self.dim = self.nrow_top+self.nblocks*self.nrow_blk+self.nrow_bot
        self.has_data = True

    def clear_data(self):
        """Nullify object's data."""
        self.nrow_top = None
        self.ncol_top = None
        self.nsiz_top = None
        self.topblk   = None
        self.nblocks  = None
        self.nrow_blk = None
        self.ncol_blk = None
        self.nsiz_blk = None
        self.blk      = None
        self.nrow_bot = None
        self.ncol_bot = None
        self.nsiz_bot = None
        self.botblk   = None
        self.has_data = False
        pass

    def set_data(self):
        """
        TODO: This method sets all of the data in an ABDMatrix. Input args will
        (necessarily) be similar to those in the Fortran ABD functions.
        """
        pass

    def to_coo_matrix(self):
        """
        Function that converts an ABDMatrix to scipy.sparse.coo_matrix.
        """
        if not self.has_data:
            raise Exception("ABDMatrix currently has no data.")

        data  = np.zeros(self.nnz)
        ivals = np.zeros(self.nnz,dtype=int)
        jvals = np.zeros(self.nnz,dtype=int)
        count = 0
        # loop over top block
        for i in range(self.nrow_top):
            for j in range(self.ncol_top):
                ivals[count] = i
                jvals[count] = j
                data[count] = self.topblk[i,j]
                count = count + 1
        # loop over interior blocks
        for n in range(self.nblocks):
            for i in range(self.nrow_blk):
                for j in range(self.ncol_blk):
                    ivals[count] = self.nrow_top + n*self.nrow_blk + i
                    jvals[count] = n*(self.ncol_blk-self.ncol_top) + j
                    data[count] = self.blk[n][i,j]
                    count = count + 1
        # loop over bot block
        for i in range(self.nrow_bot):
            for j in range(self.ncol_bot):
                ivals[count] = self.dim - self.nrow_bot + i
                jvals[count] = self.dim - self.ncol_bot + j
                data[count] = self.botblk[i,j]
                count = count + 1
        assert (count == self.nnz) # just make sure all entries accounted for
        return sp.coo_matrix((data,(ivals,jvals)),shape=[self.dim,self.dim])

    def to_dense(self):
        """
        Function that converts an ABDMatrix to a dense format.
        """
        return sp.spmatrix.todense(self.to_coo_matrix())

## END OF class EbacoliData

############################################################################
############################################################################

#### Functions that use ABDMatrix objects

def read_column_major_from_file(file_handle,nrow,ncol):
    """
    Read column major ordered data from file, single entry per line.

    - file_handle needs to be queued up to where the data begins on entry, it is
      left where the data stops on exit

    """
    data = np.zeros((nrow,ncol))
    for j in range(ncol):
        for i in range(nrow):
            number = float(file_handle.readline().split()[0])
            data[i][j] = number
    return data

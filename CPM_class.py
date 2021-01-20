'''
   A basic cellular potts class. 
   Call self.MonteCarloStep() to perform 
   an update of the lattice
'''
import numpy as np
from random import randint, sample



class Potts:

    def __init__(self, Q, LATTICE, J='None', NEIGHBORHOOD='Moore'):
        '''
        Q is an integer. The number of != cells is Q-1
        LATTICE is an 2D-array with Q integer values in [0,Q]

        J is a (2*Q+2)xQ matrix of coefficients \n
        including the lagragian coefficients, cell-cell adhesion\n
        and target area/perimeter for each cell type.
        J(i,j) is the interaction between cell type i and j.

        NEIGHBORHOOD is the type of neighborhood to use \n
        by default it is the Moore type (8 closest neighbors) 
        '''
        self.latt = LATTICE.astype(int)       
        self.shape = LATTICE.shape  
        self.q    = Q+1           # Number of cells IDs (spins). 
                                  # ID 0 is kept for the medium.
        self._neighborhood = NEIGHBORHOOD

        self._Cells = [[0,0] for _ in range(Q)]
        self._GetAreasAndPerimeters()
        self._GetParameters(J)   # Initializes the CPM's parameters
        

        self.StepsCount = 0     # Count of the MCS
        
    def _GetParameters(self, J):

        if (type(J)==str):
            A,L = self._Cells[1]
            self._J = np.ones((self.q,self.q))*A
            self._Lam_area = np.ones((self.q,))/L
            self._Lam_peri = np.ones((self.q,))/L
            self._Target   = np.array([40,20]).reshape((2,))

        else:
            print("Initializing parameters from given array.")
            self._J = J[:self.q,:]                # Cell interactions 
            self._Lam_area = J[self.q,:]          # Lagragian coefficients
                                                  # for the area of cells
            self._Lam_peri = J[self.q:self.q+2,:] # Lagragian coefficients
                                                  # for the perimeter of cells
            self._Target   = J[self.q+2:,0]       # Target dimensions [A_0, L_0]

        self._Beta = 1
        

    def _Area_Perimeter(self, L, i):
        ''' 
        Computes area and perimeter of the cell i in lattice L 
        Perimeter's formula is probably not reliable for now
        '''
        copy_L = L==i

        area = copy_L.sum()
        perimeter = np.logical_or( (copy_L[:,1:] != copy_L[:,:-1])[1:,:],
                                   (copy_L[1:,:] != copy_L[:-1,:])[:,1:] ).sum() + 1
        return area, perimeter


    def _GetAreasAndPerimeters(self):
        ''' 
        Updates the cells area and perimeter
        for the current lattice state
        '''
        for k in range(self.q-1):
            self._Cells[k] = self._Area_Perimeter(self.latt, k)          
            
        
    def _Neighbors(self, x,y):
        '''
        Return a vector of the IDs of the neighbors \n
        of the vertex (x,y) according to the method \n
        defined by self._neighborhood 
        '''
        n,m = self.shape
        if self._neighborhood == 'Moore':
            cols = np.array([(y-1)%m, y%m, (y+1)%m])
            rows = np.array([(x-1)%n, x%n, (x+1)%n])
            
            return np.delete(self.latt[rows][:,cols].flatten(), 4)

        if self._neighborhood=='VonNeumann':
            indexes = [(x,(y+1)%m), ((x-1)%n,y), ((x+1)%n,y), (x,(y-1)%m)]
            return self.latt[indexes[:][0],indexes[:][1]].flatten()            

        else:
            raise ValueError("The type of neighborhood {} is not implemented yet".format(self._neighborhood))
        
    
    



    def _Diff_Energy(self, x,y, z):
        '''
        Compute the difference in energy between the current state
        of self.latt and the proposed update.
        The proposed update is to change the ID of vertex (x,y) 
        into z.
        The energy from a neighbors is counted if it belongs to a
        different cell than the current vertex, i.e. it represents
        the boundary energy of the interacting cells

        The formula for Hamiltonian is taken from the paper:
        https://doi.org/10.3389/fphy.2018.00061

        Delta_E = new_energy - old_energy
        '''

        old_z = self.latt[x,y]
        Delta_E = 0.0
        
        # Neighbors' interactions
        neigh = self._Neighbors(x,y)

        for k in neigh:
            Delta_E -= self._J[z, k]*float(k!=z)
            Delta_E += self._J[old_z, k]*float(k!=old_z)

        # Constraints (e.g. area, perimeter, chemotaxis...)
        for k in [z, old_z]:
            A,P = self._Area_Perimeter(self.latt, k)
            Delta_E += self._Lam_area[k]*(1.0 + 2.0*A - 2.0*self._Target[0])

            # TODO: Add the perimeter constraint

        return Delta_E

    def _Accept_Or_Not(self, x,y, z):
        '''
        Accept the proposed new state z with probabiliy
        p = exp(-Delta_E*Beta) if Delta_E > 0 or
        p = 1 if Delta_E <= 0
        '''
        if self.latt[x,y]==z:
            return

        Delta_E = self._Diff_Energy(x,y,z)
        if np.exp(-Delta_E*self._Beta) > np.random.random_sample():
            self.latt[x,y] = z
            self._NbOfChanges += 1.0

    def MonteCarloStep(self):
        ''' 
        Compute the new state of the lattice
        '''

        self._NbOfChanges = 0.0
        self.StepsCount += 1

        n,m = self.shape
        N = n*m

        # Shuffles the vertices for the random walk
        randomized_vertices = sample(range(N), N)
        for vertex in randomized_vertices:
            x,y = vertex % n, vertex // n
            # Chose a trial state among the neighbors
            # and test it
            z   = np.random.choice(self._Neighbors(x,y))
            self._Accept_Or_Not(x,y, z)

        if self.StepsCount % 5 == 0:
            print("{}% of changes at step {}".format(100.0*self._NbOfChanges/N,
                                                     self.StepsCount))
                
        
        

        
        

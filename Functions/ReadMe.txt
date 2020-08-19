'MAC.m' contains the MAC routine that conducts RREF first then builds a dissimilarity matrix D based on our proposed measure, before conducting spectral clustering on D.

'HMAC' is the hierarchical version of MAC that does not involve any use of spectral clustering at all. It can start with an adjacency matrix with no prior information, or start with some prior information about data connectiveness either through labelling constraint or through sparse representation. Then at every iteration, two points / components are merged based on the dissimilarity (use the proposed measure) between them.

'MACadj.m' contains a variant of MAC where the dissimilarity matrix is constructed purely based on RREF, without using the proposed dissimilarity measure.  

'MACw.m' is a variant of 'MAC.m' that takes as inputs one extra variable -- W. We obtain W here through sparse representation of each point through other points (code can be found in 'buildW.m'). Then, the adjacency matrix is formed based on W and everything else follows as in 'MAC.m'.
   

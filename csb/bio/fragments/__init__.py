
class FragmentMatch(object):
    
    def __init__(self, id, qstart, qend, probability, rmsd, qlength):
        
        self.id = id
        self.qstart = qstart
        self.qend = qend
        self.probability = probability
        self.rmsd = rmsd
        self.qlength = qlength
        
    @property
    def length(self):
        
        return self.qend - self.qstart + 1
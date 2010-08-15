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
    
class PredictionContainer(object):
    
    def __init__(self, target, isites_prediction, hmm_prediction, combined_prediction):
        
        self.target = target
        
        self.isites = isites_prediction
        self.hmm = hmm_prediction
        self.combined = combined_prediction 
        
class Prediction(object):
    
    def __init__(self, alignment, coordinates):
        
        self.alignment = alignment
        self.coordinates = coordinates     
import numpy
import csb.pyutils


class FragmentMatch(object):
    
    def __init__(self, id, qstart, qend, probability, rmsd, tm_score, qlength):
        
        self._id = id
        self._qstart = qstart
        self._qend = qend
        self._probability = probability        
        self._rmsd = rmsd
        self._tm_score = tm_score
        self._qlength = qlength
        
    @property
    def id(self):
        return self._id

    @property
    def qstart(self):
        return self._qstart
    
    @property
    def qend(self):
        return self._qend
    
    @property
    def qlength(self):
        return self._qlength    
    
    @property
    def rmsd(self):
        return self._rmsd
    
    @property
    def tm_score(self):
        return self._tm_score
    
    @property
    def probability(self):
        return self._probability
            
    @property
    def length(self):
        return self.qend - self.qstart + 1
    
    @property
    def source_id(self):
        raise NotImplementedError()
        
    @property
    def start(self):
        raise NotImplementedError()

    @property
    def end(self):
        raise NotImplementedError()      
    
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
        
class Target(csb.pyutils.AbstractNIContainer):
    
    def __init__(self, id, length, residues, overlap=None):
    
        self._id = id
        self._accession = id[:4]
        self._chain_id = [4]
        self._length = length
        self._overlap = overlap
        
        self._assignments = csb.pyutils.ReadOnlyCollectionContainer(type=Assignment)
            
        resi = [TargetResidue(native) for native in residues]
        self._residues = csb.pyutils.ReadOnlyCollectionContainer(items=resi, 
                                            type=TargetResidue, start_index=1)
        
    @property
    def _children(self):
        return self._residues        
    
    @property
    def id(self):
        return self._id
    
    @property
    def accession(self):
        return self._accession
    
    @property
    def chain_id(self):
        return self._chain_id
    
    @property
    def max_overlap(self):
        return self._overlap
    
    @property
    def length(self):
        return self._length
    
    @property
    def sequence(self):
        return ''.join(str(r.native.type) for r in self)
    
    @property
    def matches(self):
        return self._assignments
    
    @property
    def residues(self):
        return self._residues
    
    def assign(self, fragment):

        if not 1 <= fragment.qstart <= fragment.qend <= len(self._residues):
            raise ValueError("Fragment out of range")
                    
        for rank in range(fragment.qstart, fragment.qend + 1):
            ai = ResidueAssignmentInfo(fragment, rank)
            self._assignments._append_item(fragment)
            self._residues[rank].assignments._append_item(ai)
            
class TargetResidue(object):
    
    def __init__(self, native_residue):
        
        self._type = native_residue.type
        self._native = native_residue.clone()
        self._assignments = csb.pyutils.ReadOnlyCollectionContainer(type=ResidueAssignmentInfo)

    @property
    def type(self):
        return self._type
            
    @property
    def native(self):
        return self._native
    
    @property
    def assignments(self):
        return self._assignments
    
class ResidueAssignmentInfo(object):
    
    def __init__(self, assignment, rank):
        
        if not assignment.qstart <= rank <= assignment.qend:
            raise ValueError('Rank {0} is not matched by this assignment')
        
        self._assignment = assignment
        self._rank = rank
        self._relrank = rank - assignment.qstart 
    
    @property
    def c_alpha(self):
        return self._assignment.backbone[self._relrank]
    
    @property
    def fragment(self):
        return self._assignment    

class Assignment(FragmentMatch):
    
    def __init__(self, source, start, end, id, qstart, qend, probability, rmsd, tm_score):

        sub = source.subregion(start, end, clone=True)
        calpha = [r.atoms['CA'].vector.copy() for r in sub.residues]
                    
        self._calpha = csb.pyutils.ReadOnlyCollectionContainer(items=calpha, type=numpy.ndarray)     
        
        self._source_id = source.entry_id
        self._start = start
        self._end = end
        
        super(Assignment, self).__init__(id, qstart, qend, probability, rmsd, tm_score, None)

    @property
    def backbone(self):
        return self._calpha
        
    @property
    def source_id(self):
        return self._source_id
        
    @property
    def start(self):
        return self._start

    @property
    def end(self):
        return self._end
    
    def transform(self, rotation, translation):
        
        for ca in self.backbone:
            newca = numpy.dot(ca, numpy.transpose(rotation)) + translation
            for i in range(3):
                ca[i] = newca[i]
    
class FragmentTypes(object):
    
    ISites = 'IS'
    HMMFragments = 'HH'
    HHThread = 'TH'
    Rosetta = 'NN'        
    
class BenchmarkAdapter(object):    

    def __init__(self, pdb_path, connection_string=None):
                
        self._pdb = pdb_path
        self._connection = None
        
        from csb.bio.io.wwpdb import find, StructureParser
        self._parser = StructureParser
        self._find = find     
    
        try:    
            import psycopg2.extras
        except ImportError:
            raise RuntimeError('Please install the psycopg2 module first')
        
        if connection_string is None:
            connection_string = self.connection_string()
            
        self._connection = psycopg2.extras.DictConnection(connection_string)
    
    def __del__(self):
        try:
            self.close()          
        except:
            pass
        
    @staticmethod
    def connection_string(database='FragmentBenchmarks', host='', username='', password=''):

        fields = ['dbname={0}'.format(database)]
        
        if host:
            fields.append('host={0}'.format(host))
        if username:
            fields.append('user={0}'.format(username))
            fields.append('password={0}'.format(password))

        return ' '.join(fields)
        
    def close(self):
        self._connection.close()
        
    def targets(self, benchmark_id):
        
        cmd = self._connection.cursor()
        try:
            cmd.callproc('reporting."GetTargets"', (benchmark_id,))
            return cmd.fetchall()
        finally:  
            cmd.close()
            
    def target_details(self, target_id):
        
        cmd = self._connection.cursor()
        try:
            cmd.callproc('reporting."GetTargetDetails"', (target_id,))
            return cmd.fetchall()
        finally:  
            cmd.close()            
            
    def assignments(self, target_id, type):
        
        cmd = self._connection.cursor()
        try:
            cmd.callproc('reporting."GetAssignments"', (target_id, type))
            return cmd.fetchall()
        finally:  
            cmd.close()
            
    def structure(self, accession):

        pdbfile = self._find(accession, self._pdb)
        
        if not pdbfile:
            raise IOError('{0} not found here: {1}'.format(accession, self._pdb))
        
        return self._parser(pdbfile).parse_structure()
            
    def prediction(self, target_id, type):
        
        info = self.target_details(target_id)
        if not info:
            raise ValueError('No such Target ID in the database: {0}'.format(target_id))
        row = info[0]
        
        id = row["Accession"]
        length = float(row["Length"])
        overlap = float(row["MaxOverlap"]) / (length or 1.)
        
        native = self.structure(id[:4]).chains[id[4]]
        target = Target(id, length, native.residues, overlap)
        
        for row in self.assignments(target_id, type):
            
            src_accession = row['Source'][:4]
            src_chain = row['Source'][4]
                        
            frag_chain = self.structure(src_accession).chains[src_chain]
            
            fragment = Assignment(source=frag_chain, 
                                  start=row['SourceStart'], 
                                  end=row['SourceEnd'], 
                                  id=row['FragmentName'], 
                                  qstart=row['Start'], 
                                  qend=row['End'], 
                                  probability=row['Probability'], 
                                  rmsd=row['RMSD'], 
                                  tm_score=row['TMScore'])
            
            target.assign(fragment)
        
        return target
            
    
  
if __name__ == '__main__':
    
    fb = BenchmarkAdapter('/media/DB/pdb/all')
    target = fb.prediction(1, FragmentTypes.HMMFragments)
    print target.id
    
            
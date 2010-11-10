"""
HMM fragments derived from HHpred HMM profiles. 
"""
import re
import csb.bio.fragments
import csb.bio.hmm
import csb.bio.structure as structure
import csb.bio.fragments.isites as isites

from csb.bio.hmm import ProfileHMMSegment, ProfileHMMRegion

class HMMFragment(ProfileHMMSegment):
    """
    Describes a HMM segment which can be slid over a subject profile.
    See the documentation of the base class for the signature of the 
    constructor.
    """
    
    def slide_over(self, other):
        """
        Find instances of self by sliding it over other and computing
        emission vector similarity and RSMD.
        
        @param other: the subject HMM
        @type other: L{ProfileHMM}
        
        @return: a list of L{isites.ProfileMatch}-es
        @rtype: list
        """
        return self._slide_over(self, other)  
    
    def _slide_over(self, this, that):

        query_residues = this.chain()        
        window = this.layers.length
        matches = []
        
        i = that.layers.start_index
        while True:
            if that.layers.last_index - i + 1 < window:
                break
            
            score, tm, tm_out, rmsd_ca = None, None, None, None
            start = i
            end = i + window - 1
            
            subject = ProfileHMMRegion(that, start, end)
            subject_residues = structure.Chain(that.id[4].upper(), residues=that.residues[start : start + window])  

            score = this.emission_similarity(subject)         
                        
            try:
                rmsd_ca = query_residues.rmsd(subject_residues)
            except structure.Broken3DStructureError:
                rmsd_ca = None                          
                                    
            if score is not None and (rmsd_ca is not None):
                match = isites.ProfileMatch(that.id, start, score, rmsd_ca, None, tm, tm_out)
                matches.append(match)               
            
            i += 1
                        
        return matches       
    
class HMMFragmentWrapper(HMMFragment):
    """
    Describes a HMM segment which can be slid over a subject profile.
    Wraps an existing L{ProfileHMMSegment} to decorate it with a slide_over 
    method.
    
    @param hmm_segment: the HMM segment to wrap as a fragment
    @type hmm_segment: L{ProfileHMMSegment}
    """
    
    def __init__(self, hmm_segment):
        
        self._segment = hmm_segment
        
    def slide_over(self, other):

        return self._slide_over(self._segment, other)      
    
class HMMLibraryWrapper(isites.Library):
    """
    Library wrapper for H-Sites fragments.
    
    @param fragments: a list of L{HMMFragment}s
    @type list:  
    """
    
    def __init__(self, fragments=None):
        
        self._ondemand = False
        
        self._fragments = []
        self._byid = {}
        self._byrep = {} 
        self.name = 'HSites'
        
        if fragments is not None:
            self.fragments = fragments        
        
    @property
    def fragments(self):
        return self._fragments
    @fragments.setter
    def fragments(self, fragments):
        self._fragments = []
        self._byid = {}
        self._byrep = {}    
        for fragment in fragments:
            self._fragments.append(fragment)
            self._byid[fragment.isite] = fragment
            self._byrep[fragment.id] = fragment            
    
    @property
    def clusters(self):
        return self.fragments
    @clusters.setter
    def clusters(self, clusters):
        self.fragments = clusters    
        
class HMMFragmentMatch(csb.bio.fragments.FragmentMatch):

    def __init__(self, fragment, qstart, qend, probability, rmsd, tm_score, qlength):
        
        if not isinstance(fragment, csb.bio.hmm.ProfileHMMSegment):
            raise TypeError(fragment)
        
        source_info = re.split('[\-\.]', fragment.id)
        self._source = source_info[0]
        self._start = int(source_info[1]) + fragment.source_start - 1
        self._end = int(source_info[1]) + fragment.source_end - 1
        
        assert (self._end - self._start + 1) == fragment.layers.length == (qend - qstart + 1)     

        super(HMMFragmentMatch, self).__init__(fragment.id, qstart, qend, probability, 
                                                                rmsd, tm_score, qlength)
        
    @property
    def source_id(self):
        return self._source
        
    @property
    def start(self):
        return self._start

    @property
    def end(self):
        return self._end  

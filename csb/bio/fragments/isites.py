"""
I-Sites fragment library APIs. 

@warning: This module is no longer developed and provided as is. 
"""

import sys
import collections
import operator
import cPickle
import math 

import csb.pyutils
import csb.bio.fragments
import csb.bio.structure as structure
import csb.bio.sequence
import csb.bio.hmm as hmm
 
    
class FragmentMatching(csb.pyutils.enum):
    """
    Enumeration of fragment assignment methods
    """
    Sequence=0; Bystroff=1; CrossEntropy=2; Soeding=3; Baker=4

class InternalISitesError(ValueError):
    pass
class InvalidParadigmError(InternalISitesError):
    pass
          
class ChainSequence(object):
    """
    Represents a single PDB chain sequence entry with torsion angles
    """  
    
    def __init__(self):
        self.id = None
        self.accession = None
        self.chain = None
        self.type = None        
        self.header = None
        self.sequence = None
        self.torsion = structure.TorsionAnglesCollection()        
        
    @property
    def has_torsion(self):
        if hasattr(self, 'torsion'):
            return len(self.torsion) > 0
        else:
            return False   
    
    @property
    def entry_id(self):
        id = self.accession
        if self.id is not None:
            id += self.id
        return id
        
    @property
    def length(self):
        if self.sequence is None:
            return 0
        return len(self.sequence)
        
class ProteinProfile(object):
    """
    Describes an I-Sites protein profile/PSSM.
        
    @param background: background amino acid frequencies
    @type background: dict
    @param matrix: initialization array with target frequencies
    @type matrix: dict 
    @param alpha: fixed "pseudocount" number for this cluster/profile
    @type alpha: float
    """      
    
    BackgroundFreqs = { 
                        'A': 0.077049999999999993, 'C': 0.016140000000000002, 'D': 0.058470000000000001, 
                        'E': 0.065640000000000004, 'F': 0.041090000000000002, 'G': 0.072179999999999994, 
                        'H': 0.023630000000000002, 'I': 0.058400000000000001, 'L': 0.089560000000000001, 
                        'K': 0.061460000000000001, 'M': 0.021680000000000001, 'N': 0.045400000000000003, 
                        'P': 0.045179999999999998, 'Q': 0.037080000000000002, 'R': 0.049790000000000001, 
                        'S': 0.061870000000000001, 'T': 0.055660000000000001, 'V': 0.069790000000000005, 
                        'W': 0.014319999999999999, 'Y': 0.035610000000000003, 
                        
                        'B': 0.051935000000000002, 'Z': 0.051360000000000003, 'X': 1.000000000000000000,
                         
                        'O': sys.float_info.min, 'U': sys.float_info.min 
                      } 
    
    def __init__(self, background=BackgroundFreqs, matrix=None, alpha=0.2):
        
        self._matrix = [ ]
        self._pssm = [ ]
        self._maxscore = 0
        self.background = background
        self.alpha = float(alpha)
         
        if matrix:        
            for col in matrix:
                self.add_column(**col)


    def add_column(self, **target_frequencies):
        """
        Append a new column to the profile.
        
        Keyword arguments are used to pass the individual amino acid target frequencies,
        e.g. A=0.12, B=0.0, etc. If an amino acid is omitted, its probability will be 
        set automatically to 0. 
        
        @param target_frequencies: amino acid frequencies in that column
        @type target_frequencies: dict
        
        @return: the index of the new column
        @rtype: int  
        
        @raise KeyError: if the target frequencies dict contains an unknown amino acid 
                         with respect to C{profile.background}          
        """
        if not set(target_frequencies.keys()).issubset(self.background):
            raise KeyError('Unrecognized residue name(s)')        
        
        profcol = { }
        pssmcol = { }
        
        for aa in self.background:
    
            assert self.background[aa] > 0
            
            if aa in target_frequencies:
                assert target_frequencies[aa] >= 0                
                profcol[aa] = target_frequencies[aa]                 
            else:
                profcol[aa] = 0.0
            
            Qa = self.background[aa]
            Pja = profcol[aa]
            alpha = self.alpha
            
            pssmcol[aa] = math.log( (Pja + alpha * Qa) / (Qa * (1 + alpha)) )
            
        self._matrix.append(profcol)
        self._pssm.append(pssmcol)
        self._maxscore += max(pssmcol.values())
        
        return len(self._matrix) - 1
        
    @property
    def length(self):
        return len(self._matrix)
    
    @property
    def matrix(self):
        return self._matrix        

    @property
    def pssm(self):
        return self._pssm
    
    @property
    def max_score(self):
        return self._maxscore

class ISitesFragmentMatch(csb.bio.fragments.FragmentMatch):
    
    def __init__(self, fragment, qstart, qend, probability, rmsd, tm_score, qlength):
        
        self._source = fragment.representative.accession + fragment.representative.chain
        self._start = fragment.representative.structure.residues[1].rank
        self._end = fragment.representative.structure.residues[-1].rank
        assert (self._end - self._start + 1) == fragment.length == (qend - qstart + 1)     

        super(ISitesFragmentMatch, self).__init__(fragment.id, qstart, qend, probability, 
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
    
class ProfileMatch(object):
    """
    Describes a profile-query match
    
    @param id: id of the match
    @type id: str
    @param start: start position of the match
    @type start: int
    @param score: score
    @type score: float
    @param rmsd: RMSD between the query and the representative structure of the 
                 fragment withing the matching region
    @type rmsd: float
    @param rmsd_dih: circular (torsion) RMSD between the query and the representative 
                    structure of the fragment withing the matching region
    @type rmsd_dih: float
    @param tm_score: TM-Score between the query and the representative 
                    structure of the fragment withing the matching region
    @type tm_score: float    
    @param tm_scoreout: TM-Score, but calculated with local structural alignment
    @type tm_scoreout: float   
    """ 
        
    def __init__(self, id, start, score, rmsd=None, rmsd_dih=None, tm_score=None, tm_scoreout=None):        
        self.id = id
        self.start = start
        self.score = score
        self.rmsd = rmsd
        self.rmsd_dih = rmsd_dih
        self.tm_score = tm_score
        self.tm_scoreout = tm_scoreout
    
    def __repr__(self):
        return "<ProfileMatch: {0.id} at {0.start}, s={0.score}, rms={0.rmsd}, tm={0.tm_score},{0.tm_scoreout} >".format(self)
    
    def __cmp__(self, other):
        return cmp(self.score, other.score)
    
    def confidence(self, b0, b1, b2, b3=1.0):
        """
        Convert the raw profile score to probability.
        
        @param b0: C{cluster.linearfit[0]}
        @type b0: float
        @param b1: C{cluster.linearfit[1]}
        @type b1: float
        @param b2: C{cluster.linearfit[2]}
        @type b2: float
        @param b3: C{cluster.linearfit[3]}
        @type b3: float       
        
        @return: probability for the hit being a true positive
        @rtype: float                 
        """
        
        if self.score is None:
            return 0
        return ((1 / math.pi) * math.atan(b0 + b1 * self.score) + b2) * b3
   
class ProfileMatches(object):
    """
    A collection of profile-database hits.
    
    @param size: target maximum size of the hitlist
    @type size: int
    @param cutoff: inclusion score cutoff
    @type cutoff: float
    @param rmsd_cutoff: inclusion RMSD cutoff
    @type rmsd_cutoff: float
    
    @raise ValueError: on invalid hitlist C{size}
    """
    
    def __init__(self, size=20, cutoff=0.5, rmsd_cutoff=0.5):
        if not size > 0:
            raise ValueError('The maximum size of the matches table must be a positive number.')
                
        self._matches = [ ]
        self.list_size = int(size)
        self.cutoff = float(cutoff)
        self.rmsd_cutoff = float(rmsd_cutoff)
        
        self._sorted = True        
        
    def _sort(self):
        if not self._sorted:
            self._matches.sort(key=operator.attrgetter('score'), reverse=True)
            self._sorted = True
        
    def add(self, id, start, score, rmsd=None):    
        """ 
        Add a new match to the collection of hits if its score satisfies the cutoff.
        
        @param id: id of the match
        @type id: str
        @param start: start position of the match
        @type start: int
        @param score: score
        @type score: float
        @param rmsd: RMSD between the query and the representative structure of the 
                     fragment withing the matching region
        @type rmsd: float
        
        @raise value_error: on invalid start position        
        """
        if not start > 0:
            raise ValueError('Match start position must be greater than zero.')
        score = float(score)
        start = int(start)
        
        self._sorted = False
                                
        if len(self._matches) < self.list_size: 
            if score >= self.cutoff: #or rmsd <= self.rmsd_cutoff:
                match = ProfileMatch(id, start, score, rmsd)                
                self._matches.append(match)
                # sort is useless if not all slots are full, leave it for the first time the data is accessed through self.hits              
        else:
            best_match = self._matches[-1]
            if score > best_match.score: #or rmsd < best_match.rmsd:                
                self._matches.pop()
                match = ProfileMatch(id, start, score, rmsd)                
                self._matches.append(match)
                self._sort()   # sort is needed in real time if all hit slots are already occupied 

    @property
    def length(self):
        return len(self._matches)
                
    @property
    def hits(self):
        if self.length > 0:
            self._sort()
            return self._matches
        else:
            return None
    
    @property
    def best_match(self):
        if self.length > 0:
            return self.hits[0]
        else:
            return None

class Hitlist(object):
    """
    I-Site instance hitlist.
    
    @param query: the query fragment
    @type query: L{Cluster}
    @param hits: a list of L{ProfileMatch}es
    @type hits: list
    """

    def __init__(self, query, hits):
        
        self.query = query
        self.hits = hits

    def __getitem__(self, i):
        return self.hits[i]

    def __iter__(self):
        return iter(self.hits)    
    
    def best_hit(self):
        """
        @return: the best hit in the list in terms of profile score
        @rtype: L{ProfileMatch} 
        """
        best = None
        for hit in self.hits:
            if best is None:
                best = hit
            else:
                if best.score < hit.score or (best.score == hit.score and best.rmsd > hit.rmsd):
                    best = hit
        return best
    
    def top(self, n):
        """
        @param n: number of top hits
        @type n: int
        
        @return: top C{n} hits in terms of profile score
        @rtype: list of L{ProfileMatch} 
        """
        
        hits = list(self.hits)
        hits.sort()
        return hits[-n:]
                       
class Library(object):
    """
    Representation of an I-Sites peptide fragment library.
    Provides dictionary-like access to the underlying L{Cluster} objects by 
    either C{cluster.id} or C{cluster.representative.accesion}.
    """
        
    def __init__(self):

        self.name = 'ISITES'        
        self.version = None
        self.centroids = None
        self.documentation = None
                    
        self._clusters = None           
        self._ondemand = False
        self._byid = None
        self._byrep = None
        
        self.__iscfiles = []
        
    def __getitem__(self, item):
        if self._ondemand:
            raise AttributeError("ID/rep-based access to clusters is not available in Delay-parsed libraries.")
        
        if isinstance(item, basestring):
            return self._byrep[item]
        else:
            return self._byid[item]
        
    def __contains__(self, key):
        try:
            self[key]
            return True
        except KeyError:
            return False            
        
    @property
    def clusters(self):
        return self._clusters
    
    @clusters.setter
    def clusters(self, clusters):
        if clusters is not None:
            if not isinstance(clusters, collections.Iterable):
                clusters = [clusters]
#                if not isinstance(clusters, collections.Iterator):
#                    for c in clusters:
#                        assert isinstance(c, Cluster), "The I-Sites Library is constructed of Cluster objects."
            
        self._clusters = clusters
        if isinstance(clusters, collections.Iterator):
            self._ondemand = True
        else:
            self._ondemand = False
            self._byid = { }
            self._byrep = { }

            for c in clusters:
                self._byid[c.id] = c
                
                accession = c.representative.accession[:4]
                
                if accession not in self._byrep:
                    self._byrep[accession] = [ ]
                self._byrep[accession].append(c)
                
    def compact(self):
        """
        Compact the library by excluding unimportant information:
        
            - documentation
            - clusters: file, program, covariance
            
        @raise AttributeError: on attempt to compact a delay-parsed library            
        """
        if self._ondemand:
            raise AttributeError("Delay-parsed libraries cannot be compacted.")

        self.documentation = None        
        for c in self.clusters:
            c.file = None
            c.program = None
            c.covariance = None

    def serialize(self, file):
        """
        Serialize the library to a file.
        
        @param file: output file name
        @type file: str
        
        @raise AttributeError: on attempt to serialize a delay-parsed library
        """
        if self._ondemand:
            raise AttributeError("Delay-parsed libraries cannot be serialized.")
                
        with open(file, 'wb') as output: 
            cPickle.dump(self, output, protocol=cPickle.HIGHEST_PROTOCOL)
    
    @staticmethod
    def deserialize(file):
        """
        Load I-sites from a serialization file.
        
        @param file: source file name
        @type file: str
        
        @return: the deserialized library object
        @rtype: L{Library}
        """

        with open(file, 'rb') as output: 
            library = cPickle.load(output)        
            
        return library

class Cluster(object):
    """
    Object representation of an I-Sites library peptide fragment.
    """    
       
    def __init__(self):
        
        self.id = None
        self.motiflen = None
        self.profilelen = None
        self.overhang = None
        self.file = None
        self.program = None
        self.mda = None
        self.dme = None            
        self.pseudocount = 0.2
        self.linearfit = None
        self.representative = None
        self.covarweight = None
        self.profile = ProteinProfile(ProteinProfile.BackgroundFreqs, alpha=self.pseudocount) 
        self.covariance = None

    def __cmp__(self, other):
        return cmp(self.id, other.id) 

    @property
    def isite(self):
        return self.id
    
    @property
    def length(self):
        return self.profile.length
        
    @property
    def has_structure(self):
        if self.representative.structure:
            return True
        return False
        
    def serialize(self, file):
        """
        Serialize the cluster to a file.
        
        @param file: output file name
        @type file: str        
        """
                
        with open(file, 'wb') as output: 
            cPickle.dump(self, output, protocol=cPickle.HIGHEST_PROTOCOL)
    
    @staticmethod
    def deserialize(file):
        """
        Load a cluster from a serialized object.
        
        @param file: source file name
        @type file: str
        
        @return: deserialized fragment object
        @rtype: L{Cluster}
        """

        with open(file, 'rb') as output: 
            cluster = cPickle.load(output)        
            
        return cluster            
        
    def attach_structure(self, pdb_file):
        """
        Parse the paradigm's pdb_file and attach structural data to 
        C{cluster.representative} (C{source} and C{structure} attributes).
        
        @param pdb_file: source PDB file containing the structure of the paradigm
        @type pdb_file: str
        
        @raise InternalISitesError: on parse and/or structure assignment errors
        @raise InvalidParadigmError: when the torsion angles in the cluster do not
                                     match the angles computed from the PDB file
        """
        from csb.bio.io import StructureParser
        rep = self.representative
        s = StructureParser(pdb_file).parse_structure()
        
        if not rep.chain:
            if s.chains.length > 1:
                raise InternalISitesError('Cluster {0} does not define any chain, but multiple choices exist.'.format(self.id))
            else:
                rep.chain = s.chains.keys()[0]
        if rep.chain not in s.chains:
            raise InternalISitesError('Structure {0} does not contain chain {1}'.format(rep.accession, rep.chain))
        
        s.chains[rep.chain].compute_torsion()
        start_residue = s.chains[rep.chain].find(rep.start)
           
        current_rank = start_residue.rank
        for i, rep_angles in enumerate(rep.angles, start=rep.angles.start_index):
            
            current_residue = s.chains[rep.chain].residues[current_rank]
            torsion_delta = 0
            
            if current_residue.torsion.phi:
                torsion_delta += abs(current_residue.torsion.phi - rep_angles.phi) 
            if current_residue.torsion.psi:
                torsion_delta += abs(current_residue.torsion.psi - rep_angles.psi)
            
            assert (i + rep.start - 1) == int(current_residue.id)
            
            if torsion_delta > 1:
                raise InvalidParadigmError('delta_phi + delta_psi for '
                                          '{3.accession} {3.chain}[{4}] is {0}: phi: {1.phi} vs {2.phi}, psi: {1.psi} vs {2.psi}'.format(
                                                                torsion_delta, current_residue.torsion, rep_angles, rep, current_rank))            
            current_rank += 1
        
        start = start_residue.rank - self.overhang
        stop = start + self.profile.length - 1
        
        if start < s.chains[rep.chain].residues.start_index or stop > s.chains[rep.chain].residues.last_index:
            raise InternalISitesError('With the ovehangs included, fragment '
                                      '{0} exceeds the boundaries of the chain: {1}..{2} vs {3.start_index}..{3.last_index}'.format(
                                                                                self.id, start, stop,  s.chains[rep.chain].residues))
        
        rep.structure = s.chains[rep.chain].subregion(start, stop)

        assert rep.structure.length == self.profile.length
                
        backbone = set(['N', 'CA', 'C'])
        for residue in rep.structure.residues:
            if not (residue.has_structure and set(residue.atoms).issuperset(backbone)):
                rep.structure = None
                raise InternalISitesError('Fragment {0} is built on discontinuous structure: {1}  {4} at {2}-{3}'.format(
                                                                                self.id, rep.accession, start, stop, residue))
        
        rep.source = ChainSequence()
        rep.source.sequence = s.chains[rep.chain].sequence
        rep.source.accession = s.accession
        rep.source.id = rep.chain
        rep.source.type = s.chains[rep.chain].type
        rep.source.torsion = s.chains[rep.chain].torsion
        
    def sequence_similarity(self, sequence, start=0):
        """
        Compute the log-odds score of C{sequence} matching C{self}'s profile.
        
        @param sequence: subject sequence
        @type sequence: str
        @param start: do the comparison starting with this offset in the sequence
        @type start: int
        
        @return: log-odds score
        @rtype: float
        """
        score = 0
        window = self.profile.length
        
        sequence = sequence[start : start + window]         
        assert len(sequence) ==  self.profile.length
        
        for column, aa in enumerate(sequence):
            aa = aa.upper()
            if aa in self.profile.pssm[column]:
                score += self.profile.pssm[column][aa]
            else:
                raise ValueError('Unknown residue {0}'.format(aa))     
        
        return score   
        
    def profile_similarity(self, profile, start=0):
        """
        Compute the similarity score of C{profile} matching C{self}'s profile (see
        Bystroff 2001, 2008).
        
        @param profile: subject profile, built with alpha=0.5
        @type profile: L{ProteinProfile}
        @param start: do the comparison starting with this offset in the subject
        @type start: int
        
        @return: similarity score calculated with Bystroff's metric
        @rtype: float
        
        @raise ValueError: on mismatch in the amino acid content at a given column 
                           in both profiles; or when the subject has not been created
                           with alpha=0.5
        """        
        if profile.alpha != 0.5:
            raise ValueError('The target profile must be built with alpha=0.5')

        score = 0
        window = self.profile.length
                
        for column in range(0, window): 
            
            if set(self.profile.pssm[column]) != set(profile.pssm[column + start]):
                raise ValueError('Both profiles must contain identical amino acids at a given column')
                        
            for aa in self.profile.pssm[column]:
                score += self.profile.pssm[column][aa] * profile.pssm[column + start][aa]
        
        return score
    
    def cross_entropy(self, profile, start=0):
        """
        Compute the similarity score of C{profile} matching C{self}'s profile using
        the cross entropy method.
        
        @param profile: subject profile
        @type profile: L{ProteinProfile}
        @param start: do the comparison starting with this offset in the subject
        @type start: int
        
        @return: the cross entropy between the profiles 
        @rtype: float
        """        
        score = 0
        window = self.profile.length
                
        for column in range(0, window):
            for aa in self.profile.pssm[column]:
                score += self.profile.pssm[column][aa] * math.exp( profile.pssm[column + start][aa] )
        
        return score    
    
    def sum_of_odds(self, profile, start=0):
        """
        Compute the similarity score of C{profile} matching C{self}'s profile using
        the log-sum-of-odds method (Soeding).
        
        @param profile: subject profile
        @type profile: L{ProteinProfile}
        @param start: do the comparison starting with this offset in the subject
        @type start: int
        
        @return: the log-sum-of-odds similarity between the profiles 
        @rtype: float
        """        
        score = 1
        window = self.profile.length
                
        for column in range(0, window):
            
            dot_product = 0
            
            for aa in self.profile.pssm[column]:
                q_emission = self.profile.matrix[column][aa]
                s_emission = profile.matrix[column + start][aa]
                dot_product += ((q_emission * s_emission) / self.profile.background[aa]) 
            
            score *= dot_product
            
        return math.log(score or sys.float_info.min)      
    
    def city_block(self, profile, start=0):
        """
        Compute the similarity score of C{profile} matching C{self}'s profile using
        the exponential city block metric (see Baker 1995).

        @param profile: subject profile
        @type profile: L{ProteinProfile}
        @param start: do the comparison starting with this offset in the subject
        @type start: int
        
        @return: the city block similarity between the profiles 
        @rtype: float
        """               

        score = 0
        window = self.profile.length
                
        for column in range(0, window):
            for aa in self.profile.pssm[column]:
                diff = math.exp( -self.profile.matrix[column][aa] ) - math.exp( -profile.matrix[column + start][aa] )
                score += abs(diff)
        
        return score    
        
    def fragment_comparer(self, mode=FragmentMatching.Bystroff):
        """
        Return a callable fragment comparer that implements the specified 
        comparison mode. 
        
        @param mode: specifies the desired comparison method to be implemented
                     by the callable - a member of the L{FragmentMatching} enum
        @type mode: L{csb.pyutils.EnumItem}
        
        @return: the proper comparison function which implements fragment comparison
                 with the C{mode} method
        @rtype: function
        
        @raise ValueError: on invalid C{mode} specified
        """
        if mode == FragmentMatching.Bystroff:
            return self.profile_similarity
        elif mode == FragmentMatching.CrossEntropy:
            return self.cross_entropy
        elif mode == FragmentMatching.Soeding:
            return self.sum_of_odds
        elif mode == FragmentMatching.Sequence:
            return self.sequence_similarity
        elif mode == FragmentMatching.Baker:
            return self.city_block
        else:
            raise ValueError(mode)     
                
    def slide_over(self, target, mode=FragmentMatching.Bystroff):
        """
        Find profile and structural matches by sliding the fragment over a specified 
        C{target} chain.
        
        @param target: subject chain to be scanned to fragment instances
        @type target: L{structure.Chain}
        @param mode: fragment instance searching mode - a L{FragmentMatching} member
        @type mode: L{csb.pyutils.EnumItem}
        
        @return: a list of L{ProfileMatch}es
        @rtype: list
        """      
        
        compare = self.fragment_comparer(mode)
        
        if hasattr(target, 'entry_id'):
            entry_id = target.entry_id
        else:
            entry_id = target.id
       
        window = self.profile.length
        rmsd_window = window - 2 * self.overhang        
        matches = [ ]
        query_residues = structure.Chain(self.representative.chain, 
                                        residues=self.representative.structure.residues[self.overhang + 1: -self.overhang])
        target_sequence = target.sequence        # caching these values for efficiency
        
        for i in range(0, len(target_sequence) - window + 1):
            
            start = i + 1
            rmsd_start = start + self.overhang                          
            score, tm, tm_out, rmsd_ca, rmsd_dih = 0, None, None, None, None
            
            if mode == FragmentMatching.Sequence:
                score = compare(target_sequence, start=i)
            else:
                score = compare(target.profile, start=i)
                                
            subject_residues = structure.Chain(target.id, residues=target.residues[rmsd_start : rmsd_start + rmsd_window])  
                        
            try:
                rmsd_ca = query_residues.rmsd(subject_residues)
            except structure.Broken3DStructureError:
                rmsd_ca = None                
                pass                  
                                    
            if score is not None and (rmsd_ca is not None or rmsd_dih is not None):
                match = ProfileMatch(entry_id, start, score, rmsd_ca, rmsd_dih, tm, tm_out)
                matches.append(match)   
        
        return matches
    
    def slide_over_fast(self, target, mode=FragmentMatching.Bystroff):
        """
        Find profile matches only (do not compute RMSD) by sliding the fragment
        over a specified C{target} chain.

        @param target: subject chain to be scanned to fragment instances
        @type target: L{structure.Chain}
        @param mode: fragment instance searching mode - a L{FragmentMatching} member
        @type mode: L{csb.pyutils.EnumItem}
        
        @return: a list of L{ProfileMatch}es
        @rtype: list
        """      
        if hasattr(target, 'entry_id'):
            entry_id = target.entry_id
        else:
            entry_id = target.id
                
        compare = self.fragment_comparer(mode)
       
        window = self.profile.length
        matches = [ ]
        
        for i in range(0, target.profile.length - window + 1):
            
            start = i + 1                      
            score = None
            
            if mode == FragmentMatching.Sequence:
                score = compare(target.sequence, start=i)
            else:
                score = compare(target.profile, start=i)               
                                    
            if score is not None:
                match = ProfileMatch(entry_id, start, score)
                matches.append(match)   
        
        return matches
    
    def to_hmm(self):
        """
        Convert this fragment into a profile HMM segment. Each column in the
        fragment's profile becomes a Match state with equivalent emissions.
        However, all Match states are connected with a fixed transition 
        probability of 100%.      
        
        @return: an HMM wrapper around the cluster's profile
        @rtype: L{hmm.ProfileHMM}
        """
        
        factory = hmm.StateFactory()

        p = hmm.ProfileHMM(hmm.ScoreUnits.Probability)
        p.family = 'is' + str(self.id)
        p.name = 'is' + str(self.id)
        p.id = str(self.id)   
        p.pseudocounts = False
        p.version = '1.6'
        
        background = dict([ (aa, self.profile.background[aa]) 
                           for aa in 'ACDEFGHIKLMNPQRSTVWY' ]) 
        
        for i, row in enumerate(self.profile.matrix, start=1):
            
            row = dict([ (aa, row[aa]) for aa in 'ACDEFGHIKLMNPQRSTVWY' ]) 
            residue = self.representative.structure.residues[i]
            
            match = factory.create_match(row, background)
            match.rank = i
            
            layer = hmm.HMMLayer(i, residue)
            layer.append(match)
            
            p.layers.append(layer)
    
        for layer in p.layers:
            
            match = layer[hmm.States.Match]
            
            if layer.rank < p.layers.last_index:
                next_match = p.layers[layer.rank + 1][hmm.States.Match]
                tran = hmm.Transition(match, next_match, 1.0)
                match.transitions.append(tran)
    
        last_tran = hmm.Transition(p.layers[-1][hmm.States.Match], p.end, 1.0)
        p.layers[-1][hmm.States.Match].transitions.append(last_tran)
        
        start_tran = hmm.Transition(p.start, p.layers[1][hmm.States.Match], 1.0)
        p.start.transitions.append(start_tran)        
        
        p.length.matches = p.layers.length
        p.length.layers = p.layers.length
        p.effective_matches = 10
        
        rep_sequence = ''.join([ str(l.residue.type) for l in p.layers ])
        seq = csb.bio.sequence.Sequence('dummy', '>dummy', rep_sequence,
                            type=csb.bio.sequence.SequenceTypes.Protein)
        p.alignment = csb.bio.sequence.A3MAlignment([seq])
    
        return p        

    @staticmethod
    def from_hmm(source, id, start=None, end=None, rep_accession=None,
                 rep_chain=None, alpha=0.2):
        """
        Build a fragment from a Profile HMM.
        
        @param source: the HMM to convert
        @type source: L{hmm.ProfileHMM}, L{hmm.ProfileHMMSegment} 
        @param id: ID of the new fragment
        @type id: int
        @param start: start position (optional)
        @type start: int
        @param end: end position (optional)
        @type end: int
        @param rep_accession: accession number for I{cluster.representative}
        @type rep_accession: str
        @param rep_chain: chain ID for I{cluster.representative}
        @type rep_chain: str
        @param alpha: alpha pseudocount for I{cluster.profile.pssm}
        @type alpha: float

        @return: an I-Sites cluster wrapper
        @rtype: L{Cluster}

        @raise TypeError: when the HMM's structural data is missing or
                          interrupted
        """

        if not source.has_structure:
            raise TypeError('This HMM has no structural data assigned and '
                            'cannot be converted.')

        if start is None:
            start = source.layers.start_index
        if end is None:
            end = source.layers.last_index

        score_units = source.score_units
        source.convert_scores(hmm.ScoreUnits.Probability)

        background = dict(ProteinProfile.BackgroundFreqs)
        for aa in source.layers[1][hmm.States.Match].emission:
            background[str(aa)] = source.layers[1][hmm.States.Match].background[aa]

        cluster = Cluster()
        cluster.id = id
        cluster.overhang = 0
        cluster.profile = ProteinProfile(background, alpha=alpha)

        residues = []
        for layer in source.layers:

            residues.append(csb.pyutils.deepcopy(layer.residue))

            if start <= layer.rank <= end:
                if not layer.residue.has_structure:
                    raise TypeError(
                        "The HMM's structure is interrupted at layer {0}.".format(
                            layer.rank))

                emission = {}

                for aa in layer[hmm.States.Match].emission:
                    emission[str(aa)] = layer[hmm.States.Match].emission[aa]

                cluster.profile.add_column(**dict(emission))

        cluster.profilelen = cluster.profile.length
        cluster.motiflen = cluster.profile.length - 2 * cluster.overhang

        if rep_accession is None:
            rep_accession = source.id[1:5].lower()
        if rep_chain is None:
            rep_chain = source.id[5].upper()

        chain = structure.Chain(rep_chain, structure.SequenceTypes.Protein,             
                                None, residues, rep_accession)
        chain.compute_torsion()

        src = ChainSequence()
        src.sequence = chain.sequence
        src.accession = chain.accession
        src.id = chain.id
        src.type = chain.type
        src.torsion = chain.torsion

        cluster.representative = RepStructureFragment(chain.accession,
                                                      chain.id, start)
        cluster.representative.source = src
        cluster.representative.structure = chain.subregion(start, end)
        cluster.representative.angles = cluster.representative.structure.torsion

        assert cluster.representative.angles.length == cluster.motiflen
        assert cluster.profile.length == cluster.representative.structure.length
        assert cluster.representative.angles.length == (cluster.profile.length -
                                                        2 * cluster.overhang)

        source.convert_scores(score_units)
        
        return cluster        

class RepStructureFragment(object):
    """
    Container class which describes an I-Sites paradigm structure for a cluster.
    
    @param accession: paradigm's accession number
    @type accession: str
    @param chain: chain ID
    @type chain: str
    @param start: the start position in the paradigm's chain where the fragment
                  instance is located (including the overhangs). In the original 
                  library this is a PDB sequence number, not a SEQRES rank. The
                  residue rank can be retrieved from the structure instead:
                  
                      >>> cluster.representative.structure.find(start).rank
                      residue rank
                      
    @type start: int
    @param angles: torsion angles of the paradigm, in degrees
    @type angles: L{structure.TorsionAnglesCollection}
    @param source: the entire paradigm chain with pre-computed torsion angles
    @type source: L{ChainSequence}
    @param struct: a segment of the paradigm's L{structure.Chain} which
                      contains the I-Site. Residue ranks and IDs are preserved
                      as they were in the original chain (no shifting)
    @type struct: L{structure.Chain}
    """
    
    def __init__(self, accession, chain, start, angles=None, source=None, struct=None):
        
        self.accession = accession
        self.chain = chain
        self.start = start
        self._source = None
        self._angles = structure.TorsionAnglesCollection()
        self._structure = None
        
        
        if angles is not None:
            self.angles = angles   
        if source is not None:
            self.source = source
        if struct is not None:
            self.structure = struct       
    
    @property
    def angles(self):
        return self._angles
    @angles.setter
    def angles(self, angles):        
        self._angles.update(angles)
           
    @property
    def source(self):
        return self._source
    @source.setter
    def source(self, chain):
        if type(chain) not in (structure.Chain, ChainSequence):
            raise TypeError("The source property must be a Chain or ChainSequence instance with torsion property pre-calculated.")
        self._source = chain
        
    @property
    def structure(self):
        return self._structure
    @structure.setter
    def structure(self, chain_fragment):
        if not isinstance(chain_fragment, structure.Chain):
            raise TypeError("The structure property must be a Chain instance.")
        self._structure = chain_fragment        

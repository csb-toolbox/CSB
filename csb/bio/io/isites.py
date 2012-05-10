"""
I-Sites fragment library parser.

@deprecated: legacy module.
"""

import os
import re
import csb.io
import csb.bio.structure as structure
import csb.bio.fragments.isites as isites

class Tags(object):
    """
    Enumeration of I-Sites flat file tags.
    """
    
    LIBRARY = 'LIBRARY'
    REMARK = 'REMARK'
    CENTROID = 'CENTROID'
    CLUSTER = 'CLUSTER'
    FROM = 'FROM'
    CREATEDBY = 'CREATEDBY'
    MDACUT = 'MDACUT'
    DMECUT = 'DMECUT'
    PSEUDOCOUNT = 'PSEUDOCOUNT'
    LINEARFIT = 'LINEARFIT'
    COVARWEIGHT = 'COVARWEIGHT'
    PARADIGM = 'PARADIGM'
    ANGLES = 'ANGLES'
    PROFILE = 'PROFILE'
    COVARTENSOR = 'COVARTENSOR'
    END = 'END' 
    
class ISitesParser(object):
    """
    Implements an I-Sites fragment library parser v.5.1+ (2008).
    
    @param flatfile: input *.isl I-Sites library file name
    @type flatfile: str 
    @param express: if True, speeds up the parser by ignoring the covariance tensors
    @type express: bool
    
    @raise IOError: when the source file cannot be found
    """
    
    def __init__(self, flatfile, express=False):
        if not os.path.exists(flatfile):
            raise IOError("Could not read file {0}".format(flatfile))
        
        self._flatfile = flatfile
        self.express = bool(express)
        self._streams = [ ]
        self._library = None
        
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        for s in self._streams:
            try:
                s.close()
            except:
                pass   
        
    def __del__(self):
        for s in self._streams:
            try:
                s.close()
            except:
                pass   
        
    def _newstream(self):
        
        s = open(self._flatfile, mode='r')
        self._streams.append(s)
        return s
            
    @property
    def library(self):
        """
        Return the general properties of the library. 
        Library Clusters are iterable, but read efficiently only on demand.
        """
        return self.parse()
        
    @property 
    def clusters(self):
        """
        Efficient iterator over all L{Cluster} objects.
        """     
        reader = csb.io.EntryReader(self._newstream(), Tags.CLUSTER, Tags.END)
        
        for entry in reader.entries():
            yield self.parse_entry(entry)

    def parseall(self):
        """
        Parse the whole library to end. 
        
        @return: object representation of the library with all clusters pre-parsed
        @rtype: L{Library}
        """
        library = self.parse()
        library.clusters = list(self.clusters)
        return library       
    
    def parse(self):
        """ 
        Parse I-sites library common/general properties. Clusters are not parsed, 
        but can be fetched on demand while iterating over C{library.clusters}.
        
        @return: object representation of the library with a 
                bound clusters generator
        @rtype: L{Library}
        """
        
        library = isites.Library()
        library.centroids = [ ]
        library.documentation = ''
        library.clusters = self.clusters
        
        stream = self._newstream()
        done = False
        
        while not done:
            
            line = stream.readline()
            if not line:
                done = True
                break            
            
            if line.startswith(Tags.REMARK) and line[len(Tags.REMARK):].strip() in ('===== start of library =====', '====== start of clusters ======='):
                done = True
                break
            
            elif line.startswith(Tags.LIBRARY):
                fields = line.split()[1:]
                if len(fields) > 1:
                    library.name, library.version = fields[0], fields[1]
                else:
                    library.version = fields[0]
            
            elif line.startswith(Tags.CENTROID):
                fields = line.split()
                index = int(fields[1])
                values = fields[2].split(',')
                 
                matrow = { }
                       
                for cn, aa in enumerate('ACDEFGHIKLMNPQRSTVWY'):
                    matrow[aa] = float(values[cn])
                    
                if index > len(library.centroids) - 1: # This is because in the current version CENTROIDS appears twice by mistake
                    assert index == len(library.centroids), "Centroid offset indices are consecutive numbers, starting from 0."                    
                    library.centroids.append(matrow)

            elif line.startswith(Tags.REMARK):
                library.documentation += line[len(Tags.REMARK)+1:]
                
        stream.close()
        
        return library                
    
    def parse_entry(self, entry):
        """
        Parse a single I-Sites entry. 
        
        @return: object representation of the entry
        @rtype: L{Cluster}
        """
        cluster = isites.Cluster()
        lines = iter(entry.splitlines())
        
        done = False
        in_profile = False
        in_tensor = False
        
        while not done:
            try:
                line = lines.next()
            except StopIteration:
                done = True
                break
            
            if line.startswith(Tags.CLUSTER):
                fields = line.split()[1:]
                cluster.id, cluster.motiflen, cluster.profilelen, cluster.overhang = map(int, fields) 

            elif line.startswith(Tags.FROM):
                cluster.file = line.split()[1]

            elif line.startswith(Tags.CREATEDBY):
                cluster.program = line[len(Tags.CREATEDBY)+1:].strip()
                                                                           
            elif line.startswith(Tags.MDACUT):
                field = line.split()[1]
                cluster.mda = float(field)
                
            elif line.startswith(Tags.DMECUT):
                field = line.split()[1]                
                cluster.dme = float(field)
                
            elif line.startswith(Tags.PSEUDOCOUNT):
                field = line.split()[1]   
                cluster.pseudocount = float(field)
                assert cluster.pseudocount > 0
                
            elif line.startswith(Tags.LINEARFIT):
                fields = line.split()[1:]
                cluster.linearfit = map(float, fields)
                
            elif line.startswith(Tags.COVARWEIGHT):
                field = line.split()[1]                  
                cluster.covarweight = float(field)
                
            elif line.startswith(Tags.PARADIGM):
                fields = line.split()[1:]
                cluster.representative = isites.RepStructureFragment(fields[0], fields[1], int(fields[2]))
                if fields[1] == '_':
                    cluster.representative.chain = ''        

            elif line.startswith(Tags.ANGLES): 
                rn = -1           
                while True:
                    try:
                        subline = lines.next()
                    except StopIteration:
                        break
                    if subline.startswith(Tags.PROFILE):
                        in_profile = True
                        break
                    elif  subline.startswith(Tags.END):
                        break
                    
                    rn += 1
                    fields = subline.split()
                    angles = map(float, fields[1:])                    
                    
                    torsion = structure.TorsionAngles(angles[0], angles[1], angles[2], units=structure.AngleUnits.Degrees)
                    j = cluster.representative.angles.append(torsion)
                    
                    assert rn == j-1 == int(fields[0]), "Angle offsets in a cluster are consecutive numbers, starting at 0."                    

            elif line.startswith(Tags.PROFILE):
                in_profile = True                
            
            elif in_profile:
                cluster.profile = isites.ProteinProfile(isites.ProteinProfile.BackgroundFreqs, alpha=cluster.pseudocount)
                rn = -1     
                subline = line
                                      
                while True: 
                    if subline.startswith(Tags.CREATEDBY) or subline.startswith(Tags.END):
                        in_profile = False
                        break
                    elif subline.startswith(Tags.COVARTENSOR):                        
                        in_tensor = True
                        in_profile = False
                        break
                    
                    rn += 1
                    fields = subline.split()                    
                    
                    assert rn == int(fields[0]), "ProteinProfile rows in a cluster are consecutive numbers," \
                                                 + " starting from 0 (cluster {0}, profile row {1}/{2}).".format(cluster.id, rn, fields[0])
                    column = { }                                 
                    for cn, aa in enumerate('ACDEFGHIKLMNPQRSTVWY'):
                        column[aa] = float(fields[cn+1])                           
                    cluster.profile.add_column(**column)
                    
                    try:
                        subline = lines.next()
                    except StopIteration:
                        in_profile = False
                        break
                
                assert cluster.profilelen == cluster.profile.length 
                    
            elif line.startswith(Tags.COVARTENSOR):
                in_tensor = True

            elif in_tensor:
                if self.express:
                    break
                                
                motiflen = cluster.motiflen 
                # cluster.covariance = [ [ [] ]*motiflen for ii in range(0, motiflen) ]
                cluster.covariance = [ ]
                for mi in range(0, motiflen):
                    cluster.covariance.append([])                    
                    for mj in range(0, motiflen):       #@UnusedVariable
                        cluster.covariance[mi].append([])
                
                rn = -1
                i = j = -1                     
                subline = line
                dimline = re.compile('^[0-9]+\s+[0-9]+\s*$')
                                      
                while True: 
                    if subline.startswith(Tags.END):
                        in_tensor = False
                        break
                    
                    rn += 1
                    fields = subline.split()
                                                            
                    if re.match(dimline, subline):
                        istr, jstr = subline.split()
                        i, j = int(istr) - 1, int(jstr) - 1
                        assert 0 <= i < motiflen and 0 <= j < motiflen, "Covariance is a [motiflen x motiflen] matrix."
                    else:
                        values = map(float, subline.split())
                        cluster.covariance[i][j].append(values)
                    
                    try:
                        subline = lines.next()
                    except StopIteration:
                        in_tensor = False                        
                        break                    
            
            elif line.startswith(Tags.END):
                done = True
                break
                                
        return cluster

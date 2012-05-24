__version__ = '0.0.1.{revision}'


class Version(object):
    """
    CSB version number.
    """
    
    def __init__(self):
        
        version = __version__.split('.')
        
        if not len(version) in (3, 4):
            raise ValueError(version)
        
        self._major = version[0]
        self._minor = version[1]
        self._micro = version[2]
        self._revision = None
        
        if len(version) == 4:
            self._revision = version[3]

    def __str__(self):
        return self.short
    
    @property
    def major(self):
        return int(self._major)  
     
    @property
    def minor(self):
        return int(self._minor)
    
    @property
    def micro(self):
        return int(self._micro)  
    
    @property
    def revision(self):
        try:
            return int(self._revision)
        except:
            return self._revision
    
    @property
    def short(self):
        """
        Canonical three-part version number.
        """
        return '{0.major}.{0.minor}.{0.micro}'.format(self)
    
    @property
    def full(self):
        """
        Full version, including the repository revision number.
        """        
        return '{0.major}.{0.minor}.{0.micro}.{0.revision}'.format(self)            
    

"""
Root package, containing all test cases as sub-packages.
"""

import csb.test as test


def updateFiles():
    """
    Refresh the pickled structures in csb/test/data. This might be needed when
    the internal representation of some classes has changed.
    """
    import os
    from csb.io import Pickle
    from csb.bio.io.wwpdb import get
    from csb.bio.structure import Ensemble, ChemElements

    model1 = get('1nz9', model=1)
    model2 = get('1nz9', model=2)
    
    ensemble = Ensemble()
    ensemble.models.append(model1)
    ensemble.models.append(model2)
    Pickle.dump(ensemble, open(os.path.join(test.Config.DATA, '1nz9.full.pickle'), 'wb'))
    
    mse = model1.chains['A'].find(164)
    mse._pdb_name = 'MSE'
    mse.atoms['SD']._element = ChemElements.Se
    mse.atoms['SD']._full_name = 'SE  '
    Pickle.dump(model1, open(os.path.join(test.Config.DATA, '1nz9.model1.pickle'), 'wb'))

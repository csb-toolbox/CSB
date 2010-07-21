def hellinger(p, q):
    """
    Hellingers distance  between probability distributions p and q

    @param p: discrete probability distribution
    @type p: numpy array
    
    @param q: discrete probability distribution
    @type q: numpy array
    """

    from numpy import sqrt, sum

    x = sqrt(p)
    y = sqrt(q)
    hel = ((x-y)**2).sum(1)
    
    return sum(hel/len(hel))

def kl(p, q):
    """
    Kullback-Leiber divergence between probability distributions p and q

    @param p: discrete probability distribution
    @type p: numpy array
    
    @param q: discrete probability distribution
    @type q: numpy array
    """
    
    from csbtbx import log
    from numpy import sum

    lp = log(p)
    lq = log(q)

    kl = sum(p * (lp - lq),1)
    return sum(kl/len(kl))

def symmetrizedkl(p,q):
    """
    Symmetrized variant of the Kullback-Leiber divergence between probability distributions p and q

    @param p: discrete probability distribution
    @type p: numpy array
    
    @param q: discrete probability distribution
    @type q: numpy array
    """
    
    from csbtbx import log
    from numpy import sum

    lp = log(p)
    lq = log(q)

    klp = sum(p * (lp - lq),1)
    klq = sum(q* (lq - lp),1)
    return sum(0.5*klp + 0.5 * klq)/len(p)


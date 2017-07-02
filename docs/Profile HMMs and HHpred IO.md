
## Hidden Markov Model API

The _csb.bio.hmm_ package defines abstractions for working with HHpred's 
HMMs and hit lists. _ProfileHMM_ is the most important object of this 
module. It describes a sequence profile hidden Markov model in the 
way HHpred sees this concept:

* a profile is composed of a list of _HMMLayers_, which contain a 
  number of _States_
* these States can be of different types: _Match_, _Insertion_, 
  _Deletion_, etc.
* a profile contains a multiple alignment, from which it is derived
* this multiple alignment is an A3M (condensed) _Alignment_, where the 
  first sequence is a master sequence
* the match states in all layers correspond to the residues of the master 
  sequence

_ProfileHMM_ objects provide list-like access to their layers:

```python
>>> hmm.layers[1](1)
<HMMLayer>    # first layer: layer at master residue=1
```
 
Every layer provides dictionary-like access to its states:

```python
>>> layer[States.Match](States.Match)
<Match State>
```
 
and every state provides dictionary-like access to its transitions to 
other states:

```python
>>> state = hmm.layers[1](1)[States.match](States.match)
>>> state.transitions[States.Insertion](States.Insertion)
<Transition>       # Match > Insertion
>>> transition.predecessor
<Match State>      # source state
>>> transition.successor
<Insertion State>  # target state
```

Whether this transition points to a state at the same (i) or the next 
layer (i+1) depends on the semantics of the source and the target states.
Building HMMs from scratch is supported through a number of append methods 
at various places:

```python
>>> layer = HMMLayer(...)
>>> layer.append(State(...))
>>> hmm.layers.append(layer)
```

See _HMMLayersCollection_, _HMMLayer_, _EmissionTable_ and _TransitionTable_ 
in our API docs for details.

## HHpred I/O

CSB provides python bindings for working with HHpred's .hhm (HMM) and .hhr 
(HHsearch result) files. These are part of the _csb.bio.io.hhpred_ module, 
which exposes two HHpred format parsers: _HHProfileParser_ and 
_HHOutputParser_. The first one is used to read HMM (.hhm) files:

```python
>>> HHProfileParser("profile.hhm", "profile.pdb").parse()
<ProfileHMM>            # ProfileHMM object
```

while the latter parses HHsearch results files (.hhr):

```python
>>> HHOutputParser().parse_file("hits.hhr"):
<HHpredHitList>        # collection of HHpredHit-s
```

See _ProfileHMM_, _HHpredHitList_ and _HHpredHit_ from _csb.bio.hmm_ 
for details. For text serialization of HMM profiles, see _HHMFileBuilder_ 
from _csb.bio.io.hhpred_ in our API docs.

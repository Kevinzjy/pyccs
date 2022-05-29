# Pyccs

Detecting cyclic consensus sequences in nanopore long reads

```
pip install pyccs
```

## Usage 

```
>>> from pyccs import find_consensus
>>>
>>> segment = "TCCCGGTCATCATAACCCCGATCGTACCCTCTGTCATAATAGTCTCGGCGGCGAGAACTGCCACTGTAAATCTGATCCCTGTCTTGAGCTGCTCTCCATCCACCTCCCTCCACCACCTCCTCCTCTGTATGATCTGCTGTAATAG"
>>>
>>> segment_str, ccs_seq = find_consensus(segment * 3)
>>> segment_str
'10-155;155-300;300-434'
>>>
>>> ccs_seq
'CATAACCCCGATCGTACCCTCTGTCATAATAGTCTCGGCGGCGAGAACTGCCACTGTAAATCTGATCCCTGTCTTGAGCTGCTCTCCATCCACCTCCCTCCACCACCTCCTCCTCTGTATGATCTGCTGTAATAGTCCCGGTCAT'
>>>
>>> len(segment), len(ccs_seq)
(145, 145)
```

## License

The code is released under the `MIT` License. See the LICENSE file for more detail
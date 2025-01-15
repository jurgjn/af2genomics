
import json

def load_json_(file):
    with open(file) as fh:
        return json.load(fh)

def load_seq_(file, chain=None):
    with open(file) as fh:
        lines_ = fh.readlines()
        l_ = None
        r_ = None
        for i, l in enumerate(lines_):
            if l == '    {\n':
                l_ = i
            elif l == '    }\n':
                r_ = i + 1
        lines_ = lines_[l_:r_]
        if not(chain is None):
            assert lines_[2].startswith('        "id": "')
            lines_[2] = '        "id": "%s",\n' % (chain,)
        return lines_

def write_multimer_json_(file, fileB, fileC):
    json_B = load_json_(fileB)
    json_C = load_json_(fileC)

    seq_B = load_seq_(fileB, chain='A')
    seq_C = load_seq_(fileC, chain='B')

    with open(file, 'w') as fh:
        fh.write("""{
  "dialect": "alphafold3",
  "version": 2,
  "name": "%s_%s",
  "sequences": [
%s
%s
  ],
  "modelSeeds": [
    1
  ],
  "bondedAtomPairs": null,
  "userCCD": null
}""" % (json_B['name'], json_C['name'], ''.join(seq_B).rstrip('\n') + ',', ''.join(seq_C).rstrip('\n')))

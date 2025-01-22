
import inspect, os, os.path, urllib.parse

import tqdm.contrib.concurrent

def projectpath(path):
    #dir_ = os.path.dirname(__file__)
    dir_ = '/cluster/project/beltrao/jjaenes'
    return os.path.join(dir_, path)

def workpath(path):
    #dir_ = os.path.dirname(__file__)
    dir_ = '/cluster/work/beltrao/jjaenes'
    return os.path.join(dir_, path)

def localprojectpath(path):
    dir_ = '/Users/jjaenes/euler-home/project'
    return os.path.join(dir_, path)

def localworkpath(path):
    dir_ = '/Users/jjaenes/euler-home/work'
    return os.path.join(dir_, path)

def af2genomicspath(path):
    dir_ = '/cluster/work/beltrao/jjaenes/24.06.10_af2genomics'
    return os.path.join(dir_, path)

def uf(x):
    return '{:,}'.format(x)

def ul(x):
    return uf(len(x))

def printsrc(*args, **kwargs):
    """
        https://stackoverflow.com/questions/3056048/filename-and-line-number-of-python-script
        https://stackoverflow.com/questions/3711184/how-to-use-inspect-to-get-the-callers-info-from-callee-in-python
        https://github.com/snakemake/snakemake/blob/main/snakemake/exceptions.py#L17
    """
    #pprint(dir(inspect.currentframe().f_back))
    #pprint(dir(inspect.getframeinfo(inspect.currentframe().f_back)))
    frameinfo_ = inspect.getframeinfo(inspect.currentframe().f_back)
    #pprint(frameinfo_)
    #pprint(dir(frameinfo_))
    filename = frameinfo_.filename
    lineno = frameinfo_.lineno
    #lineno = workflow.linemaps[filename][ frameinfo_.lineno ]
    print(f'{os.path.basename(filename)}:{lineno}', *args, **kwargs)

def printlen(x, *args, **kwargs):
    name_ = inspect.stack()[1][3] #https://stackoverflow.com/questions/5067604/determine-function-name-from-within-that-function-without-using-traceback
    if name_ != '<module>':
        print(f'{name_}:', uf(len(x)), *args, **kwargs)
    else:
        print(uf(len(x)), *args, **kwargs)

def penc(s):
    """Encode exotic characters within identifiers using percent-encoding, e.g.:
        BF2.649 => BF2%2e649
        UH-232-(+) => UH-232-%28%2b%29
        bis(maltolato)oxovanadium(IV) => bis%28maltolato%29oxovanadium%28IV%29
    """
    #return urllib.parse.quote_plus(s) # Does not decode dots, e.g. BF2.649
    def penc_char(c):
        if c.isalnum() or c in ['_', '-']:
            return c
        else:
            return "%{0:0>2}".format(format(ord(c), "x"))

    return "".join(map(penc_char, s))

def pdec(s):
    """Opposite of pdec(). Sanity check:
    for id_ in ['BF2.649', 'UH-232-(+)', 'bis(maltolato)oxovanadium(IV)']:
        print(id_, '=>', penc(id_))
        assert id_ == pdec(penc(id_))
    """
    # unquote_plus seems to accept characters that aren't encoded by urllib.parse.quote_plus
    return urllib.parse.unquote_plus(s)

def pfile(asset_id=None, compound_id=None, struct_id=None, screen_id=None, step='{prev_steps}', base='{base}', suffix='', v=False):
    """
        Possible alternative to pf() above "asset id"
        For a small, finite set of "entity identifiers" (e.g. model; drug, model, cavity)
        define adhoc manual prefixes (using wildcard restrictions to ease parsing)
        - swiss_model_id:
        - swiss_protein_id:
    """
    if asset_id is None:
        l_asset_id = []
        if not(compound_id is None):
            if compound_id != '{}':
                (compound_ns, compound_name) = compound_id.split('_', maxsplit=1)
                if compound_ns == 'prism' and compound_name[0] != '%':
                    compound_pref = f'{compound_ns}_{compound_name[0]}'
                elif compound_ns == 'prism' and compound_name[0] == '%':
                    compound_pref = f'{compound_ns}_'
                else:
                    compound_pref = compound_ns
            else:
                compound_pref = '{compound_pref}'
                compound_id = '{compound_id}'
            l_asset_id.append(compound_pref)
            l_asset_id.append(compound_id)

        if not(struct_id is None):
            if struct_id != '{}':
                struct_pref = os.path.join(struct_id[:2], struct_id[2:4], struct_id[4:6])
                #struct_pref = struct_id[:2]
            else:
                struct_pref = '{struct_pref}'
                struct_id = '{struct_id}'
            l_asset_id.append(struct_pref)
            l_asset_id.append(struct_id)

        if not(screen_id is None):
            l_asset_id.append(screen_id)

        asset_id = os.path.join(*l_asset_id)

    # add remaining elements: base, step, suffix
    filepath = os.path.join(base, step, f'{asset_id}{suffix}')
    if v: print('pfile: ', filepath)
    return filepath


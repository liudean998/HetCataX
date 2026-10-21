parameters_all_default = {
    'prec': "Medium vasp4 Normal vasp5",
    'nbands': "NELECT/2 + NIONS/2 (non-spinpolarized) = 0.6*NELECT + NMAG (spin-polarized)",
    'nelect': '– (number of valence electrons)',
    'nupdown': None,
    'emin': '– (lowest KS-eigenvalue - ∆)',
    'emax': '– (highest KS-eigenvalue + ∆)',
    'encut': None,
    'nedos': 301,
    'ismear': 1,
    'sigma': 0.2,
    'ispin': 1,
    'istart': 0,
    'icharg': 2,
    'nsw': 0,
    'ibrion': '-1 nsw=0/1 else 0',
    'isif': '0 if ibrion = 0 else 2',
    'isym': '=1 if VASP runs with US-PP’s =2 if PAW data sets are used',
    'lorbit': 'F',  # F 3.2；0 4+
    'voskown': 0,
    'lnoncollinear': 'F',  # 'only for vasp4.5'
    'lsorbit': 'F',  # 'only for vasp4.5'
    'saxis': None,
    'magmom': None,
    'hfscreen': 0,
    'lhfcalc': 'F',
    'aexx': '=0.25 for LHFCALC=.TRUE. =0.0 for LHFCALC=.FALSE.',
    'algo': 'Normal',
    'gga': None,
    'ldipol': 'F',
    'idipol': None,
    'efield': None,
    'ldau': 'F',  # only for vasp4.6
    'ldaul': None,
    'ldauu': None,
    'lmaxmix': 2,
    'ivdw': 0,
    'luse_vdw': 'F',
    'lsol': None,  # 手册中没找到，vasprun中存在
    'lwave': 'T',
    'lcharg': 'T',
    'lvtot': 'F',
    'lelf': 'F',
    'loptics': 'F',
    'solution': 'null',
    'excharge': None,
    'potcar': None,
    'divisions': None,
    'kpointlist': None,
    'weights': None,
    'ediff': 10 ** -4,
    'ediffg': None,  # 10*ediff
    'lreal': ' F',
    'axee': None
}
parameters_all_5 = {
    'prec': "Normal",
    'nbands': None,
    'nelect': None,
    'nupdown': None,
    'emin': None,
    'emax': None,
    'encut': None,
    'nedos': 301,
    'ismear': 1,
    'sigma': 0.2,
    'ispin': 1,
    'istart': 0,
    'icharg': 2,
    'nsw': 0,
    'ibrion': 0,
    'isif': 0,
    'isym': 2,  # '=1 if VASP runs with US-PP’s =2 if PAW data sets are used'
    'lorbit': 'F',
    'voskown': 0,
    'lnoncollinear': 'F',
    'lsorbit': 'F',
    'saxis': None,
    'magmom': None,
    'hfscreen': 0,
    'lhfcalc': 'F',
    'aexx': None,
    'algo': 'Normal',
    'gga': None,
    'ldipol': 'F',
    'idipol': None,
    'efield': None,
    'ldau': 'F',
    'ldaul': None,
    'ldauu': None,
    'lmaxmix': 2,
    'ivdw': 0,
    'luse_vdw': 'F',
    'lsol': None,  # 手册中没找到，vasprun中存在
    'lwave': 'T',
    'lcharg': 'T',
    'lvtot': 'F',
    'lelf': 'F',
    'loptics': 'F',
    'solution': 'null',
    'excharge': None,
    'potcar': None,
    'divisions': None,
    'kpointlist': None,
    'weights': None,
    'ediff': 10 ** -4,
    'ediffg': None,
    'lreal': 'F',
}


def init_vasp5_param(params, pot='PAW'):
    params_5 = parameters_all_5.copy()
    for _i, _p in params.items():
        if _i in list(params_5.keys()):
            params_5[_i] = _p
    if params_5['aexx'] is None:
        if params_5['lhfcalc'] == 'T':
            params_5['aexx'] = 0.25
        else:
            params_5['aexx'] = 0.0
    if params_5['nsw'] in [0, 1]:
        params_5['ibrion'] = -1
        params_5['isif'] = 2
    if pot == 'US':
        params_5['isym'] = 1
    if params_5['ediffg'] is None:
        params_5['ediffg'] = 10 * params_5['ediff']
    return params_5

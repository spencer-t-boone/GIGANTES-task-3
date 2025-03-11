def add_ctr(TT, id, mjd, x0_loc, point = 'Enceladus'):
    ctr = {
        'type': 'control',
        'epoch': str(float(mjd)) + ' TDB',
        'name': 'ctr' + str(id),
        'state': [
            {'name': 'SC_center',
            'point': point,
            'coi': 'Enceladus',
            'axes': 'SEROT',
            'project': False,
            'dynamics': 'SaturnSystem',
            'value': {
                'pos_x': str(x0_loc[0]) + ' km',
                'pos_y': str(x0_loc[1]) + ' km',
                'pos_z': str(x0_loc[2]) + ' km',
                'vel_x': str(x0_loc[3]) + ' km/s',
                'vel_y': str(x0_loc[4]) + ' km/s',
                'vel_z': str(x0_loc[5]) + ' km/s'
            }},
            {'name': 'SC_mass',
             'value': '1000 kg'},
             {'name': 'SC_dv',
              'value': '0 m/s'}
        ]}
    TT.append(ctr)
    return TT

def add_match(TT, id, dt):
    match = {'type': 'match',
             'name': 'match' + str(id),
             'input': 'SC',
             'left': {
                 'reference': 'ctr'+str(id),
                 'dt': str(dt/2) + ' day'},
             'right': {
                 'reference': 'ctr'+str(id+1),
                 'dt': str(dt/2) + ' day'},
             'body': 'Enceladus',
             'vars': 'cart'}
    TT.append(match)
    return TT

def add_man(TT, id):
    man = {'type': 'manoeuvre',
            'name': 'man' + str(id),
            'model': 'impulsive',
            'input': 'SC',
            'thruster': 'main',
            'config': {'point': {'reference': 'ctr'+str(id), 'dt': '0 s'},
                       'direction': {
                           'axes': 'SEROT',
                           'body': 'Enceladus',
                            'dv1': '0.0001 m/s',
                            'dv2': '0.0001 m/s',
                            'ras': '0 deg',
                            'dec': '0 deg'
                       }}}
    TT.append(man)
    return TT
    
def add_end_point(TT, id, dt):
    point = {'type': 'point',
             'input': 'SC',
             'name': 'end_point',
             'point': {'reference': 'ctr' + str(id),
                       'dt': str(dt)+' s'}}
    TT.append(point)
    return TT

def create_trajectory_config():
    config_traj = {'settings': {'relTol': 1e-12, 'steps': 10000000},
               'setup': [
                   {'name': 'SC',
                    'type': 'group',
                    'spacecraft': 'SC', 
                    'input': [
                    {'name': 'center',
                        'type': 'point',
                        'axes': 'ICRF'},
                    {'name': 'mass',
                        'type': 'scalar',
                        'unit': 'kg'},
                    {'name': 'dv',
                        'type': 'scalar',
                        'unit': 'm/s'}
                        ]}]
                    }
    return config_traj

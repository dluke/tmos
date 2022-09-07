
from collections import OrderedDict 
# define the types and formats of a set of header names
# make sure the the header name is smaller than the format length

ldouble = '{:13.8f}'
sdouble = '{:12.9f}'
# Dictionary based data format
# format is {<name>: <type>, <format_str>}
txtform = {
        'time': (float, '{:16.10f}'),
        'process': (str, '{:>10s}'),
        'trigger': (str, '{:>10s}'),
        'x': (float, ldouble),
        'y': (float, ldouble),
        'z': (float, ldouble),
        'ax_x': (float, ldouble),
        'ax_y': (float, ldouble),
        'ax_z': (float, ldouble),
        'rax_x': (float, ldouble),
        'rax_y': (float, ldouble),
        'rax_z': (float, ldouble),
        'trail_x': (float, ldouble),
        'trail_y': (float, ldouble),
        'trail_z': (float, ldouble),
        'center_x': (float, ldouble),
        'center_y': (float, ldouble),
        'center_z': (float, ldouble),
        'd_x': (float, sdouble),
        'd_y': (float, sdouble),
        'd_z': (float, sdouble),
        'd_tx': (float, sdouble),
        'd_ty': (float, sdouble),
        'd_tz': (float, sdouble),
        'd_cx': (float, sdouble),
        'd_cy': (float, sdouble),
        'd_cz': (float, sdouble),
        'd_ax_x': (float, sdouble),
        'd_ax_y': (float, sdouble),
        'd_ax_z': (float, sdouble),
        'd_rax_x': (float, sdouble),
        'd_rax_y': (float, sdouble),
        'd_rax_z': (float, sdouble),
        'dx': (int, '{:>2d}'),
        'dt': (int, '{:>2d}'),
        'dp': (int, '{:>2d}'),
        'pl_avg': (float, '{:11.6f}'),
        'l_total': (float, '{:11.6f}'),
        'npili': (int, '{:>6d}'),
        'nbound': (int, '{:>6d}'),
        'ntaut': (int, '{:>6d}'),
        'fluor_ntfp': (int, '{:>6d}'),
        'iscat_ntfp': (int, '{:>6d}'),
        'nsteps': (int, '{:>6d}'),
        'ncontacts': (int, '{:>10d}'),
        'pidx': (int, '{:06d}'),
        'plength': (float, '{:11.6f}'),
        'inter_excess': (float, '{:13.6f}'),
        'last_shortening': (float, '{:13.6f}'),
        'pleq': (float, '{:11.6f}'),
        'nseg': (int, '{:>4d}'),
        'lasta': (float, '{:11.6f}'),
        '|force|': (float, '{:13.4f}'), # dep this for rms
        'rms': (float, '{:13.4f}'),
        'en_pili': (float, '{:13.4f}'),
        'en_surface': (float, '{:13.4f}'),
        'surfacepz': (float, '{:13.4f}'),
        'moverlap': (float, '{:11.6f}'),

        'pbrf': (float, '{:7.3f}'),
        'statetime': (float, '{:11.6f}'),
        'pcount': (int, '{:>6d}'),
        'npili': (int, '{:>6d}'),
        'anchor_x' : (float, '{:11.6f}'),
        'anchor_y' : (float, '{:11.6f}'),
        'anchor_z' : (float, '{:11.6f}'),
        'paxis_x' : (float, '{:11.6f}'),
        'paxis_y' : (float, '{:11.6f}'),
        'paxis_z' : (float, '{:11.6f}'),
        'attach_x' : (float, '{:11.6f}'),
        'attach_y' : (float, '{:11.6f}'),
        'attach_z' : (float, '{:11.6f}'),
        'isbound' : (int, '{:>7d}'),
        'istaut' : (int, '{:>7d}'),
        'ret_motor' : (int, '{:>10d}'),
        'ext_motor' : (int, '{:>10d}'),
        'cycles' : (int, '{:>7d}')
        }

prettynames = {
        'k_spawn': r'$k_{\mathrm{spawn}}$',
        'pilivar': r'$\kappa$',
        'bound_pili_participation': r'bound pili participation',
        'pole_travel_score': r'Pole Travel Ratio'
}
longnames = {
        'k_ext_off': r'Extension motor unbinding rate ($s^{-1}$)',
        'pilivar': r'Width of Pili Distribution',
        'k_spawn': r'Pili Spawn Rate ($s^{-1}$)',
        'k_resample': r'Pili Resampling Rate ($s^{-1}$)',
        'dwell_time': r'Pili sticking time (s)',
        'anchor_angle_threshold': r'Anchor Angle Threshold (radians)',
        'anchor_angle_smoothing_fraction': r'Anchor Angle Fraction',
}

shortnames = {
        'k_ext_off': r'$k_{\mathrm{ext,off}}$',
        'pilivar': r'$\kappa$',
        'k_spawn': r'$k_\mathrm{spawn}$',
        'k_resample': r'$k_\mathrm{resample}$',
        'dwell_time': r'$\tau_\mathrm{dwell}$',
        'anchor_angle_threshold': r'$t_\mathrm{anchor}$',
        'anchor_angle_smoothing_fraction': r'$a_\mathrm{smooth}$'
}

# define useful constructions from this data 
keys = list(txtform.keys())
htype = {head: txtform[head][0] for head in txtform}
hform = {head: txtform[head][1] for head in txtform}

import re
def format_get_width(formstr):
    return int(re.search(':>?\d+', formstr).group().lstrip(':>'))
def hwform(head, form):
    return '{:>%ds}' % max(format_get_width(form),len(head)) 
hwidth = {head: hwform(head, form) for head, form in list(hform.items())}

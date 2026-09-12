"""fill_manuscript.py

Substitute every {{placeholder}} in manuscript/manuscript.md with a value taken
from results/results.json, and refuse to write anything if a placeholder has no
value behind it.

This is the mechanism that keeps the manuscript honest: a number can only reach
the text by coming out of a run.  If a placeholder is unresolved the script
exits non-zero and names it; if a value is computed and never used it says so,
which catches a claim that was dropped from the text but not from the code.
"""
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
RES = os.path.join(ROOT, 'results', 'results.json')
SRC = os.path.join(ROOT, 'manuscript', 'manuscript.md')
DST = os.path.join(ROOT, 'manuscript', 'manuscript_filled.md')

# The only two values not produced by a run.  Filled once the archive is
# deposited; overridable from the environment so no edit is needed.
GITHUB_URL = os.environ.get(
    'FK_GITHUB_URL', 'https://github.com/dulianzeqiraj/fushe-kuqe-oed')
ZENODO_DOI = os.environ.get(
    'FK_ZENODO_DOI', 'https://doi.org/10.5281/zenodo.PLACEHOLDER')


def rng(vals, fmt='%.3f', sep=' to '):
    return (fmt % min(vals)) + sep + (fmt % max(vals))


def pct(x, fmt='%.1f'):
    return (fmt % (100.0 * x)) + ' per cent'


def sci(x):
    """Render a very small fraction the way it should be read aloud."""
    return '%.1e' % x


def _mean(v):
    v = v if isinstance(v, list) else [v]
    return sum(v) / len(v)


def _median(v):
    v = sorted(v if isinstance(v, list) else [v])
    n = len(v)
    return v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2])



def table_inventory(d, g, pr):
    rows = [
        ('Monitoring points with a DRASTIC index', '%d' % d['n_wells_total']),
        ('Of which inside the digitised boundary', '%d' % d['n_wells_inside']),
        ('With a pumping-test hydraulic conductivity', '%d' % d['n_with_K']),
        ('With a nitrate concentration', '%d' % d['n_with_NO3']),
        ('Inside the boundary, with conductivity', '%d' % d['n_inside_with_K']),
        ('Inside the boundary, with nitrate', '%d' % d['n_inside_with_NO3']),
        ('Aquifer boundary vertices', '%d' % d['n_boundary_vertices']),
        ('Digitised area (km²)', '%.1f' % g['area_polygon_km2']),
        ('Hydraulic conductivity (m d⁻¹)', '%.0f to %.0f, mean %.1f, sd %.1f'
         % (d['K_min'], d['K_max'], d['K_mean'], d['K_sd'])),
        ('Nitrate (mg L⁻¹)', '%.2f to %.2f, mean %.2f, sd %.2f'
         % (d['NO3_min'], d['NO3_max'], d['NO3_mean'], d['NO3_sd'])),
        ('DRASTIC index', '%.0f to %.0f, mean %.1f, sd %.1f'
         % (d['DRASTIC_min'], d['DRASTIC_max'], d['DRASTIC_mean'], d['DRASTIC_sd'])),
        ('Lithological classes of the K wells (coarse/medium/fine)',
         '%d / %d / %d' % tuple(d['facies_counts'])),
        ('ANOVA of K on the lithological label',
         'η² = %.3f, F(2,%d) = %.2f, p = %.2f'
         % (d['facies_K_eta2'], d['n_with_K'] - 3, d['facies_K_F'], d['facies_K_p'])),
        ('Exponential variogram of log₁₀K',
         'nugget %.4f, partial sill %.4f, range %.0f m'
         % (pr['variogram_nugget'], pr['variogram_sill'], pr['variogram_range_m'])),
    ]
    out = ['| Quantity | Value |', '|---|---|']
    out += ['| %s | %s |' % r for r in rows]
    return chr(10).join(out)


def table_loo(p, cal):
    out = ['| Method | Leave-one-out R2 |', '|---|---|']
    for q, pw in enumerate(p['idw_power']):
        out.append('| Inverse distance, power %d | %+.3f |' % (pw, p['idw_R2'][q]))
    for q, rg in enumerate(p['kriging_range_m']):
        out.append('| Ordinary kriging, range %g km | %+.3f |'
                   % (rg / 1000.0, p['kriging_R2'][q]))
    for q, nm in enumerate(p['linear_names']):
        out.append('| Linear regression on %s | %+.3f |' % (nm, p['linear_R2'][q]))
    out.append('| Transient transport model | %+.3f |' % cal['loo_R2'])
    out.append('| Network mean (baseline) | 0.000 |')
    return chr(10).join(out)


def table_design(des, loc, mr):
    rows = [
        ('Transport observations, Fushe-Kuqe',
         '%.1f' % (100 * (1 - (des['sd_goal_after_D_days'] / des['sd_goal_prior_days']) ** 2)),
         '%.1f' % (100 * (1 - (des['sd_goal_after_goal_days'] / des['sd_goal_prior_days']) ** 2)),
         '%.0f' % des['median_separation_m']),
        ('Characterisation points, Fushe-Kuqe',
         '%.1f' % (100 * loc['var_reduction_D']),
         '%.1f' % (100 * loc['var_reduction_goal']),
         '%.0f' % loc['median_separation_m']),
        ('Characterisation points, MRVA',
         '%.1f' % (100 * mr['var_reduction_D']),
         '%.1f' % (100 * mr['var_reduction_goal']),
         '%.0f' % mr['median_separation_m']),
    ]
    out = ['| Observation type and site | D-optimal (%) | Goal-oriented (%) | '
           'Median separation (m) |', '|---|---|---|---|']
    out += ['| %s | %s | %s | %s |' % r for r in rows]
    return chr(10).join(out)


def table_sweep(sw):
    out = ['| Assumption varied | Range | R | Other effect |', '|---|---|---|---|']
    out.append('| Nitrate measurement sd (mg L⁻¹) | %g to %g | %.3f to %.3f | '
               'observation share %.1e to %.1e |'
               % (min(sw['sigma_obs']), max(sw['sigma_obs']),
                  min(sw['R_vs_sigma']), max(sw['R_vs_sigma']),
                  min(sw['frac_obs_vs_sigma']), max(sw['frac_obs_vs_sigma'])))
    out.append('| Pumping-test relative error | %g to %g | %.3f to %.3f | - |'
               % (min(sw['sigma_K_rel']), max(sw['sigma_K_rel']),
                  min(sw['R_vs_sigmaK']), max(sw['R_vs_sigmaK'])))
    out.append('| Longitudinal dispersivity (m) | %g to %g | %.3f to %.3f | '
               'median separation %.0f to %.0f m |'
               % (min(sw['alphaL']), max(sw['alphaL']),
                  min(sw['R_vs_alphaL']), max(sw['R_vs_alphaL']),
                  min(sw['sep_vs_alphaL']), max(sw['sep_vs_alphaL'])))
    out.append('| Vulnerability exponent | %g to %g | - | median separation '
               '%.0f to %.0f m, goal advantage %.2f to %.2f%% |'
               % (min(sw['gammaV']), max(sw['gammaV']),
                  min(sw['sep_vs_gamma']), max(sw['sep_vs_gamma']),
                  100 * min(sw['adv_vs_gamma']), 100 * max(sw['adv_vs_gamma'])))
    return chr(10).join(out)


def build(r):
    d, p, g = r['data'], r['predictability'], r['grid']
    pr, fw, cal = r['prior'], r['forward'], r['calibration']
    fi, des, loc = r['fisher'], r['design'], r['local']
    sam, mr, sw = r['sampling'], r['mrva'], r['sweep']
    m = {}

    # ---------------------------------------------------------------- data
    m['n_total'] = '%d' % d['n_wells_total']
    m['n_inside'] = '%d' % d['n_wells_inside']
    m['n_K'] = '%d' % d['n_with_K']
    m['n_NO3'] = '%d' % d['n_with_NO3']
    m['n_inside_K'] = '%d' % d['n_inside_with_K']
    m['n_inside_NO3'] = '%d' % d['n_inside_with_NO3']
    m['n_boundary'] = '%d' % d['n_boundary_vertices']
    m['K_min'] = '%.0f' % d['K_min']
    m['K_max'] = '%.0f' % d['K_max']
    m['K_mean'] = '%.1f' % d['K_mean']
    m['K_sd'] = '%.1f' % d['K_sd']
    m['NO3_min'] = '%.2f' % d['NO3_min']
    m['NO3_max'] = '%.2f' % d['NO3_max']
    m['NO3_mean'] = '%.2f' % d['NO3_mean']
    m['NO3_sd'] = '%.2f' % d['NO3_sd']
    m['DR_min'] = '%.0f' % d['DRASTIC_min']
    m['DR_max'] = '%.0f' % d['DRASTIC_max']
    m['DR_mean'] = '%.1f' % d['DRASTIC_mean']
    fc = d['facies_counts']
    m['facies_c'], m['facies_m'], m['facies_f'] = ('%d' % fc[0], '%d' % fc[1],
                                                   '%d' % fc[2])
    m['eta2'] = '%.3f' % d['facies_K_eta2']
    m['Fstat'] = '%.2f' % d['facies_K_F']
    m['pval'] = '%.2f' % d['facies_K_p']
    m['df1'] = '2'
    m['df2'] = '%d' % (d['n_with_K'] - 3)

    # ------------------------------------------------------ predictability
    m['best_pred_R2'] = '%+.2f' % p['best_R2']
    m['krig_min'] = '%+.2f' % min(p['kriging_R2'])
    m['krig_max'] = '%+.2f' % max(p['kriging_R2'])
    m['idw_min'] = '%+.2f' % min(p['idw_R2'])
    m['idw_max'] = '%+.2f' % max(p['idw_R2'])
    lr = dict(zip(p['linear_names'], p['linear_R2']))
    m['lin_xy'] = '%+.2f' % lr['X,Y']
    m['lin_drastic'] = '%+.2f' % lr['DRASTIC']
    m['corr_north'] = '%+.2f' % p['corr_northing']
    m['corr_along'] = '%+.2f' % p['corr_alongflow']
    m['corr_drastic'] = '%+.2f' % p['corr_drastic']
    m['corr_K'] = '%+.2f' % p['corr_K']

    # ------------------------------------------------------- grid and prior
    m['nx'] = '%d' % g['nx']
    m['ny'] = '%d' % g['ny']
    m['dx'] = '%g' % g['dx_m']
    m['n_active'] = '%d' % g['n_active']
    m['area_grid'] = '%.1f' % g['area_km2']
    m['area_poly'] = '%.1f' % g['area_polygon_km2']
    m['ne_mean'] = '%.3f' % pr['network_mean']
    m['ne_sd'] = '%.3f' % pr['network_sd']
    m['sd_iv'] = '%.4f' % pr['sd_independent']
    sdf = pr['sd_facies_correlated']
    m['sd_c'], m['sd_m'], m['sd_f'] = ('%.4f' % sdf[0], '%.4f' % sdf[1],
                                       '%.4f' % sdf[2])
    m['sd_used'] = '%.4f' % pr['sd_used']
    m['sd_inflation'] = '%.1f' % pr['sd_inflation_pct']
    m['nugget'] = '%.4f' % pr['variogram_nugget']
    m['sill'] = '%.4f' % pr['variogram_sill']
    m['vrange'] = '%.0f' % pr['variogram_range_m']
    m['n_kl'] = '%d' % pr['n_kl']
    m['kl_var'] = '%.0f' % (100 * pr['kl_variance_captured'])

    # ------------------------------------------------------------- forward
    m['gradient'] = '%.2e' % fw['gradient']
    m['q0'] = '%.4f' % fw['darcy_flux_m_per_day']
    m['gam_min'] = '%.2f' % fw['Gamma_min']
    m['gam_max'] = '%.2f' % fw['Gamma_max']

    # --------------------------------------------------------- calibration
    m['T0_yr'] = '%.1f' % cal['T0_years']
    m['cal_rmse'] = '%.2f' % cal['rmse_mg_per_L']
    m['cal_R2'] = '%.2f' % cal['R2']
    m['cal_loo_rmse'] = '%.2f' % cal['loo_rmse_mg_per_L']
    m['cal_loo_R2'] = '%+.2f' % cal['loo_R2']
    m['drastic_contrib'] = '%.3f' % (cal['R2'] - cal['R2_without_DRASTIC'])
    m['null_R2'] = '%.2f' % cal['null_best_R2']
    m['null_kappa'] = '%.2f' % cal['null_best_kappa']
    m['frac_neg'] = pct(cal['frac_source_negative_before_clipping'], '%.0f')

    # -------------------------------------------------- Fisher information
    m['tr_obs'] = '%.3e' % fi['trace_observations']
    m['tr_prior'] = '%.3e' % fi['trace_prior']
    m['frac_obs_sci'] = sci(fi['observation_share'])
    m['dofs'] = '%.4f' % fi['dofs']
    m['n_modes'] = '%d' % fi['n_modes']
    m['info_gain'] = '%.3f' % fi['information_gain_nats']
    m['sd_red_pct'] = '%.2f' % (100 * fi['sd_reduction_fraction'])
    m['sd_prior_mean'] = '%.5f' % fi['sd_prior_mean']
    m['sd_post_mean'] = '%.5f' % fi['sd_post_mean']
    m['fd_check'] = '%.1e' % fi['finite_difference_check']
    m['fd_order'] = '%.2f' % fi['fd_convergence_order']

    # -------------------------------------------------------------- design
    m['horizon_frac'] = '%.0f' % (100 * des['horizon_fraction'])
    m['n_new'] = '%d' % des['n_new_wells']
    m['sd_local'] = '%.3f' % des['sd_local']
    m['tau_mean_yr'] = '%.1f' % (des['tau_weighted_mean_days'] / 365.25)
    m['sd_goal_prior'] = '%.1f' % des['sd_goal_prior_days']
    m['sd_goal_D'] = '%.1f' % des['sd_goal_after_D_days']
    m['sd_goal_G'] = '%.1f' % des['sd_goal_after_goal_days']
    m['sd_goal_red_pct'] = '%.1f' % (
        100 * (1 - des['sd_goal_after_goal_days'] / des['sd_goal_prior_days']))
    m['median_sep_m'] = '%.0f' % des['median_separation_m']
    m['V_goal'] = '%.1f' % _mean(des['goal_DRASTIC'])
    m['V_D'] = '%.1f' % _mean(des['D_DRASTIC'])
    m['V_domain'] = '%.1f' % des['domain_mean_DRASTIC']
    adv = (des['sd_goal_after_D_days'] ** 2 - des['sd_goal_after_goal_days'] ** 2) \
        / des['sd_goal_after_D_days'] ** 2
    m['goal_adv_pct'] = pct(adv, '%.1f')
    hor = des['information_horizon_m']
    m['horizon_med'] = '%.0f' % _median(hor)
    m['horizon_min'] = '%.0f' % min(hor)
    m['horizon_max'] = '%.0f' % max(hor)
    m['corr_constrained'] = '%.3f' % des['information_correlation_max']
    m['corr_horizon'] = '%.3f' % des['information_correlation_horizon_max']
    m['horizon_over_range'] = '%.1f' % (des['information_horizon_median_m']
                                        / pr['variogram_range_m'])
    m['sd_goal_D_hor'] = '%.1f' % des['sd_goal_after_D_horizon_days']
    m['sd_goal_G_hor'] = '%.1f' % des['sd_goal_after_goal_horizon_days']
    m['corr_free'] = '%.3f' % des['information_correlation_unconstrained_max']
    m['free_span'] = '%.0f' % des['unconstrained_span_m']
    m['min_sep'] = '%.0f' % des['min_separation_m']
    m['R_equiv'] = '%.2f' % des['R_equivalence']
    m['R_inverse'] = '%.0f' % (1.0 / des['R_equivalence'])
    m['sd_pump_ne'] = '%.4f' % des['sd_pump_test_on_ne']
    m['gain_conc'] = '%.3g' % des['gain_conc_best']
    m['gain_pump'] = '%.3g' % des['gain_pump_best']

    # ------------------------------------------------------ local operator
    m['loc_D_pct'] = pct(loc['var_reduction_D'], '%.1f')
    m['loc_goal_pct'] = pct(loc['var_reduction_goal'], '%.1f')
    m['loc_ratio'] = '%.1f' % (loc['var_reduction_goal'] /
                               max(loc['var_reduction_D'], 1e-12))
    m['loc_sep_m'] = '%.0f' % loc['median_separation_m']
    m['loc_sep_km'] = '%.1f' % (loc['median_separation_m'] / 1000.0)

    # ------------------------------------------------------------ sampling
    m['gn_iters'] = '%d' % sam['gn_iterations']
    m['cond_H'] = '%.1f' % sam['condition_number']
    m['iters_nat'] = '%d' % sam['iters_natural']
    m['iters_euc'] = '%d' % sam['iters_euclidean']
    m['sqrt_kappa'] = '%.1f' % sam['sqrt_condition']
    m['M_imp'] = '%d' % sam['M_implicit']
    m['ess_imp'] = '%.3f' % sam['ess_implicit']
    m['M_sir'] = '%d' % sam['M_sir']
    m['ess_sir'] = '%.3f' % sam['ess_sir']

    # ---------------------------------------------------------------- MRVA
    m['mrva_n'] = '%d' % mr['n_wells']
    m['mrva_K_min'] = '%.2f' % mr['K_min']
    m['mrva_K_max'] = '%.1f' % mr['K_max']
    m['mrva_K_mean'] = '%.1f' % mr['K_mean']
    m['mrva_V_min'] = '%.0f' % mr['V_min']
    m['mrva_V_max'] = '%.0f' % mr['V_max']
    m['mrva_mismatch'] = '%d' % mr['drastic_mismatch']
    m['mrva_sep_km'] = '%.1f' % (mr['median_separation_m'] / 1000.0)
    m['mrva_D_pct'] = pct(mr['var_reduction_D'], '%.1f')
    m['mrva_goal_pct'] = pct(mr['var_reduction_goal'], '%.1f')
    m['mrva_V_goal'] = '%.1f' % mr['mean_DRASTIC_goal']
    m['mrva_V_D'] = '%.1f' % mr['mean_DRASTIC_D']

    # --------------------------------------------------------------- sweep
    m['sw_frac_range'] = rng([v for v in sw['frac_obs_vs_sigma']], '%.1e')
    m['sw_R_sigma_range'] = rng(sw['R_vs_sigma'], '%.2f')
    m['sw_R_sigmaK_range'] = rng(sw['R_vs_sigmaK'], '%.2f')
    m['sw_sep_range'] = rng(sw['sep_vs_gamma'], '%.0f')
    m['sw_R_alpha_range'] = rng(sw['R_vs_alphaL'], '%.2f')

    m['TABLE_INVENTORY'] = table_inventory(d, g, pr)
    m['TABLE_LOO'] = table_loo(p, cal)
    m['TABLE_DESIGN'] = table_design(des, loc, mr)
    m['TABLE_SWEEP'] = table_sweep(sw)

    m['GITHUB_URL'] = GITHUB_URL
    m['ZENODO_DOI'] = ZENODO_DOI
    return m


def main():
    with open(RES, encoding='utf-8') as f:
        r = json.load(f)
    m = build(r)
    with open(SRC, encoding='utf-8') as f:
        txt = f.read()

    used, missing = set(), set()

    def sub(mo):
        k = mo.group(1)
        if k in m:
            used.add(k)
            return str(m[k])
        missing.add(k)
        return mo.group(0)

    out = re.sub(r'\{\{(\w+)\}\}', sub, txt)
    if missing:
        print('UNRESOLVED PLACEHOLDERS: ' + ', '.join(sorted(missing)),
              file=sys.stderr)
        return 1
    left = re.findall(r'\{\{.*?\}\}', out)
    if left:
        print('STILL PRESENT: %s' % left, file=sys.stderr)
        return 1
    with open(DST, 'w', encoding='utf-8') as f:
        f.write(out)
    print('wrote %s' % DST)
    print('%d placeholders resolved' % len(used))
    unused = sorted(set(m) - used)
    if unused:
        print('computed but unused: %s' % ', '.join(unused))
    return 0


if __name__ == '__main__':
    sys.exit(main())

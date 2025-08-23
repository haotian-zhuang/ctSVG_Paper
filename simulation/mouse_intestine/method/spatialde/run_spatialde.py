import pandas as pd
import NaiveDE
import SpatialDE

f = 'mouse_intestine'
signal = 'v1'

counts = pd.read_csv('svg/simu/' + f + '/' + signal + '/method/spatialde_global/counts.csv.gz', index_col=0)
sample_info = pd.read_csv('svg/simu/' + f + '/' + signal + '/method/spatialde_global/sample_info.csv.gz', index_col=0)

norm_expr = NaiveDE.stabilize(counts.T).T
resid_expr = NaiveDE.regress_out(sample_info, norm_expr.T, 'np.log(total_counts)').T

all_results = []
for k in sorted(sample_info['cc'].unique(), key=lambda x: int(x)):
    print(f"Processing cluster {k}")
    
    i = sample_info['cc'] == k
    X_clu = sample_info.loc[i, ['row', 'col']]
    counts_clu = counts.loc[i, :]
    gene_clu = (counts_clu > 0).sum(axis=0) >= (0.01 * counts_clu.shape[0])
    resid_expr_clu = resid_expr.loc[i, gene_clu]
    
    results = SpatialDE.run(X_clu, resid_expr_clu)
    results['cluster'] = k
    all_results.append(results)

res_spatialde = pd.concat(all_results, ignore_index=True)
res_spatialde.to_csv('svg/simu/' + f + '/' + signal + '/method/spatialde/raw_res.csv')

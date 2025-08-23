import pandas as pd
import glob
import os
import NaiveDE
import SpatialDE

for filepath in glob.glob('svg/crop/full/spatialde/data/counts/*.csv.gz'):
    f = os.path.basename(filepath)
    
    print(f'Processing file {f}')

    counts = pd.read_csv('svg/crop/full/spatialde/data/counts/' + f, index_col=0)
    sample_info = pd.read_csv('svg/crop/full/spatialde/data/sample_info/' + f, index_col=0)

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
    res_spatialde.to_csv('svg/crop/full/spatialde/raw_res/' + f)

import pandas as pd
import NaiveDE
import SpatialDE

f = 'human_pancreas'
signal = 'v1'

counts = pd.read_csv('svg/simu/' + f + '/' + signal + '/method/spatialde_global/counts.csv.gz', index_col=0)
sample_info = pd.read_csv('svg/simu/' + f + '/' + signal + '/method/spatialde_global/sample_info.csv.gz', index_col=0)

norm_expr = NaiveDE.stabilize(counts.T).T
resid_expr = NaiveDE.regress_out(sample_info, norm_expr.T, 'np.log(total_counts)').T

X = sample_info[['row', 'col']]
results = SpatialDE.run(X, resid_expr)
results.to_csv('svg/simu/' + f + '/' + signal + '/method/spatialde_global/raw_res.csv')

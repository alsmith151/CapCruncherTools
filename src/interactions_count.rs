use std::{error::Error, path::{PathBuf, Path}};

use anyhow::Ok;
use itertools::Itertools;
use log::{debug, info, warn};
use polars::prelude::*;
use pyo3::prelude::*;
use pyo3_polars::PyDataFrame;
use rayon::prelude::*;
use std::collections::HashMap;

fn get_viewpoint(
    parquet: &Path,
    viewpoint: &str,
    remove_exclusions: bool,
    remove_viewpoints: bool,
) -> LazyFrame{
    let mut df_lazy = LazyFrame::scan_parquet(parquet, Default::default())
        .expect("couldnt read parquet")
        .filter(col("viewpoint").eq(lit(viewpoint)));

    if remove_exclusions == true {
        df_lazy = df_lazy.filter(col("exclusion").ne(&col("viewpoint")).into());
    }

    if remove_viewpoints == true {
        df_lazy = df_lazy.filter(col("capture_count").eq(0));
    }

    df_lazy = df_lazy.select(&[col("parent_id"), col("restriction_fragment")]);

    df_lazy
}

#[pyclass]
pub struct RestrictionFragmentCounter {
    parquet: PathBuf,
    viewpoint: String,
    remove_exclusions: bool,
    remove_viewpoints: bool,
    subsample: bool,
    subsample_frac: f64,
}

#[pymethods]
impl RestrictionFragmentCounter {
    #[new]
    #[pyo3(text_signature = "(parquet, viewpoint, remove_exclusions, remove_viewpoints, subsample, subsample_frac)")]
    pub fn new(
        parquet: PathBuf,
        viewpoint: String,
        remove_exclusions: Option<bool>,
        remove_viewpoints: Option<bool>,
        subsample: Option<bool>,
        subsample_frac: Option<f64>,
    ) -> Self {
        let remove_exclusions = match remove_exclusions {
            Some(b) => b,
            None => false,
        };

        let remove_viewpoints = match remove_viewpoints {
            Some(b) => b,
            None => false,
        };

        let subsample = match subsample {
            Some(b) => b,
            None => false,
        };

        let subsample_frac = match subsample_frac {
            Some(f) => f,
            None => 1.0,
        };

        Self {
            parquet,
            viewpoint,
            remove_exclusions,
            remove_viewpoints,
            subsample,
            subsample_frac,
        }
    }

    fn count_interactions(&self) -> PyDataFrame {
        // Group the dataframe by parent_id
        let df = 
            get_viewpoint(&self.parquet, &self.viewpoint, self.remove_exclusions, self.remove_viewpoints)
            .collect()
            .expect("couldnt collect lazyframe");

        let interaction_counts: HashMap<(u64, u64), i32> = df
            .partition_by(vec!["parent_id"])
            .expect("couldnt partition by parent_id")
            .into_iter()
            .map(|df| {
                let ser = &df
                    .select_series(["restriction_fragment"])
                    .expect("couldnt extract restriction_fragment column")[0];
                let rf = ser.u64().expect("should be u64");

                let mut rf_combs = HashMap::new();
                for comb in rf.into_iter().combinations(2) {
                    let (a, b) = (comb[0], comb[1]);

                    match (a, b) {
                        (Some(a), Some(b)) => {
                            let (a, b) = vec![a, b].into_iter().sorted().collect_tuple().unwrap();
                            *rf_combs.entry((a, b)).or_insert(0) += 1;
                        }
                        _ => {}
                    }
                }

                rf_combs
            })
            .fold(
                HashMap::new(),
                |mut a: HashMap<(u64, u64), i32>, b: HashMap<(u64, u64), i32>| {
                    for (k, v) in b {
                        *a.entry(k).or_insert(1) += v;
                    }
                    a
                },
            );

        let df_counts = DataFrame::new(vec![
            Series::new(
                "bin1_id",
                interaction_counts
                    .keys()
                    .map(|(a, _)| *a)
                    .collect::<Vec<_>>(),
            ),
            Series::new(
                "bin2_id",
                interaction_counts
                    .keys()
                    .map(|(_, b)| *b)
                    .collect::<Vec<_>>(),
            ),
            Series::new(
                "count",
                interaction_counts.values().map(|v| *v).collect::<Vec<_>>(),
            ),
        ])
        .expect("couldnt create dataframe");

        PyDataFrame(df_counts)
    }
}
